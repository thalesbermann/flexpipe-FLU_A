#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# ======================
# CONFIGURAÇÕES
# ======================

INPUT = "sequences_h1n1_HA_QC.tsv"
OUTPUT = "metadata_H1N1_HA_subsampled.tsv"

START_YEAR = 2015
MAX_PER_COUNTRY_YEAR = 15
MAX_PER_COUNTRY_GLOBAL = 100
SEED = 42

np.random.seed(SEED)

# ======================
# CHECAGEM DE ARQUIVO
# ======================

input_path = Path(INPUT)
if not input_path.exists():
    sys.exit(f"Erro: arquivo não encontrado -> {INPUT}")

# ======================
# LEITURA
# ======================

df = pd.read_csv(INPUT, sep="\t", dtype=str)
print(f"Dataset inicial: {len(df):,}")

# ======================
# CHECAGEM DE COLUNAS
# ======================

required_cols = ["Accession", "Collection_Date", "Country", "clade"]
missing_cols = [col for col in required_cols if col not in df.columns]

if missing_cols:
    sys.exit(f"Erro: colunas ausentes no input -> {missing_cols}")

# ======================
# FILTROS
# ======================

# manter apenas datas completas YYYY-MM-DD
df = df[df["Collection_Date"].str.match(r"^\d{4}-\d{2}-\d{2}$", na=False)].copy()
print(f"Após filtro de data completa: {len(df):,}")

# criar ano
df["year"] = pd.to_numeric(df["Collection_Date"].str[:4], errors="coerce")
df = df[df["year"].notna()].copy()
df["year"] = df["year"].astype(int)

# filtrar ano inicial
df = df[df["year"] >= START_YEAR].copy()
print(f"Após filtro >= {START_YEAR}: {len(df):,}")

# remover sem país ou sem clado
df = df.dropna(subset=["Country", "clade"]).copy()
df = df[(df["Country"].str.strip() != "") & (df["clade"].str.strip() != "")].copy()
print(f"Após remover vazios em Country/clade: {len(df):,}")

# clado principal
df["clade_major"] = df["clade"].str.split(".").str[0]

# marcar Brasil
df["is_brazil"] = df["Country"].str.contains("Brazil", case=False, na=False)

df_br = df[df["is_brazil"]].copy()
df_global = df[~df["is_brazil"]].copy()

print(f"Brasil mantido integralmente: {len(df_br):,}")
print(f"Global antes do subsampling: {len(df_global):,}")

# guardar colunas originais para manter consistência
original_cols = df_global.columns.tolist()

# ======================
# FUNÇÃO DE SUBSAMPLING GLOBAL
# ======================

def sample_country_year(group: pd.DataFrame) -> pd.DataFrame:
    """
    Subsample por Country + year, preservando diversidade de clade_major.
    Regras:
    - Se grupo <= MAX_PER_COUNTRY_YEAR: mantém tudo
    - Se grupo > MAX_PER_COUNTRY_YEAR:
        1. tenta garantir 1 sequência por clade_major
        2. completa o restante aleatoriamente
    """
    country, year = group.name
    group = group.copy()

    n_total = len(group)

    if n_total <= MAX_PER_COUNTRY_YEAR:
        group["Country"] = country
        group["year"] = year
        return group

    sampled_parts = []

    # garante 1 por clade_major
    for _, sub in group.groupby("clade_major"):
        sampled_parts.append(sub.sample(1, random_state=SEED))

    sampled = pd.concat(sampled_parts).drop_duplicates()

    # se já passou do limite, corta
    if len(sampled) >= MAX_PER_COUNTRY_YEAR:
        sampled = sampled.sample(MAX_PER_COUNTRY_YEAR, random_state=SEED)
        sampled["Country"] = country
        sampled["year"] = year
        return sampled

    # completa com o restante
    remaining = group.loc[~group.index.isin(sampled.index)]
    n_needed = MAX_PER_COUNTRY_YEAR - len(sampled)

    if len(remaining) > 0 and n_needed > 0:
        extra = remaining.sample(min(n_needed, len(remaining)), random_state=SEED)
        sampled = pd.concat([sampled, extra]).drop_duplicates()

    sampled["Country"] = country
    sampled["year"] = year
    return sampled

# ======================
# SUBSAMPLING GLOBAL
# ======================

df_global_sampled = (
    df_global
    .drop(columns=["Country", "year"])
    .groupby([df_global["Country"], df_global["year"]], group_keys=False)
    .apply(sample_country_year, include_groups=False)
    .reset_index(drop=True)
)

# garantir ordem de colunas
for col in original_cols:
    if col not in df_global_sampled.columns:
        df_global_sampled[col] = pd.NA
df_global_sampled = df_global_sampled[original_cols]

print(f"Global após subsampling por Country+year: {len(df_global_sampled):,}")

# ======================
# CAP GLOBAL POR PAÍS
# ======================

def cap_country_global(group: pd.DataFrame) -> pd.DataFrame:
    """
    Aplica cap máximo por país apenas no conjunto global.
    Se passar do limite, mantém as mais recentes.
    """
    country = group.name
    group = group.copy()

    if len(group) <= MAX_PER_COUNTRY_GLOBAL:
        group["Country"] = country
        return group

    group = group.sort_values("Collection_Date", ascending=False).head(MAX_PER_COUNTRY_GLOBAL)
    group["Country"] = country
    return group

df_global_capped = (
    df_global_sampled
    .drop(columns=["Country"])
    .groupby(df_global_sampled["Country"], group_keys=False)
    .apply(cap_country_global, include_groups=False)
    .reset_index(drop=True)
)

# garantir ordem de colunas
for col in original_cols:
    if col not in df_global_capped.columns:
        df_global_capped[col] = pd.NA
df_global_capped = df_global_capped[original_cols]

print(f"Global após cap por país ({MAX_PER_COUNTRY_GLOBAL}): {len(df_global_capped):,}")

# ======================
# JUNTAR BRASIL + GLOBAL
# ======================

df_final = pd.concat([df_br, df_global_capped], ignore_index=True)

# remover duplicatas de acesso se houver
df_final = df_final.drop_duplicates(subset=["Accession"]).copy()

print(f"Total final: {len(df_final):,}")

# ======================
# RELATÓRIOS RÁPIDOS
# ======================

print("\nTop 20 países no dataset final:")
print(df_final["Country"].value_counts().head(20).to_string())

print("\nTop 20 clades no dataset final:")
print(df_final["clade"].value_counts().head(20).to_string())

# ======================
# SALVAR
# ======================

df_final.to_csv(OUTPUT, sep="\t", index=False)
print(f"\nArquivo salvo em: {OUTPUT}")
