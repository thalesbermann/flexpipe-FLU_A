# flexpipe-FLU_A
Nextstrain pipeline for genomic epidemiology of Influenza A virus (H1N1 or H3N2), segments HA and NA. This tool is derived from [flexpipe](https://github.com/InstitutoTodosPelaSaude/flexpipe) and adapted for Influenza A, similarly to how flexpipe-RSV was created for respiratory syncytial virus.

This repository contains all essential files to generate Influenza A Nextstrain builds for the HA and NA segments of H1N1 and H3N2. Using this pipeline, users can perform genomic epidemiology analyses, visualize phylogeographic results, and track Influenza A spread based on genomic data and associated metadata.

![Nextstrain panel with RSV overview](nextstrain_results.png)

### Getting started
To run this pipeline for Influenza A projects, see the instructions available in the original [flexpipe repository](https://github.com/InstitutoTodosPelaSaude/flexpipe), which covers Unix CLI navigation, installation of a Nextstrain environment with conda/mamba, and a step-by-step tutorial on generating a Nextstrain build (preparing, aligning, and visualizing genomic data).

---

### Influenza A subsampling (NCBI Virus metadata)
This repository includes a custom QC and subsampling pipeline designed for Influenza A sequence datasets derived from NCBI Virus.

### Usage
Place your sequences and metadata files in the working directory, then run:
```bash
python3 subsample_FLU_A.py
```

### Dataset retrieval (NCBI Virus)
Sequences were retrieved from NCBI Virus using the following filters:
- Virus: Influenza A
- Subtype: H1N1 or H3N2
- Segment: HA or NA
- Minimum length: ≥ 1000 bp
- Collection date: ≥ 2010

### Quality control (viralQC)
Sequences were processed using [viralQC](https://github.com/InstitutoTodosPelaSaude/viralQC):
```bash
vqc run \
  --input sequences.fasta \
  --output-dir qc_output \
  --output-file results.tsv \
  --datasets-dir viralqc_datasets \
  --blast-database viralqc_datasets/blast.fasta \
  --blast-database-metadata viralqc_datasets/blast.tsv \
  --cores 8 \
  -v
```
Sequences were retained only if they met:
- `genomeQuality` = A or B
- Coverage ≥ 70%

### Filtering and preprocessing
- Retains only complete dates (`YYYY-MM-DD`)
- Includes only sequences collected from **year 2015 onwards**
- Metadata is merged with viralQC results and matched by accession
- A `clade` column is added for downstream diversity-aware subsampling

### Subsampling strategy
The dataset is reduced to approximately 2,000 sequences while preserving geographic, temporal, and evolutionary diversity.

**Brazil (focal country)**
- No subsampling applied — all sequences retained

**Global dataset**
- Sampling unit: country + year
- Maximum of **15 sequences per country per year**
- Clade diversity: stratified by `clade_major` (first hierarchical level of clade); at least one sequence per `clade_major` is preserved when possible, with remaining slots filled randomly
- Maximum of **100 sequences per country** (excluding focal country); if exceeded, the most recent sequences are kept based on `Collection_Date`

### Outputs
The script generates:
- `metadata_H1N1_HA_subsampled.tsv` → filtered and subsampled metadata #example if H1N1_HA

### Sequence extraction
You can extract the corresponding sequences using:
```bash
cut -f1 metadata_H1N1_HA_subsampled.tsv | tail -n +2 > keep_ids.txt
seqkit grep -f keep_ids.txt sequences.fasta -o subsampled.fasta
```

---

### Adjustments for Influenza A runs
The Snakefile provided here is pre-configured for Influenza A - H1N1 and H3N2 - HA and NA segments. Since these are individual gene segments (not full genomes), no UTR trimming is applied — masking is set to 1 bp on both ends as a minimal placeholder. The clock rate is left undefined so that Nextstrain (augur) estimates it automatically from the data. Example:

```python
rule parameters:
    params:
        mask_5prime = 1,
        mask_3prime = 1,
        bootstrap = 1000,
        model = "MFP",
        root = "least-squares",
```

---


### Author
* **Thales Bermann, Instituto Todos pela Saúde (ITpS)** - thalesbermann@gmail.com

## License
This project is licensed under the MIT License.
