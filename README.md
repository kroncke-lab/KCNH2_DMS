# Deep Mutational Scan of KCNH2

## Overview
This repository contains the analysis code supporting the publication:

> Ng et al., *The American Journal of Human Genetics* (2022). 109(7):1208–1216. doi: 10.1016/j.ajhg.2022.05.003; PMID: 35688148  PMCID: PMC9300756

The scripts implement a deep mutational scanning (DMS) workflow for KCNH2 missense variants. The pipeline:

1. Builds barcode consensus sequences from raw reads.
2. Identifies protein-coding variants by comparing reads to a supplied wild-type reference.
3. Links each barcode to its associated variant to produce a barcode–variant index.
4. Quantifies barcode abundance in sort-seq experiments to measure enrichment across cell populations.

All code and data referenced by the manuscript are included in this repository.

## Repository layout

- `main.py` – entry point that orchestrates barcode identification and enrichment quantification.
- `Contents/` – supplementary notebooks, figures, and tables associated with the manuscript.
- `R_src/` – R utilities used for downstream analysis and visualization.
- `models/` – trained models used for DMS interpretation.

## Requirements

- Python 3.7 or later
- [`virtualenv`](https://virtualenv.pypa.io/)
- `matplotlib==3.4.2`
- `biopython==1.77` *(required for the `Alphabet.IUPAC` module)*
- `pandas`
- `numpy`

> **Note:** The Biopython dependency is pinned to 1.77. Later versions remove `Alphabet.IUPAC`, which is required by the current scripts.

## Recommended environment setup

```bash
# create a working directory and clone the repository
cd /path/to/working/directory

# (optional) place input FASTQ/FASTA files in a dedicated folder
mkdir DNA

# create and activate a virtual environment
python3 -m venv KCNH2_DMS-env
source KCNH2_DMS-env/bin/activate

# install dependencies
pip install matplotlib==3.4.2 biopython==1.77 pandas numpy
```

If you are working on a managed cluster that provides environment modules, the same steps can be adapted to use `module load python3` prior to creating the virtual environment.

## Input data expectations

Create a directory (e.g., `/path/to/working/directory/DNA/`) containing only the FASTQ files generated from the sort-seq experiment and a reference FASTA named `plasmid.fasta`. The FASTA must have **two** entries:

```text
>dna_of_full_plasmid
ATCGATCG...
>dna_of_protein_coding_region
ATCGATCG...
```

The barcode region must have identical flanking sequences across all reads because the consensus step performs stringent filtering on these positions.

## Configuring and running the pipeline

Edit the `folder` variable at the top of `main.py` to point to the directory containing your FASTQ and FASTA files:

```python
folder = "/path/to/working/directory/DNA/"
```

Then execute the pipeline from the virtual environment:

```bash
python main.py
```

If you only need to build the barcode–variant index and do not want to run the enrichment quantification, comment out the commands after line 15 in `main.py`.

### Adjusting the barcode length

The default barcode length is 18 base pairs. To work with barcodes of a different size, update the constants in `consensus.py` (lines 65, 81, 105, and 158) to reflect the desired length.

## Tips for new Python users

Directory structure matters when working with virtual environments. One robust layout is:

```text
/working/directory/
├── DNA/                    # FASTQ and plasmid.fasta files
├── KCNH2_DMS/              # cloned repository
└── KCNH2_DMS-env/          # virtual environment
```

Activate the environment from `/working/directory/` and run the scripts via:

```bash
source KCNH2_DMS-env/bin/activate
python KCNH2_DMS/main.py
```

This layout ensures all relative paths in the repository resolve correctly.

## Citation

If you use this code or data in your work, please cite the original manuscript referenced above.
