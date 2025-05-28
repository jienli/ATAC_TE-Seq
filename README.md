# ATAC_TE-Seq

## Environment Setup

This project uses two Conda environments: `ATAC_Core` and `snakemake`. You can recreate them using the provided files in the `envs/` directory.

### 1. ATAC_Core Environment

#### Using Conda
```bash
conda env create -f envs/ATAC_Core_Environment.yml
conda activate ATAC_Core
```

#### Using pip (if needed)
```bash
pip install -r envs/ATAC_Core_Environment.txt
```

### 2. Snakemake Environment

#### Using Conda
```bash
conda env create -f envs/snakemakeEnvironment.yml
conda activate snakemake
```

#### Using pip (if needed)
```bash
pip install -r envs/snakemakeEnvironment.txt
```
