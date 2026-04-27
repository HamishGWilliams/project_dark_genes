#!/bin/bash
#SBATCH --mem 48G
#SBATCH --partition uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

module load miniconda3/latest

# install
conda create -n eggnog-mapper -c conda-forge -c bioconda eggnog-mapper -y
conda activate eggnog-mapper

# check
which emapper.py
which download_eggnog_data.py

emapper.py --version
download_eggnog_data.py --help | head

# download databases
EGGNOG_DATA_DIR="/uoa/scratch/users/r02hw22/project_dark_genes/eggnog_data"
mkdir -p "$EGGNOG_DATA_DIR"

download_eggnog_data.py \
  --data_dir "$EGGNOG_DATA_DIR" \
  -y

