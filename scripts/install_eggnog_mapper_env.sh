#!/bin/bash
#SBATCH --job-name=install_eggnog
#SBATCH --mem=48G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

mkdir -p logs/outputs logs/errors

module load miniconda3/latest

source "$(conda info --base)/etc/profile.d/conda.sh"

export CONDA_ENVS_DIRS="/uoa/scratch/users/r02hw22/project_dark_genes/conda_envs"
export CONDA_PKGS_DIRS="/uoa/scratch/users/r02hw22/project_dark_genes/conda_pkgs"

mkdir -p "$CONDA_ENVS_DIRS" "$CONDA_PKGS_DIRS"

ENV_PREFIX="${CONDA_ENVS_DIRS}/eggnog-mapper-2.1.13"

echo "Conda base: $(conda info --base)"
echo "Environment prefix: $ENV_PREFIX"
echo "Package cache: $CONDA_PKGS_DIRS"

conda config --show channels
conda config --set channel_priority strict

# Clean possibly corrupted downloads from previous failed attempt.
conda clean -a -y

# Remove partial environment if present.
rm -rf "$ENV_PREFIX"

# Use Python 3.10 because current Bioconda eggNOG-mapper requires Python >=3.7,<3.12.
conda create \
  -p "$ENV_PREFIX" \
  -c conda-forge \
  -c bioconda \
  python=3.10 \
  eggnog-mapper=2.1.13 \
  -y

conda activate "$ENV_PREFIX"

echo "Checking installed tools"
which emapper.py
which download_eggnog_data.py
which diamond || true
which hmmscan || true
which mmseqs || true

emapper.py --version || true
download_eggnog_data.py --help | head

echo "eggNOG-mapper environment installation completed successfully."