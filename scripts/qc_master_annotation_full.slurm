#!/bin/bash
#SBATCH --job-name=qc_master_full
#SBATCH --mem=64G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=06:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

mkdir -p logs/outputs logs/errors 02_annotation/master/full/qc

module load python/3.11.9 || module load python/3.9.12 || true

echo "========================================"
echo "Running full master annotation QC"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Project: $PROJECT_DIR"
echo "========================================"

bash scripts/qc_master_annotation_full.sh

echo
echo "Finished: $(date)"
echo "Full master annotation QC completed."