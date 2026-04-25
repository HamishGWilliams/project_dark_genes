#!/bin/bash
#SBATCH --mem 128G
#SBATCH --partition uoa-compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

cd /uoa/home/r02hw22/sharedscratch/project_dark_genes/

# copy trimmed reads
cp -v /uoa/scratch/users/r02hw22/MS_RNAseq/trimmed/*.fastq.gz 05_rnaseq/00_raw_trimmed/

# copy fastq files
cp -rv /uoa/scratch/users/r02hw22/MS_RNAseq/trimmed/fastqc_out/* 05_rnaseq/01_qc/fastqc/

# copy fastqc outputs
cp -rv /uoa/scratch/users/r02hw22/MS_RNAseq/trimmed/fastqc_out/* 05_rnaseq/01_qc/fastqc/

# copy multiqc outputs
cp -rv /uoa/scratch/users/r02hw22/MS_RNAseq/trimmed/multiqc_out/* 05_rnaseq/01_qc/multiqc/