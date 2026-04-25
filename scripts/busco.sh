#!/bin/bash
#SBATCH --mem 96G
#SBATCH --partition uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

cd /uoa/home/r02hw22/sharedscratch/project_dark_genes/

module load busco/5.3.2

busco -i 00_raw/equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa \
      -m proteins \
      -l metazoa_odb10 \
      -o equina_busco_proteins \
      -c 32

mkdir -p 04_reports/busco

find equina_busco_proteins -maxdepth 2 -type f | sort

cp equina_busco_proteins/short_summary*.txt 04_reports/busco/ 2>/dev/null || true
cp equina_busco_proteins/short_summary*.json 04_reports/busco/ 2>/dev/null || true
cp equina_busco_proteins/run_metazoa_odb10/full_table.tsv 04_reports/busco/ 2>/dev/null || true
