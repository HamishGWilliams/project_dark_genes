#!/bin/bash
#SBATCH --mem 128G
#SBATCH --partition uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

cd /uoa/home/r02hw22/sharedscratch/project_dark_genes/

module load diamond/2.1.10
# module load blast-plus/2.16.0

QUERY="00_raw/equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa"
REF_FASTA="/path/to/uniprot_reference_proteins.fasta"
DB_PREFIX="02_annotation/db/uniprot_ref"
OUT_TSV="02_annotation/diamond/equina_vs_uniprot_ref.diamond.tsv"

diamond makedb \
  --in "$REF_FASTA" \
  --db "$DB_PREFIX"

diamond blastp \
  --query "$QUERY" \
  --db "$DB_PREFIX" \
  --out "$OUT_TSV" \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
  --threads 32

awk -F '\t' '!seen[$1]++' "$OUT_TSV" > 02_annotation/diamond/equina_vs_uniprot_ref.top_hits.tsv

wc -l "$OUT_TSV" 02_annotation/diamond/equina_vs_uniprot_ref.top_hits.tsv
head 02_annotation/diamond/equina_vs_uniprot_ref.top_hits.tsv
