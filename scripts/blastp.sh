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

module load blast-plus/2.16.0

QUERY="00_raw/equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa"
REF_FASTA="/path/to/uniprot_reference_proteins.fasta"
BLAST_DB="02_annotation/db/uniprot_ref_blast"
OUT_TSV="02_annotation/blastp/equina_vs_uniprot_ref.blastp.tsv"

makeblastdb \
  -in "$REF_FASTA" \
  -dbtype prot \
  -out "$BLAST_DB"

blastp \
  -query "$QUERY" \
  -db "$BLAST_DB" \
  -out "$OUT_TSV" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -num_threads 32

awk -F '\t' '!seen[$1]++' "$OUT_TSV" > 02_annotation/blastp/equina_vs_uniprot_ref.top_hits.tsv

wc -l "$OUT_TSV" 02_annotation/blastp/equina_vs_uniprot_ref.top_hits.tsv
head 02_annotation/blastp/equina_vs_uniprot_ref.top_hits.tsv