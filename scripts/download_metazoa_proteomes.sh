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

# downloaded the proteome lists name from uniprot
# metazoa proteomes only (fallback)

# capture UPIDs to download TrEMBL proteomes
awk 'NR>1 {print $1}' \
  02_annotation/reference_dbs/metazoa_proteomes.tsv \
  > 02_annotation/reference_dbs/metazoa_selected_upids.txt

# Download
OUT=02_annotation/reference_dbs/metazoa_selected.fasta
: > "$OUT"

while read -r UPID; do
  echo "Downloading $UPID"
  curl -L -G 'https://rest.uniprot.org/uniprotkb/stream' \
    --data-urlencode 'format=fasta' \
    --data-urlencode "query=(proteome:${UPID}) AND (reviewed:false)" \
    >> "$OUT"
  echo >> "$OUT"
done < 02_annotation/reference_dbs/metazoa_selected_upids.txt



