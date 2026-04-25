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

# Archive results

mkdir -p 04_reports/busco

BUSCO_DIR="/uoa/scratch/users/r02hw22/project_dark_genes/equina_busco_proteins"

# Inspect what BUSCO produced
find "$BUSCO_DIR" -maxdepth 3 -type f | sort

# Copy summary files if present
cp "$BUSCO_DIR"/short_summary*.txt  04_reports/busco/ 2>/dev/null || true
cp "$BUSCO_DIR"/short_summary*.json 04_reports/busco/ 2>/dev/null || true

# Locate and copy full_table.tsv wherever BUSCO placed it
FULL_TABLE=$(find "$BUSCO_DIR" -name "full_table.tsv" | head -n 1)

if [ -n "$FULL_TABLE" ]; then
    cp "$FULL_TABLE" 04_reports/busco/
else
    echo "full_table.tsv not found under $BUSCO_DIR"
fi

# Make a lightweight plain-text snapshot of the key metrics
SUMMARY_TXT=$(find "$BUSCO_DIR" -name "short_summary*.txt" | head -n 1)

if [ -n "$SUMMARY_TXT" ]; then
    grep -E 'C:|Complete BUSCOs|Single-copy|Duplicated|Fragmented|Missing|Total BUSCO' \
        "$SUMMARY_TXT" > 04_reports/busco/busco_metrics_snapshot.txt
fi

# Check what was archived
ls -lh 04_reports/busco