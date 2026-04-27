#!/bin/bash
#SBATCH --job-name=manual_eggnog_db
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
EGGNOG_DATA_DIR="/uoa/scratch/users/r02hw22/project_dark_genes/eggnog_data"
BASE="http://eggnog5.embl.de/download/emapperdb-5.0.2"

cd "$PROJECT_DIR"
mkdir -p "$EGGNOG_DATA_DIR"
cd "$EGGNOG_DATA_DIR"

echo "Started: $(date)"
echo "Working directory: $(pwd)"
echo "Database source: $BASE"

df -h "$EGGNOG_DATA_DIR"
df -i "$EGGNOG_DATA_DIR"

echo "Removing failed/partial previous downloads..."
rm -f eggnog.db eggnog.db.gz
rm -f eggnog_proteins.dmnd eggnog_proteins.dmnd.gz
rm -f eggnog.taxa.db eggnog.taxa.tar.gz

echo "Checking remote files exist..."
wget --spider "${BASE}/eggnog.db.gz"
wget --spider "${BASE}/eggnog.taxa.tar.gz"
wget --spider "${BASE}/eggnog_proteins.dmnd.gz"

echo "Downloading eggnog.db.gz..."
wget --progress=dot:giga --tries=10 --timeout=60 --waitretry=30 \
  -O eggnog.db.gz "${BASE}/eggnog.db.gz"

test -s eggnog.db.gz
ls -lh eggnog.db.gz

echo "Decompressing eggnog.db.gz..."
gunzip -f eggnog.db.gz
test -s eggnog.db
ls -lh eggnog.db

echo "Downloading eggnog.taxa.tar.gz..."
wget --progress=dot:giga --tries=10 --timeout=60 --waitretry=30 \
  -O eggnog.taxa.tar.gz "${BASE}/eggnog.taxa.tar.gz"

test -s eggnog.taxa.tar.gz
ls -lh eggnog.taxa.tar.gz

echo "Extracting eggnog.taxa.tar.gz..."
tar -zxf eggnog.taxa.tar.gz
rm eggnog.taxa.tar.gz
test -s eggnog.taxa.db
ls -lh eggnog.taxa.db

echo "Downloading eggnog_proteins.dmnd.gz..."
wget --progress=dot:giga --tries=10 --timeout=60 --waitretry=30 \
  -O eggnog_proteins.dmnd.gz "${BASE}/eggnog_proteins.dmnd.gz"

test -s eggnog_proteins.dmnd.gz
ls -lh eggnog_proteins.dmnd.gz

echo "Decompressing eggnog_proteins.dmnd.gz..."
gunzip -f eggnog_proteins.dmnd.gz
test -s eggnog_proteins.dmnd
ls -lh eggnog_proteins.dmnd

echo "Final database directory:"
ls -lh "$EGGNOG_DATA_DIR"
du -sh "$EGGNOG_DATA_DIR"

echo "Verifying eggNOG-mapper database recognition..."

export EGGNOG_DATA_DIR="$EGGNOG_DATA_DIR"

EMAPPER="/uoa/scratch/users/r02hw22/project_dark_genes/conda_envs/eggnog-mapper-2.1.13/bin/emapper.py"

if [[ ! -x "$EMAPPER" ]]; then
    echo "ERROR: emapper.py not found or not executable: $EMAPPER" >&2
    exit 1
fi

"$EMAPPER" --version