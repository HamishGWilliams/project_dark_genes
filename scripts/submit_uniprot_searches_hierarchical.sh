#!/bin/bash

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

mkdir -p logs/outputs logs/errors
mkdir -p 02_annotation/db
mkdir -p 02_annotation/diamond
mkdir -p 02_annotation/blastp
mkdir -p 02_annotation/intermediate
mkdir -p scripts

QUERY="00_raw/equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa"

# Hierarchical strategy:
# 1. Run DIAMOND for all three UniProt databases.
# 2. After Swiss-Prot DIAMOND completes, extract proteins with no Swiss-Prot DIAMOND hit.
# 3. Run BLASTp only for the no-hit subset against the selected databases.
#
# Toggle BLASTp fallback databases below. Usually keep Swiss-Prot off because
# BLASTp is intended here as a secondary rescue step after DIAMOND Swiss-Prot.
RUN_DIAMOND="true"
RUN_BLASTP_NOHIT="true"
BLASTP_FALLBACK_DBS=(
  "trembl_cnidaria_selected"
)
	# Add   "metazoa_selected" back when proteomes are done downloading

DB_NAMES=(
  "swissprot_all"
  "trembl_cnidaria_selected"
)
	# Add   "metazoa_selected" back when proteomes are done downloading


DB_FASTAS=(
  "02_annotation/reference_dbs/swissprot_all.fasta"
  "02_annotation/reference_dbs/trembl_cnidaria_selected.fasta"
)
	# add   "02_annotation/reference_dbs/metazoa_selected.fasta" back after download

WORKER_DIAMOND="scripts/run_uniprot_diamond_single_db.sh"
WORKER_BLASTP="scripts/run_uniprot_blastp_single_db.sh"
WORKER_EXTRACT="scripts/extract_nohit_from_swissprot_diamond.sh"

cat > "$WORKER_DIAMOND" <<'EOF'
#!/bin/bash
#SBATCH --mem=128G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err


PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

module load diamond/2.1.10

THREADS="${SLURM_NTASKS:-8}"

: "${DB_NAME:?Need DB_NAME}"
: "${REF_FASTA:?Need REF_FASTA}"
: "${QUERY:?Need QUERY}"

mkdir -p 02_annotation/db 02_annotation/diamond

if [[ ! -f "$QUERY" ]]; then
    echo "ERROR: Query FASTA not found: $QUERY" >&2
    exit 1
fi
if [[ ! -f "$REF_FASTA" ]]; then
    echo "ERROR: Reference FASTA not found: $REF_FASTA" >&2
    exit 1
fi

DIAMOND_DB="02_annotation/db/${DB_NAME}"
DIAMOND_OUT="02_annotation/diamond/equina_vs_${DB_NAME}.diamond.tsv"
DIAMOND_TOP="02_annotation/diamond/equina_vs_${DB_NAME}.diamond.top_hits.tsv"

echo "Running DIAMOND for $DB_NAME"
diamond makedb --in "$REF_FASTA" --db "$DIAMOND_DB"

diamond blastp \
  --query "$QUERY" \
  --db "$DIAMOND_DB" \
  --out "$DIAMOND_OUT" \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
  --threads "$THREADS"

awk -F '\t' '!seen[$1]++' "$DIAMOND_OUT" > "$DIAMOND_TOP"

wc -l "$DIAMOND_OUT" "$DIAMOND_TOP"
head "$DIAMOND_TOP" || true
EOF

cat > "$WORKER_EXTRACT" <<'EOF'
#!/bin/bash
#SBATCH --mem=32G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=02:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

: "${QUERY:?Need QUERY}"
: "${SWISSPROT_TOP:?Need SWISSPROT_TOP}"
: "${NOHIT_FASTA:?Need NOHIT_FASTA}"
: "${NOHIT_LIST:?Need NOHIT_LIST}"

python3 - <<'PY'
from pathlib import Path

query = Path("00_raw/equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa")
top_hits = Path("02_annotation/diamond/equina_vs_swissprot_all.diamond.top_hits.tsv")
nohit_fasta = Path("02_annotation/intermediate/equina_no_swissprot_diamond_hit.fa")
nohit_list = Path("02_annotation/intermediate/equina_no_swissprot_diamond_hit.ids.txt")

hit_ids = set()
with top_hits.open() as fh:
    for line in fh:
        if line.strip():
            hit_ids.add(line.split('\t', 1)[0])

with query.open() as fh_in, nohit_fasta.open('w') as fh_out, nohit_list.open('w') as fh_ids:
    keep = False
    current_id = None
    for line in fh_in:
        if line.startswith('>'):
            current_id = line[1:].strip().split()[0]
            keep = current_id not in hit_ids
            if keep:
                fh_out.write(line)
                fh_ids.write(current_id + '\n')
        else:
            if keep:
                fh_out.write(line)

print(f"Swiss-Prot DIAMOND hits: {len(hit_ids)}")
print(f"No-hit FASTA written to: {nohit_fasta}")
print(f"No-hit ID list written to: {nohit_list}")
PY

wc -l "$NOHIT_LIST"
grep -c '^>' "$NOHIT_FASTA"
head "$NOHIT_LIST" || true
EOF

cat > "$WORKER_BLASTP" <<'EOF'
#!/bin/bash
#SBATCH --mem=128G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

module load blast+/2.14.1 || true

THREADS="${SLURM_NTASKS:-8}"

: "${DB_NAME:?Need DB_NAME}"
: "${REF_FASTA:?Need REF_FASTA}"
: "${QUERY_NOHIT:?Need QUERY_NOHIT}"

mkdir -p 02_annotation/db 02_annotation/blastp

if [[ ! -f "$QUERY_NOHIT" ]]; then
    echo "ERROR: No-hit query FASTA not found: $QUERY_NOHIT" >&2
    exit 1
fi
if [[ ! -f "$REF_FASTA" ]]; then
    echo "ERROR: Reference FASTA not found: $REF_FASTA" >&2
    exit 1
fi
if ! command -v makeblastdb >/dev/null 2>&1 || ! command -v blastp >/dev/null 2>&1; then
    echo "ERROR: BLAST+ tools not available in PATH" >&2
    exit 1
fi

BLAST_DB="02_annotation/db/${DB_NAME}_blast"
BLAST_OUT="02_annotation/blastp/equina_no_swissprot_diamond_hit_vs_${DB_NAME}.blastp.tsv"
BLAST_TOP="02_annotation/blastp/equina_no_swissprot_diamond_hit_vs_${DB_NAME}.blastp.top_hits.tsv"

echo "Running BLASTp fallback for $DB_NAME"
makeblastdb -in "$REF_FASTA" -dbtype prot -out "$BLAST_DB"

blastp \
  -query "$QUERY_NOHIT" \
  -db "$BLAST_DB" \
  -out "$BLAST_OUT" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -num_threads "$THREADS"

awk -F '\t' '!seen[$1]++' "$BLAST_OUT" > "$BLAST_TOP"

wc -l "$BLAST_OUT" "$BLAST_TOP"
head "$BLAST_TOP" || true
EOF

chmod +x "$WORKER_DIAMOND" "$WORKER_BLASTP" "$WORKER_EXTRACT"

# Validate query and reference FASTAs
if [[ ! -f "$QUERY" ]]; then
    echo "ERROR: Query FASTA not found: $QUERY" >&2
    exit 1
fi

for i in "${!DB_NAMES[@]}"; do
    if [[ ! -f "${DB_FASTAS[$i]}" ]]; then
        echo "WARNING: Missing FASTA for ${DB_NAMES[$i]}: ${DB_FASTAS[$i]}" >&2
    fi
done

DIAMOND_JOB_IDS=()
SWISSPROT_DIAMOND_JOB_ID=""

# Submit DIAMOND for all databases first
if [[ "$RUN_DIAMOND" == "true" ]]; then
    for i in "${!DB_NAMES[@]}"; do
        DB_NAME="${DB_NAMES[$i]}"
        REF_FASTA="${DB_FASTAS[$i]}"

        if [[ ! -f "$REF_FASTA" ]]; then
            echo "Skipping DIAMOND for $DB_NAME because FASTA is missing"
            continue
        fi

        echo "Submitting DIAMOND job for $DB_NAME"
        JOB_ID=$(sbatch --parsable \
          --job-name="diamond_${DB_NAME}" \
          --export=ALL,DB_NAME="$DB_NAME",REF_FASTA="$REF_FASTA",QUERY="$QUERY" \
          "$WORKER_DIAMOND")

        echo "  submitted job $JOB_ID"
        DIAMOND_JOB_IDS+=("$JOB_ID")

        if [[ "$DB_NAME" == "swissprot_all" ]]; then
            SWISSPROT_DIAMOND_JOB_ID="$JOB_ID"
        fi
    done
fi

# After Swiss-Prot DIAMOND finishes, extract no-hit proteins
EXTRACT_JOB_ID=""
if [[ "$RUN_BLASTP_NOHIT" == "true" ]]; then
    if [[ -z "$SWISSPROT_DIAMOND_JOB_ID" ]]; then
        echo "ERROR: Swiss-Prot DIAMOND job was not submitted; cannot derive no-hit subset for BLASTp fallback" >&2
        exit 1
    fi

    echo "Submitting Swiss-Prot no-hit extraction job"
    EXTRACT_JOB_ID=$(sbatch --parsable \
      --dependency=afterok:${SWISSPROT_DIAMOND_JOB_ID} \
      --job-name="extract_nohit_swissprot" \
      --export=ALL,QUERY="$QUERY",SWISSPROT_TOP="02_annotation/diamond/equina_vs_swissprot_all.diamond.top_hits.tsv",NOHIT_FASTA="02_annotation/intermediate/equina_no_swissprot_diamond_hit.fa",NOHIT_LIST="02_annotation/intermediate/equina_no_swissprot_diamond_hit.ids.txt" \
      "$WORKER_EXTRACT")

    echo "  submitted job $EXTRACT_JOB_ID"
fi

# Submit BLASTp only for the no-hit subset against fallback databases
if [[ "$RUN_BLASTP_NOHIT" == "true" ]]; then
    for DB_NAME in "${BLASTP_FALLBACK_DBS[@]}"; do
        REF_FASTA=""
        for i in "${!DB_NAMES[@]}"; do
            if [[ "${DB_NAMES[$i]}" == "$DB_NAME" ]]; then
                REF_FASTA="${DB_FASTAS[$i]}"
                break
            fi
        done

        if [[ -z "$REF_FASTA" || ! -f "$REF_FASTA" ]]; then
            echo "Skipping BLASTp fallback for $DB_NAME because FASTA is missing"
            continue
        fi

        echo "Submitting BLASTp fallback job for $DB_NAME"
        JOB_ID=$(sbatch --parsable \
          --dependency=afterok:${EXTRACT_JOB_ID} \
          --job-name="blastp_${DB_NAME}_nohit" \
          --export=ALL,DB_NAME="$DB_NAME",REF_FASTA="$REF_FASTA",QUERY_NOHIT="02_annotation/intermediate/equina_no_swissprot_diamond_hit.fa" \
          "$WORKER_BLASTP")

        echo "  submitted job $JOB_ID"
    done
fi

echo "Submission complete."
echo "DIAMOND databases searched: ${DB_NAMES[*]}"
echo "BLASTp fallback databases: ${BLASTP_FALLBACK_DBS[*]}"
