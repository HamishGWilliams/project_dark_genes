#!/bin/bash
#SBATCH --mem 48G
#SBATCH --partition uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

# ------------------------------------------------------------
# qc_cds_protein_pairing.sh
#
# Purpose:
#   1. Confirm CDS/protein pairing by FASTA ID
#   2. Generate QC tables for CDS and protein lengths
#   3. Flag short proteins for review
#
# Usage:
#   bash qc_cds_protein_pairing.sh <cds_fasta> <protein_fasta> [output_dir]
#
# Example:
#   bash qc_cds_protein_pairing.sh \
#     /full/path/to/equina...nucl.fa \
#     /full/path/to/equina...aa.fa \
#     /full/path/to/01_qc
#
# Notes:
#   - Uses the first token after ">" as the sequence ID
#   - Resolves all paths to absolute paths
#   - Safe to run from any working directory
# ------------------------------------------------------------

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: bash $0 <cds_fasta> <protein_fasta> [output_dir]" >&2
    exit 1
fi

# Resolve an absolute path for a file or directory path.
# Works whether or not the target already exists.
abs_path() {
    local target="$1"

    if [[ -d "$target" ]]; then
        (cd "$target" && pwd)
    else
        local parent
        parent="$(dirname "$target")"
        local base
        base="$(basename "$target")"
        mkdir -p "$parent"
        (cd "$parent" && printf "%s/%s\n" "$(pwd)" "$base")
    fi
}

# Resolve file path without creating the file itself
abs_file_path() {
    local target="$1"
    local parent
    parent="$(dirname "$target")"
    local base
    base="$(basename "$target")"
    mkdir -p "$parent"
    (cd "$parent" && printf "%s/%s\n" "$(pwd)" "$base")
}

CDS_FASTA="$(abs_file_path "$1")"
PROT_FASTA="$(abs_file_path "$2")"
OUTDIR_RAW="${3:-01_qc}"
OUTDIR="$(abs_path "$OUTDIR_RAW")"

if [[ ! -f "$CDS_FASTA" ]]; then
    echo "ERROR: CDS FASTA not found: $CDS_FASTA" >&2
    exit 1
fi

if [[ ! -f "$PROT_FASTA" ]]; then
    echo "ERROR: Protein FASTA not found: $PROT_FASTA" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

CDS_LEN_TSV="$OUTDIR/cds_lengths.tsv"
PROT_LEN_TSV="$OUTDIR/protein_lengths.tsv"

CDS_SORTED="$OUTDIR/cds_lengths.sorted.tsv"
PROT_SORTED="$OUTDIR/protein_lengths.sorted.tsv"

ONLY_CDS="$OUTDIR/ids_only_in_cds.tsv"
ONLY_PROT="$OUTDIR/ids_only_in_protein.tsv"
SHARED="$OUTDIR/shared_ids.tsv"
PAIRING_TSV="$OUTDIR/cds_protein_pairing_qc.tsv"

LT50="$OUTDIR/proteins_lt50aa.tsv"
LT100="$OUTDIR/proteins_lt100aa.tsv"
LT150="$OUTDIR/proteins_lt150aa.tsv"

REPORT="$OUTDIR/qc_summary.txt"

extract_lengths() {
    local fasta="$1"
    local outfile="$2"

    awk '
    BEGIN {
        OFS="\t"
        id=""
        seqlen=0
    }
    /^>/ {
        if (id != "") {
            print id, seqlen
        }
        id = substr($1, 2)
        seqlen = 0
        next
    }
    {
        gsub(/[[:space:]]/, "", $0)
        seqlen += length($0)
    }
    END {
        if (id != "") {
            print id, seqlen
        }
    }' "$fasta" > "$outfile"
}

echo "Extracting CDS lengths from: $CDS_FASTA"
extract_lengths "$CDS_FASTA" "$CDS_LEN_TSV"

echo "Extracting protein lengths from: $PROT_FASTA"
extract_lengths "$PROT_FASTA" "$PROT_LEN_TSV"

sort -k1,1 "$CDS_LEN_TSV" > "$CDS_SORTED"
sort -k1,1 "$PROT_LEN_TSV" > "$PROT_SORTED"

echo "Comparing IDs..."
join -t $'\t' -1 1 -2 1 "$CDS_SORTED" "$PROT_SORTED" > "$SHARED" || true
join -t $'\t' -1 1 -2 1 -v 1 "$CDS_SORTED" "$PROT_SORTED" > "$ONLY_CDS" || true
join -t $'\t' -1 1 -2 1 -v 2 "$CDS_SORTED" "$PROT_SORTED" > "$ONLY_PROT" || true

echo "Building pairing QC table..."
awk '
BEGIN {
    OFS="\t"
    print "id","cds_len_nt","protein_len_aa","cds_divisible_by_3","expected_aa_if_no_stop","expected_aa_if_terminal_stop_removed","length_match_exact","length_match_minus_stop"
}
{
    id = $1
    cds = $2
    aa  = $3

    div3 = (cds % 3 == 0 ? "yes" : "no")
    exp_exact = int(cds / 3)
    exp_minus = int(cds / 3) - 1

    exact = (exp_exact == aa ? "yes" : "no")
    minus = (exp_minus == aa ? "yes" : "no")

    print id, cds, aa, div3, exp_exact, exp_minus, exact, minus
}' "$SHARED" > "$PAIRING_TSV"

echo "Flagging short proteins..."
awk -F'\t' '$2 < 50  {print}' "$PROT_LEN_TSV" > "$LT50"
awk -F'\t' '$2 < 100 {print}' "$PROT_LEN_TSV" > "$LT100"
awk -F'\t' '$2 < 150 {print}' "$PROT_LEN_TSV" > "$LT150"

CDS_COUNT=$(wc -l < "$CDS_LEN_TSV")
PROT_COUNT=$(wc -l < "$PROT_LEN_TSV")
SHARED_COUNT=$(wc -l < "$SHARED")
ONLY_CDS_COUNT=$(wc -l < "$ONLY_CDS")
ONLY_PROT_COUNT=$(wc -l < "$ONLY_PROT")
LT50_COUNT=$(wc -l < "$LT50")
LT100_COUNT=$(wc -l < "$LT100")
LT150_COUNT=$(wc -l < "$LT150")

DIV3_COUNT=$(awk -F'\t' 'NR>1 && $4=="yes" {c++} END{print c+0}' "$PAIRING_TSV")
EXACT_MATCH_COUNT=$(awk -F'\t' 'NR>1 && $7=="yes" {c++} END{print c+0}' "$PAIRING_TSV")
MINUS_STOP_MATCH_COUNT=$(awk -F'\t' 'NR>1 && $8=="yes" {c++} END{print c+0}' "$PAIRING_TSV")
EITHER_MATCH_COUNT=$(awk -F'\t' 'NR>1 && ($7=="yes" || $8=="yes") {c++} END{print c+0}' "$PAIRING_TSV")

sort -k2,2n "$CDS_LEN_TSV" > "$OUTDIR/cds_lengths.by_length.tsv"
sort -k2,2n "$PROT_LEN_TSV" > "$OUTDIR/protein_lengths.by_length.tsv"

{
    echo "QC summary"
    echo "=========="
    echo
    echo "Input files"
    echo "  CDS FASTA:     $CDS_FASTA"
    echo "  Protein FASTA: $PROT_FASTA"
    echo "  Output dir:    $OUTDIR"
    echo
    echo "Counts"
    echo "  CDS entries:           $CDS_COUNT"
    echo "  Protein entries:       $PROT_COUNT"
    echo "  Shared IDs:            $SHARED_COUNT"
    echo "  IDs only in CDS:       $ONLY_CDS_COUNT"
    echo "  IDs only in protein:   $ONLY_PROT_COUNT"
    echo
    echo "Pairing / translation plausibility"
    echo "  Shared IDs with CDS divisible by 3:        $DIV3_COUNT"
    echo "  Shared IDs with exact length match:        $EXACT_MATCH_COUNT"
    echo "  Shared IDs with CDS/3 - 1 match:           $MINUS_STOP_MATCH_COUNT"
    echo "  Shared IDs matching either expectation:    $EITHER_MATCH_COUNT"
    echo
    echo "Short protein flags"
    echo "  Proteins <50 aa:      $LT50_COUNT"
    echo "  Proteins <100 aa:     $LT100_COUNT"
    echo "  Proteins <150 aa:     $LT150_COUNT"
    echo
    echo "Length summaries"
    echo
    awk '
    BEGIN{count=0; sum=0}
    {a[++count]=$2; sum+=$2}
    END{
        if(count>0){
            print "CDS lengths (nt)"
            print "  count: " count
            print "  min:   " a[1]
            print "  max:   " a[count]
            print "  mean:  " sum/count
            if(count%2==1){print "  median:" a[(count+1)/2]}
            else{print "  median:" (a[count/2]+a[count/2+1])/2}
        }
    }' "$OUTDIR/cds_lengths.by_length.tsv"
    echo
    awk '
    BEGIN{count=0; sum=0}
    {a[++count]=$2; sum+=$2}
    END{
        if(count>0){
            print "Protein lengths (aa)"
            print "  count: " count
            print "  min:   " a[1]
            print "  max:   " a[count]
            print "  mean:  " sum/count
            if(count%2==1){print "  median:" a[(count+1)/2]}
            else{print "  median:" (a[count/2]+a[count/2+1])/2}
        }
    }' "$OUTDIR/protein_lengths.by_length.tsv"
    echo
    echo "Key output files"
    echo "  $CDS_LEN_TSV"
    echo "  $PROT_LEN_TSV"
    echo "  $PAIRING_TSV"
    echo "  $ONLY_CDS"
    echo "  $ONLY_PROT"
    echo "  $LT50"
    echo "  $LT100"
    echo "  $LT150"
} > "$REPORT"

echo
echo "Done."
echo "Main summary: $REPORT"
echo "Inspect:"
echo "  $REPORT"
echo "  $PAIRING_TSV"
echo "  $ONLY_CDS"
echo "  $ONLY_PROT"