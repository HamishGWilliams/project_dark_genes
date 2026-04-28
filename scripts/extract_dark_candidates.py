#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys
from collections import defaultdict


# Allow very large TSV fields.
try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    limit = sys.maxsize
    while True:
        limit = limit // 10
        try:
            csv.field_size_limit(limit)
            break
        except OverflowError:
            continue


DEFAULT_INPUT_TSV = (
    "/uoa/scratch/users/r02hw22/project_dark_genes/"
    "02_annotation/master/full/equina_representative_full.master_annotation.tsv"
)

DEFAULT_PROTEIN_FASTA = (
    "/uoa/scratch/users/r02hw22/project_dark_genes/"
    "02_annotation/input/equina_representative_longest_per_gene.no_stop.fa"
)

DEFAULT_OUTDIR = (
    "/uoa/scratch/users/r02hw22/project_dark_genes/"
    "03_dark_candidates"
)

TARGET_CLASSES = {
    "function_dark_no_current_annotation",
    "function_dark_but_signalp_secretory_candidate",
}

EMPTY_VALUES = {
    "",
    ".",
    "-",
    "NA",
    "N/A",
    "na",
    "n/a",
    "None",
    "none",
    "NULL",
    "null",
    "no_hit",
    "no hit",
    "No hit",
    "NO_HIT",
}

AMBIGUOUS_TERMS = re.compile(
    r"\b("
    r"putative|hypothetical|uncharacteri[sz]ed|unknown|uncertain|ambiguous|"
    r"partial|fragment|low.?confidence|weak|borderline|possible|probable|"
    r"domain.?of.?unknown|DUF|no.?hit|no.?annotation|not.?classified|"
    r"conflict|contradictory|mixed|predicted protein|unnamed protein"
    r")\b",
    re.IGNORECASE,
)

GENERIC_FASTA_SEQ_RE = re.compile(r"^[A-Za-z\*\.\-]+$")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract dark-gene candidate records from a master annotation TSV, "
            "retain ambiguous/manual-review evidence, and extract matching FASTA sequences."
        )
    )

    parser.add_argument(
        "--input",
        default=DEFAULT_INPUT_TSV,
        help="Input master annotation TSV.",
    )

    parser.add_argument(
        "--protein-fasta",
        default=DEFAULT_PROTEIN_FASTA,
        help=(
            "FASTA file used to extract candidate sequences. "
            "Default: equina_representative_longest_per_gene.no_stop.fa"
        ),
    )

    parser.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
        help="Output directory.",
    )

    parser.add_argument(
        "--classification-column",
        default=None,
        help="Optional explicit annotation class column name.",
    )

    parser.add_argument(
        "--id-column",
        default=None,
        help="Optional explicit candidate/protein ID column name.",
    )

    return parser.parse_args()


def open_tsv(path):
    return open(path, "r", newline="", encoding="utf-8", errors="replace")


def is_informative(value):
    return (value or "").strip() not in EMPTY_VALUES


def clean_sequence(seq):
    return re.sub(r"\s+", "", seq or "")


def looks_like_fasta_sequence(seq):
    seq = clean_sequence(seq)
    return len(seq) > 0 and bool(GENERIC_FASTA_SEQ_RE.match(seq))


def unique_preserve_order(values):
    seen = set()
    out = []

    for value in values:
        if value not in seen:
            out.append(value)
            seen.add(value)

    return out


def detect_classification_column(path, explicit=None, scan_limit=10000):
    with open_tsv(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames or []

        if not headers:
            sys.exit(f"ERROR: no header detected in input TSV: {path}")

        if explicit:
            if explicit not in headers:
                sys.exit(f"ERROR: classification column '{explicit}' not found.")
            return explicit

        preferred = [
            "annotation_class",
            "function_classification",
            "functional_classification",
            "dark_classification",
            "annotation_classification",
            "final_classification",
            "classification",
            "class",
            "dark_status",
            "function_class",
            "functional_status",
        ]

        lower_to_original = {h.lower(): h for h in headers}

        for name in preferred:
            if name in lower_to_original:
                return lower_to_original[name]

        scores = defaultdict(int)

        for i, row in enumerate(reader):
            if i >= scan_limit:
                break

            for header in headers:
                value = row.get(header, "")
                if value in TARGET_CLASSES:
                    scores[header] += 1

        if not scores:
            sys.exit(
                "ERROR: could not auto-detect the classification column. "
                "Re-run with --classification-column COLUMN_NAME."
            )

        return max(scores, key=scores.get)


def detect_id_column(headers, explicit=None):
    if explicit:
        if explicit not in headers:
            sys.exit(f"ERROR: ID column '{explicit}' not found.")
        return explicit

    preferred = [
        "protein_id",
        "protein",
        "protein_name",
        "protein_accession",
        "query",
        "query_id",
        "seqid",
        "seq_id",
        "sequence_id",
        "transcript_id",
        "mrna_id",
        "gene_id",
        "locus_tag",
        "id",
        "accession",
        "name",
    ]

    lower_to_original = {h.lower(): h for h in headers}

    for name in preferred:
        if name in lower_to_original:
            return lower_to_original[name]

    return headers[0]


def evidence_columns(headers, class_col, id_col):
    evidence_patterns = re.compile(
        r"("
        r"blast|diamond|hmmer|hmm|pfam|interpro|ipr|eggnog|ortholog|"
        r"kegg|ko|go|cog|signalp|targetp|tmhmm|deeploc|secret|"
        r"transmembrane|tm|domain|motif|annotation|description|product|"
        r"hit|subject|evalue|bitscore|identity|coverage|qcov|scov|"
        r"note|comment|evidence|confidence|status|rank|source|ambiguity"
        r")",
        re.IGNORECASE,
    )

    excluded = {class_col, id_col}

    return [
        header
        for header in headers
        if header not in excluded and evidence_patterns.search(header)
    ]


def summarise_review_evidence(row, ev_cols):
    parts = []
    flags = set()

    for col in ev_cols:
        val = (row.get(col, "") or "").strip()

        if not is_informative(val):
            continue

        parts.append(f"{col}={val}")

        if AMBIGUOUS_TERMS.search(val) or AMBIGUOUS_TERMS.search(col):
            flags.add("ambiguous_or_low_confidence_evidence")

    if parts:
        flags.add("non_classifying_evidence_retained")
    else:
        flags.add("no_extra_evidence_detected")

    return ";".join(sorted(flags)), " | ".join(parts) if parts else "none"


def alias_variants(identifier):
    raw = str(identifier or "").strip()

    if not raw:
        return []

    variants = set()

    def add(value):
        value = str(value or "").strip()
        if value:
            variants.add(value)
            variants.add(value.replace(" ", "_"))

            if "." in value:
                variants.add(value.rsplit(".", 1)[0])

    add(raw)

    first_token = raw.split()[0]
    add(first_token)

    for token in re.split(r"[|,;]", raw):
        add(token)

    for match in re.findall(
        r"(?:^|[\s;])(?:ID|Name|protein_id|transcript_id|gene_id|Parent)=([^;\s]+)",
        raw,
    ):
        add(match)

    expanded = set()

    for value in variants:
        expanded.add(value)

        for prefix in [
            "protein:",
            "transcript:",
            "gene:",
            "cds:",
            "mRNA:",
            "rna:",
            "ID=",
            "Name=",
            "Parent=",
        ]:
            if value.startswith(prefix):
                expanded.add(value[len(prefix):])

    return list(expanded)


def read_fasta_index(path):
    index = {}
    records = {}
    duplicate_aliases = 0

    current_header = None
    current_id = None
    chunks = []

    def store_record(header, seq):
        nonlocal duplicate_aliases

        if header is None:
            return

        seq = clean_sequence(seq)

        if not looks_like_fasta_sequence(seq):
            return

        primary_id = header.split()[0]
        records[primary_id] = {
            "primary_id": primary_id,
            "header": header,
            "sequence": seq,
        }

        aliases = alias_variants(header) + alias_variants(primary_id)

        for alias in aliases:
            if alias in index:
                duplicate_aliases += 1
                continue

            index[alias] = records[primary_id]

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    store_record(current_header, "".join(chunks))

                current_header = line[1:].strip()
                current_id = current_header.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())

    if current_header is not None:
        store_record(current_header, "".join(chunks))

    return index, records, duplicate_aliases


def find_sequence_for_id(identifier, fasta_index):
    for alias in alias_variants(identifier):
        if alias in fasta_index:
            return fasta_index[alias]

    return None


def fasta_safe_id(raw_id):
    clean = str(raw_id).strip()
    clean = re.sub(r"\s+", "_", clean)
    return clean


def wrap_fasta(seq, width=60):
    seq = clean_sequence(seq)
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def write_fasta_record(handle, candidate_id, classification, fasta_record):
    source_id = fasta_record["primary_id"]
    seq = fasta_record["sequence"]

    handle.write(
        f">{fasta_safe_id(candidate_id)} "
        f"dark_class={classification} "
        f"source_fasta_id={source_id}\n"
    )
    handle.write(wrap_fasta(seq))
    handle.write("\n")


def get_candidate_rows(path, class_col):
    with open_tsv(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row_num, row in enumerate(reader, start=2):
            classification = row.get(class_col, "")

            if classification in TARGET_CLASSES:
                yield row_num, row


def main():
    args = parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"ERROR: input TSV not found: {args.input}")

    if not os.path.exists(args.protein_fasta):
        sys.exit(f"ERROR: FASTA file not found: {args.protein_fasta}")

    os.makedirs(args.outdir, exist_ok=True)

    class_col = detect_classification_column(
        args.input,
        explicit=args.classification_column,
    )

    with open_tsv(args.input) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames or []

    id_col = detect_id_column(headers, explicit=args.id_column)
    ev_cols = evidence_columns(headers, class_col, id_col)

    fasta_index, fasta_records, duplicate_aliases = read_fasta_index(args.protein_fasta)

    out_tsv = os.path.join(args.outdir, "equina_dark_candidates.tsv")
    out_review = os.path.join(args.outdir, "equina_dark_candidates.manual_review_evidence.tsv")
    out_missing = os.path.join(args.outdir, "equina_dark_candidates.missing_sequences.tsv")

    out_all_fasta = os.path.join(args.outdir, "equina_dark_candidates.all.fa")

    class_fasta_paths = {
        "function_dark_no_current_annotation": os.path.join(
            args.outdir,
            "equina_dark_candidates.function_dark_no_current_annotation.fa",
        ),
        "function_dark_but_signalp_secretory_candidate": os.path.join(
            args.outdir,
            "equina_dark_candidates.function_dark_but_signalp_secretory_candidate.fa",
        ),
    }

    added_cols = [
        "dark_candidate_category",
        "manual_review_flags",
        "nonclassifying_or_ambiguous_evidence",
        "sequence_status",
        "source_fasta_id",
        "source_fasta_header",
    ]

    out_headers = unique_preserve_order(headers + added_cols)

    review_headers = unique_preserve_order(
        [
            id_col,
            class_col,
            "manual_review_flags",
            "nonclassifying_or_ambiguous_evidence",
        ] + ev_cols
    )

    total_rows = 0
    candidate_rows = 0
    fasta_written = 0
    missing_seq_count = 0
    candidates_with_review_evidence = 0
    class_counts = defaultdict(int)
    class_fasta_counts = defaultdict(int)

    with open_tsv(args.input) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for _ in reader:
            total_rows += 1

    with open(out_tsv, "w", newline="", encoding="utf-8") as tsv_out, \
            open(out_review, "w", newline="", encoding="utf-8") as review_out, \
            open(out_missing, "w", newline="", encoding="utf-8") as missing_out, \
            open(out_all_fasta, "w", encoding="utf-8") as all_fasta:

        class_handles = {
            cls: open(path, "w", encoding="utf-8")
            for cls, path in class_fasta_paths.items()
        }

        try:
            writer = csv.DictWriter(
                tsv_out,
                delimiter="\t",
                fieldnames=out_headers,
                extrasaction="ignore",
            )
            writer.writeheader()

            review_writer = csv.DictWriter(
                review_out,
                delimiter="\t",
                fieldnames=review_headers,
                extrasaction="ignore",
            )
            review_writer.writeheader()

            missing_writer = csv.writer(missing_out, delimiter="\t")
            missing_writer.writerow([id_col, class_col, "reason"])

            for row_num, row in get_candidate_rows(args.input, class_col):
                classification = row.get(class_col, "")
                candidate_rows += 1
                class_counts[classification] += 1

                candidate_id = (row.get(id_col, "") or "").strip()

                if not candidate_id:
                    candidate_id = f"missing_id_tsv_row_{row_num}"

                flags, evidence_summary = summarise_review_evidence(row, ev_cols)

                if "non_classifying_evidence_retained" in flags:
                    candidates_with_review_evidence += 1

                fasta_record = find_sequence_for_id(candidate_id, fasta_index)

                if fasta_record:
                    row["sequence_status"] = "sequence_found_in_external_fasta"
                    row["source_fasta_id"] = fasta_record["primary_id"]
                    row["source_fasta_header"] = fasta_record["header"]

                    write_fasta_record(all_fasta, candidate_id, classification, fasta_record)
                    write_fasta_record(class_handles[classification], candidate_id, classification, fasta_record)

                    fasta_written += 1
                    class_fasta_counts[classification] += 1

                else:
                    row["sequence_status"] = "missing_from_external_fasta"
                    row["source_fasta_id"] = "none"
                    row["source_fasta_header"] = "none"

                    missing_writer.writerow(
                        [candidate_id, classification, "candidate_id_not_found_in_external_fasta"]
                    )
                    missing_seq_count += 1

                row["dark_candidate_category"] = classification
                row["manual_review_flags"] = flags
                row["nonclassifying_or_ambiguous_evidence"] = evidence_summary

                writer.writerow(row)

                review_row = {
                    id_col: candidate_id,
                    class_col: classification,
                    "manual_review_flags": flags,
                    "nonclassifying_or_ambiguous_evidence": evidence_summary,
                }

                for col in ev_cols:
                    review_row[col] = row.get(col, "")

                review_writer.writerow(review_row)

        finally:
            for handle in class_handles.values():
                handle.close()

    summary_path = os.path.join(args.outdir, "equina_dark_candidates.summary.txt")

    with open(summary_path, "w", encoding="utf-8") as summary:
        summary.write("Dark candidate extraction summary\n")
        summary.write("=================================\n")
        summary.write(f"Input TSV: {args.input}\n")
        summary.write(f"Source FASTA: {args.protein_fasta}\n")
        summary.write(f"Output directory: {args.outdir}\n")
        summary.write(f"Classification column: {class_col}\n")
        summary.write(f"Candidate ID column: {id_col}\n")
        summary.write(f"Total input rows scanned: {total_rows}\n")
        summary.write(f"Total FASTA records indexed: {len(fasta_records)}\n")
        summary.write(f"Duplicate FASTA aliases ignored: {duplicate_aliases}\n")
        summary.write(f"Total dark candidates: {candidate_rows}\n")

        for cls in sorted(TARGET_CLASSES):
            summary.write(f"{cls}: {class_counts[cls]}\n")

        summary.write(f"Candidates with retained review evidence: {candidates_with_review_evidence}\n")
        summary.write(f"Evidence columns retained for manual review: {len(ev_cols)}\n")
        summary.write(f"FASTA records written: {fasta_written}\n")

        for cls in sorted(TARGET_CLASSES):
            summary.write(f"FASTA records written [{cls}]: {class_fasta_counts[cls]}\n")

        summary.write(f"Candidates missing sequences: {missing_seq_count}\n")
        summary.write("\nOutputs:\n")
        summary.write(f"- Candidate TSV: {out_tsv}\n")
        summary.write(f"- Manual review evidence TSV: {out_review}\n")
        summary.write(f"- Combined FASTA: {out_all_fasta}\n")

        for cls, path in sorted(class_fasta_paths.items()):
            summary.write(f"- Class FASTA [{cls}]: {path}\n")

        summary.write(f"- Missing sequence report: {out_missing}\n")

    print("Done.")
    print(f"Summary: {summary_path}")
    print(f"Candidate TSV: {out_tsv}")
    print(f"Manual review TSV: {out_review}")
    print(f"Combined FASTA: {out_all_fasta}")
    print(f"FASTA records written: {fasta_written}")
    print(f"Missing sequences: {missing_seq_count}")

    if missing_seq_count > 0:
        print(
            "NOTE: some candidates were not found in the external FASTA. "
            "Check equina_dark_candidates.missing_sequences.tsv.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()