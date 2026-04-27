#!/usr/bin/env bash

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

MASTER="02_annotation/master/test100/equina_representative_test100.master_annotation.tsv"
SUMMARY="02_annotation/master/test100/equina_representative_test100.master_annotation.summary.txt"

FASTA="02_annotation/eggnog/test_input/equina_representative_test100.no_stop.fa"
REP_LOOKUP="01_qc/isoforms/representative_longest_protein_per_gene.tsv"

DIAMOND_SWISSPROT="02_annotation/diamond/equina_vs_swissprot_all_representative.top_hits.tsv"
DIAMOND_TREMBL_CNIDARIA="02_annotation/diamond/equina_vs_trembl_cnidaria_selected_representative.top_hits.tsv"

BLASTP_SWISSPROT="02_annotation/blastp/equina_vs_swissprot_all.best_hit_representative.tsv"
BLASTP_TREMBL_CNIDARIA="02_annotation/blastp/equina_vs_trembl_cnidaria_selected.best_hit_representative.tsv"

INTERPROSCAN="02_annotation/interproscan/raw/equina_representative_test100.interproscan.tsv"
EGGNOG="02_annotation/eggnog/raw/equina_representative_test100.emapper.annotations"
SIGNALP="02_annotation/signalp/summary/equina_representative_test100.signalp5_summary.tsv"

QC_DIR="02_annotation/master/test100/qc"
REPORT="${QC_DIR}/equina_representative_test100.master_annotation.qc_report.txt"

mkdir -p "$QC_DIR"

python3 - <<'PY'
from pathlib import Path
import csv
from collections import defaultdict, Counter
import math


PROJECT_DIR = Path("/uoa/home/r02hw22/sharedscratch/project_dark_genes")

MASTER = Path("02_annotation/master/test100/equina_representative_test100.master_annotation.tsv")
SUMMARY = Path("02_annotation/master/test100/equina_representative_test100.master_annotation.summary.txt")

FASTA = Path("02_annotation/eggnog/test_input/equina_representative_test100.no_stop.fa")
REP_LOOKUP = Path("01_qc/isoforms/representative_longest_protein_per_gene.tsv")

DIAMOND_SWISSPROT = Path("02_annotation/diamond/equina_vs_swissprot_all_representative.top_hits.tsv")
DIAMOND_TREMBL_CNIDARIA = Path("02_annotation/diamond/equina_vs_trembl_cnidaria_selected_representative.top_hits.tsv")

BLASTP_SWISSPROT = Path("02_annotation/blastp/equina_vs_swissprot_all.best_hit_representative.tsv")
BLASTP_TREMBL_CNIDARIA = Path("02_annotation/blastp/equina_vs_trembl_cnidaria_selected.best_hit_representative.tsv")

INTERPROSCAN = Path("02_annotation/interproscan/raw/equina_representative_test100.interproscan.tsv")
EGGNOG = Path("02_annotation/eggnog/raw/equina_representative_test100.emapper.annotations")
SIGNALP = Path("02_annotation/signalp/summary/equina_representative_test100.signalp5_summary.tsv")

QC_DIR = Path("02_annotation/master/test100/qc")
REPORT = QC_DIR / "equina_representative_test100.master_annotation.qc_report.txt"

ANNOTATION_CLASS_ORDER = [
    "annotated_swissprot_supported",
    "sequence_supported_trembl_cnidaria",
    "domain_supported_interpro",
    "orthology_supported_eggnog",
    "function_dark_but_signalp_secretory_candidate",
    "function_dark_no_current_annotation",
]


def exists_nonempty(path):
    return path.exists() and path.stat().st_size > 0


def write_mismatch_file(name, rows, header):
    path = QC_DIR / name
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def read_fasta_ids(path):
    ids = []
    lengths = {}
    current = None
    chunks = []

    if not path.exists():
        return ids, lengths

    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = len("".join(chunks))
                current = line[1:].strip().split()[0]
                ids.append(current)
                chunks = []
            else:
                chunks.append(line.strip())

        if current is not None:
            lengths[current] = len("".join(chunks))

    return ids, lengths


def load_master(path):
    if not path.exists():
        raise SystemExit(f"ERROR: Master table not found: {path}")

    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)

    by_id = {}
    duplicates = []

    for row in rows:
        pid = row["protein_id"]
        if pid in by_id:
            duplicates.append(pid)
        by_id[pid] = row

    return rows, by_id, duplicates


def load_representative_lookup(path):
    rep_to_gene = {}
    rep_to_isoform_count = {}

    if not exists_nonempty(path):
        return rep_to_gene, rep_to_isoform_count

    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rep_id = row.get("representative_mrna_id", "")
            if not rep_id:
                continue
            rep_to_gene[rep_id] = row.get("gene_id", "")
            rep_to_isoform_count[rep_id] = row.get("isoform_count", "")

    return rep_to_gene, rep_to_isoform_count


def load_top_hits(path):
    hits = {}

    if not exists_nonempty(path):
        return hits

    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue

            qid = parts[0]

            if qid not in hits:
                hits[qid] = {
                    "subject": parts[1],
                    "pident": parts[2],
                    "evalue": parts[10],
                    "bitscore": parts[11],
                }

    return hits


def join_set(values):
    values = sorted(v for v in values if v and v != "-")
    return "|".join(values) if values else ""


def load_interproscan(path):
    data = defaultdict(lambda: {
        "interpro_hit": "no",
        "interpro_member_databases": set(),
        "interpro_signature_accessions": set(),
        "interpro_accessions": set(),
        "interpro_descriptions": set(),
        "go_terms_interpro": set(),
        "pathways_interpro": set(),
    })

    if not exists_nonempty(path):
        return data

    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue

            pid = parts[0]
            rec = data[pid]
            rec["interpro_hit"] = "yes"

            analysis = parts[3] if len(parts) > 3 else "-"
            sig_acc = parts[4] if len(parts) > 4 else "-"
            ipr_acc = parts[11] if len(parts) > 11 else "-"
            ipr_desc = parts[12] if len(parts) > 12 else "-"
            go_terms = parts[13] if len(parts) > 13 else "-"
            pathways = parts[14] if len(parts) > 14 else "-"

            if analysis != "-":
                rec["interpro_member_databases"].add(analysis)
            if sig_acc != "-":
                rec["interpro_signature_accessions"].add(sig_acc)
            if ipr_acc != "-":
                rec["interpro_accessions"].add(ipr_acc)
            if ipr_desc != "-":
                rec["interpro_descriptions"].add(ipr_desc)

            if go_terms != "-":
                for go in go_terms.split("|"):
                    if go:
                        rec["go_terms_interpro"].add(go)

            if pathways != "-":
                for pathway in pathways.split("|"):
                    if pathway:
                        rec["pathways_interpro"].add(pathway)

    final = {}
    for pid, rec in data.items():
        final[pid] = {
            "interpro_hit": rec["interpro_hit"],
            "interpro_member_databases": join_set(rec["interpro_member_databases"]),
            "interpro_signature_accessions": join_set(rec["interpro_signature_accessions"]),
            "interpro_accessions": join_set(rec["interpro_accessions"]),
            "interpro_descriptions": join_set(rec["interpro_descriptions"]),
            "go_terms_interpro": join_set(rec["go_terms_interpro"]),
            "pathways_interpro": join_set(rec["pathways_interpro"]),
        }

    return final


def detect_eggnog_header(path):
    if not path.exists():
        return None

    with path.open() as fh:
        for line in fh:
            if line.startswith("#query"):
                return line.lstrip("#").rstrip("\n").split("\t")

    return None


def load_eggnog(path):
    data = {}

    if not exists_nonempty(path):
        return data

    header = detect_eggnog_header(path)

    fallback_names = [
        "query",
        "seed_ortholog",
        "evalue",
        "score",
        "eggNOG_OGs",
        "max_annot_lvl",
        "COG_category",
        "Description",
        "Preferred_name",
        "GOs",
        "EC",
        "KEGG_ko",
        "KEGG_Pathway",
        "KEGG_Module",
        "KEGG_Reaction",
        "KEGG_rclass",
        "BRITE",
        "KEGG_TC",
        "CAZy",
        "BiGG_Reaction",
        "PFAMs",
    ]

    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if header and len(parts) == len(header):
                row = dict(zip(header, parts))
                qid = row.get("query", parts[0])
            else:
                row = {}
                qid = parts[0]
                for i, name in enumerate(fallback_names):
                    row[name] = parts[i] if i < len(parts) else ""

            data[qid] = {
                "eggnog_hit": "yes",
                "eggnog_seed_ortholog": row.get("seed_ortholog", ""),
                "eggnog_evalue": row.get("evalue", ""),
                "eggnog_score": row.get("score", ""),
                "eggnog_ogs": row.get("eggNOG_OGs", ""),
                "eggnog_max_annot_lvl": row.get("max_annot_lvl", ""),
                "eggnog_cog_category": row.get("COG_category", ""),
                "eggnog_description": row.get("Description", ""),
                "eggnog_preferred_name": row.get("Preferred_name", ""),
                "go_terms_eggnog": row.get("GOs", ""),
                "eggnog_ec": row.get("EC", ""),
                "eggnog_kegg_ko": row.get("KEGG_ko", ""),
                "eggnog_kegg_pathway": row.get("KEGG_Pathway", ""),
                "eggnog_pfam": row.get("PFAMs", ""),
            }

    return data


def load_signalp(path):
    data = {}

    if not exists_nonempty(path):
        return data

    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            raw = line.rstrip("\n")
            parts = raw.split()
            if not parts:
                continue

            pid = parts[0]
            lower = raw.lower()

            is_positive = (
                "signal peptide" in lower
                or "\tsp\t" in lower
                or " sp " in lower
                or " sp(" in lower
                or lower.endswith("\tsp")
                or lower.endswith(" sp")
            )

            data[pid] = {
                "signalp_result_present": "yes",
                "signalp_positive": "yes" if is_positive else "no",
                "signalp_raw": raw,
            }

    return data


def expected_class(row):
    if row["diamond_swissprot_hit"] == "yes" or row["blastp_swissprot_hit"] == "yes":
        return "annotated_swissprot_supported"

    if row["diamond_trembl_cnidaria_hit"] == "yes" or row["blastp_trembl_cnidaria_hit"] == "yes":
        return "sequence_supported_trembl_cnidaria"

    if row["interpro_hit"] == "yes":
        return "domain_supported_interpro"

    if row["eggnog_hit"] == "yes":
        return "orthology_supported_eggnog"

    if row["signalp_positive"] == "yes":
        return "function_dark_but_signalp_secretory_candidate"

    return "function_dark_no_current_annotation"


def compare_hit_source(master_by_id, source_hits, source_name, master_prefix):
    """
    Compare homology source hits against the master table.

    Important:
    For test100 validation, only proteins present in the master table are checked.
    Source-file IDs outside the test100 FASTA are ignored because the DIAMOND/BLASTp
    source files may come from a larger representative/full run.
    """
    mismatches = []

    for pid in sorted(master_by_id):
        row = master_by_id[pid]
        in_source = pid in source_hits

        expected_hit = "yes" if in_source else "no"
        observed_hit = row.get(f"{master_prefix}_hit", "")

        if observed_hit != expected_hit:
            mismatches.append([pid, source_name, "hit_status", expected_hit, observed_hit])
            continue

        if in_source:
            source = source_hits[pid]
            field_map = {
                "subject": f"{master_prefix}_subject",
                "pident": f"{master_prefix}_pident",
                "evalue": f"{master_prefix}_evalue",
                "bitscore": f"{master_prefix}_bitscore",
            }

            for source_field, master_field in field_map.items():
                expected_value = source.get(source_field, "")
                observed_value = row.get(master_field, "")

                if expected_value != observed_value:
                    mismatches.append([
                        pid,
                        source_name,
                        master_field,
                        expected_value,
                        observed_value,
                    ])

    return mismatches

def compare_interpro(master_by_id, interpro):
    fields = [
        "interpro_hit",
        "interpro_member_databases",
        "interpro_signature_accessions",
        "interpro_accessions",
        "interpro_descriptions",
        "go_terms_interpro",
        "pathways_interpro",
    ]

    mismatches = []
    all_ids = sorted(set(master_by_id) | set(interpro))

    for pid in all_ids:
        if pid not in master_by_id:
            mismatches.append([pid, "interproscan", "source_has_id_not_in_master", "", ""])
            continue

        row = master_by_id[pid]
        source = interpro.get(pid, {})

        for field in fields:
            if field == "interpro_hit":
                expected_value = source.get(field, "no")
            else:
                expected_value = source.get(field, "")

            observed_value = row.get(field, "")

            if expected_value != observed_value:
                mismatches.append([pid, "interproscan", field, expected_value, observed_value])

    return mismatches


def compare_eggnog(master_by_id, eggnog):
    fields = [
        "eggnog_hit",
        "eggnog_seed_ortholog",
        "eggnog_evalue",
        "eggnog_score",
        "eggnog_ogs",
        "eggnog_max_annot_lvl",
        "eggnog_cog_category",
        "eggnog_description",
        "eggnog_preferred_name",
        "go_terms_eggnog",
        "eggnog_ec",
        "eggnog_kegg_ko",
        "eggnog_kegg_pathway",
        "eggnog_pfam",
    ]

    mismatches = []
    all_ids = sorted(set(master_by_id) | set(eggnog))

    for pid in all_ids:
        if pid not in master_by_id:
            mismatches.append([pid, "eggnog", "source_has_id_not_in_master", "", ""])
            continue

        row = master_by_id[pid]
        source = eggnog.get(pid, {})

        for field in fields:
            if field == "eggnog_hit":
                expected_value = source.get(field, "no")
            else:
                expected_value = source.get(field, "")

            observed_value = row.get(field, "")

            if expected_value != observed_value:
                mismatches.append([pid, "eggnog", field, expected_value, observed_value])

    return mismatches


def compare_signalp(master_by_id, signalp):
    fields = [
        "signalp_result_present",
        "signalp_positive",
        "signalp_raw",
    ]

    mismatches = []
    all_ids = sorted(set(master_by_id) | set(signalp))

    for pid in all_ids:
        if pid not in master_by_id:
            mismatches.append([pid, "signalp", "source_has_id_not_in_master", "", ""])
            continue

        row = master_by_id[pid]
        source = signalp.get(pid, {})

        for field in fields:
            if field == "signalp_result_present":
                expected_value = source.get(field, "no")
            elif field == "signalp_positive":
                expected_value = source.get(field, "no")
            else:
                expected_value = source.get(field, "")

            observed_value = row.get(field, "")

            if expected_value != observed_value:
                mismatches.append([pid, "signalp", field, expected_value, observed_value])

    return mismatches


def parse_summary_counts(summary_path):
    counts = {}

    if not exists_nonempty(summary_path):
        return counts

    in_table = False

    with summary_path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line.startswith("annotation_class\tcount\tproportion\tpercentage"):
                in_table = True
                continue

            if in_table:
                if not line.strip():
                    continue

                parts = line.split("\t")
                if len(parts) >= 2:
                    try:
                        counts[parts[0]] = int(parts[1])
                    except ValueError:
                        pass

    return counts


def main():
    QC_DIR.mkdir(parents=True, exist_ok=True)

    report_lines = []
    failed_checks = 0

    def add_check(name, passed, detail):
        nonlocal failed_checks
        status = "PASS" if passed else "FAIL"
        if not passed:
            failed_checks += 1
        report_lines.append(f"{status}\t{name}\t{detail}")

    # Load master
    master_rows, master_by_id, duplicate_master_ids = load_master(MASTER)

    # FASTA row universe
    fasta_ids, fasta_lengths = read_fasta_ids(FASTA)

    add_check(
        "master_file_exists",
        exists_nonempty(MASTER),
        str(MASTER),
    )

    add_check(
        "fasta_file_exists",
        exists_nonempty(FASTA),
        str(FASTA),
    )

    add_check(
        "master_row_count_equals_fasta_count",
        len(master_rows) == len(fasta_ids),
        f"master_rows={len(master_rows)} fasta_ids={len(fasta_ids)}",
    )

    add_check(
        "master_has_no_duplicate_protein_ids",
        len(duplicate_master_ids) == 0,
        f"duplicate_count={len(duplicate_master_ids)}",
    )

    master_ids = set(master_by_id)
    fasta_id_set = set(fasta_ids)

    missing_from_master = sorted(fasta_id_set - master_ids)
    extra_in_master = sorted(master_ids - fasta_id_set)

    write_mismatch_file(
        "fasta_ids_missing_from_master.tsv",
        [[x] for x in missing_from_master],
        ["protein_id"],
    )

    write_mismatch_file(
        "master_ids_not_in_fasta.tsv",
        [[x] for x in extra_in_master],
        ["protein_id"],
    )

    add_check(
        "all_fasta_ids_present_in_master",
        len(missing_from_master) == 0,
        f"missing_from_master={len(missing_from_master)}",
    )

    add_check(
        "no_master_ids_outside_fasta",
        len(extra_in_master) == 0,
        f"extra_in_master={len(extra_in_master)}",
    )

    # Protein lengths
    length_mismatches = []

    for pid, row in master_by_id.items():
        expected = str(fasta_lengths.get(pid, ""))
        observed = row.get("protein_length", "")

        if expected != observed:
            length_mismatches.append([pid, expected, observed])

    write_mismatch_file(
        "protein_length_mismatches.tsv",
        length_mismatches,
        ["protein_id", "expected_length_from_fasta", "observed_master_length"],
    )

    add_check(
        "protein_lengths_match_fasta",
        len(length_mismatches) == 0,
        f"mismatches={len(length_mismatches)}",
    )

    # Representative lookup
    rep_to_gene, rep_to_isoform_count = load_representative_lookup(REP_LOOKUP)
    rep_mismatches = []

    if exists_nonempty(REP_LOOKUP):
        for pid, row in master_by_id.items():
            expected_gene = rep_to_gene.get(pid, pid)
            expected_isoform_count = rep_to_isoform_count.get(pid, "1")

            if row.get("gene_id", "") != expected_gene:
                rep_mismatches.append([
                    pid,
                    "gene_id",
                    expected_gene,
                    row.get("gene_id", ""),
                ])

            if row.get("isoform_count", "") != expected_isoform_count:
                rep_mismatches.append([
                    pid,
                    "isoform_count",
                    expected_isoform_count,
                    row.get("isoform_count", ""),
                ])

        write_mismatch_file(
            "representative_lookup_mismatches.tsv",
            rep_mismatches,
            ["protein_id", "field", "expected_from_lookup", "observed_master"],
        )

        add_check(
            "representative_lookup_matches_master",
            len(rep_mismatches) == 0,
            f"mismatches={len(rep_mismatches)}",
        )
    else:
        add_check(
            "representative_lookup_file_exists",
            False,
            f"missing_or_empty={REP_LOOKUP}",
        )

    # Source files presence
    source_files = {
        "diamond_swissprot": DIAMOND_SWISSPROT,
        "diamond_trembl_cnidaria": DIAMOND_TREMBL_CNIDARIA,
        "blastp_swissprot": BLASTP_SWISSPROT,
        "blastp_trembl_cnidaria": BLASTP_TREMBL_CNIDARIA,
        "interproscan": INTERPROSCAN,
        "eggnog": EGGNOG,
        "signalp": SIGNALP,
    }

    for label, path in source_files.items():
        add_check(
            f"{label}_file_exists_nonempty",
            exists_nonempty(path),
            str(path),
        )

    # DIAMOND and BLASTp
    diamond_swiss = load_top_hits(DIAMOND_SWISSPROT)
    diamond_trembl = load_top_hits(DIAMOND_TREMBL_CNIDARIA)
    blastp_swiss = load_top_hits(BLASTP_SWISSPROT)
    blastp_trembl = load_top_hits(BLASTP_TREMBL_CNIDARIA)

    homology_mismatches = []
    homology_mismatches.extend(
        compare_hit_source(
            master_by_id,
            diamond_swiss,
            "diamond_swissprot",
            "diamond_swissprot",
        )
    )
    homology_mismatches.extend(
        compare_hit_source(
            master_by_id,
            diamond_trembl,
            "diamond_trembl_cnidaria",
            "diamond_trembl_cnidaria",
        )
    )
    homology_mismatches.extend(
        compare_hit_source(
            master_by_id,
            blastp_swiss,
            "blastp_swissprot",
            "blastp_swissprot",
        )
    )
    homology_mismatches.extend(
        compare_hit_source(
            master_by_id,
            blastp_trembl,
            "blastp_trembl_cnidaria",
            "blastp_trembl_cnidaria",
        )
    )

    write_mismatch_file(
        "homology_source_mismatches.tsv",
        homology_mismatches,
        ["protein_id", "source", "field", "expected_from_source", "observed_master"],
    )

    add_check(
        "diamond_blastp_fields_match_sources",
        len(homology_mismatches) == 0,
        f"mismatches={len(homology_mismatches)}",
    )

    # InterProScan
    interpro = load_interproscan(INTERPROSCAN)
    interpro_mismatches = compare_interpro(master_by_id, interpro)

    write_mismatch_file(
        "interproscan_mismatches.tsv",
        interpro_mismatches,
        ["protein_id", "source", "field", "expected_from_source", "observed_master"],
    )

    add_check(
        "interproscan_fields_match_source",
        len(interpro_mismatches) == 0,
        f"mismatches={len(interpro_mismatches)}",
    )

    # eggNOG
    eggnog = load_eggnog(EGGNOG)
    eggnog_mismatches = compare_eggnog(master_by_id, eggnog)

    write_mismatch_file(
        "eggnog_mismatches.tsv",
        eggnog_mismatches,
        ["protein_id", "source", "field", "expected_from_source", "observed_master"],
    )

    add_check(
        "eggnog_fields_match_source",
        len(eggnog_mismatches) == 0,
        f"mismatches={len(eggnog_mismatches)}",
    )

    # SignalP
    signalp = load_signalp(SIGNALP)
    signalp_mismatches = compare_signalp(master_by_id, signalp)

    write_mismatch_file(
        "signalp_mismatches.tsv",
        signalp_mismatches,
        ["protein_id", "source", "field", "expected_from_source", "observed_master"],
    )

    add_check(
        "signalp_fields_match_source",
        len(signalp_mismatches) == 0,
        f"mismatches={len(signalp_mismatches)}",
    )

    # Annotation hierarchy
    class_mismatches = []

    for pid, row in master_by_id.items():
        expected = expected_class(row)
        observed = row.get("annotation_class", "")

        if expected != observed:
            class_mismatches.append([pid, expected, observed])

    write_mismatch_file(
        "annotation_class_mismatches.tsv",
        class_mismatches,
        ["protein_id", "expected_annotation_class", "observed_master_annotation_class"],
    )

    add_check(
        "annotation_class_matches_hierarchy",
        len(class_mismatches) == 0,
        f"mismatches={len(class_mismatches)}",
    )

    # Summary counts
    observed_class_counts = Counter(row["annotation_class"] for row in master_rows)
    summary_counts = parse_summary_counts(SUMMARY)

    summary_mismatches = []

    for annotation_class in ANNOTATION_CLASS_ORDER:
        expected = observed_class_counts.get(annotation_class, 0)
        observed = summary_counts.get(annotation_class, None)

        if observed is None:
            summary_mismatches.append([annotation_class, expected, "missing_from_summary"])
        elif observed != expected:
            summary_mismatches.append([annotation_class, expected, observed])

    write_mismatch_file(
        "summary_count_mismatches.tsv",
        summary_mismatches,
        ["annotation_class", "expected_from_master", "observed_in_summary"],
    )

    add_check(
        "summary_counts_match_master",
        len(summary_mismatches) == 0,
        f"mismatches={len(summary_mismatches)}",
    )

    # Coverage summary
    coverage_rows = []

    def yes_count(field):
        return sum(1 for row in master_rows if row.get(field) == "yes")

    total = len(master_rows)

    coverage_specs = [
        ("diamond_swissprot_hit", yes_count("diamond_swissprot_hit")),
        ("blastp_swissprot_hit", yes_count("blastp_swissprot_hit")),
        ("diamond_trembl_cnidaria_hit", yes_count("diamond_trembl_cnidaria_hit")),
        ("blastp_trembl_cnidaria_hit", yes_count("blastp_trembl_cnidaria_hit")),
        ("interpro_hit", yes_count("interpro_hit")),
        ("eggnog_hit", yes_count("eggnog_hit")),
        ("signalp_positive", yes_count("signalp_positive")),
    ]

    for label, count in coverage_specs:
        proportion = count / total if total else 0
        coverage_rows.append([label, count, f"{proportion:.4f}", f"{proportion * 100:.2f}%"])

    write_mismatch_file(
        "source_coverage_summary.tsv",
        coverage_rows,
        ["source_field", "count_yes", "proportion", "percentage"],
    )

    # Annotation class counts
    class_rows = []
    for annotation_class in ANNOTATION_CLASS_ORDER:
        count = observed_class_counts.get(annotation_class, 0)
        proportion = count / total if total else 0
        class_rows.append([annotation_class, count, f"{proportion:.4f}", f"{proportion * 100:.2f}%"])

    write_mismatch_file(
        "annotation_class_summary.tsv",
        class_rows,
        ["annotation_class", "count", "proportion", "percentage"],
    )

    # Final report
    with REPORT.open("w") as fh:
        fh.write("Master annotation QC report\n")
        fh.write("===========================\n\n")
        fh.write(f"Master table: {MASTER}\n")
        fh.write(f"Summary file: {SUMMARY}\n")
        fh.write(f"FASTA row universe: {FASTA}\n")
        fh.write(f"Rows in master: {len(master_rows)}\n")
        fh.write(f"IDs in FASTA: {len(fasta_ids)}\n")
        fh.write(f"Failed checks: {failed_checks}\n\n")

        fh.write("Checks\n")
        fh.write("------\n")
        fh.write("status\tcheck\tdetail\n")
        for line in report_lines:
            fh.write(line + "\n")

        fh.write("\nAnnotation class counts\n")
        fh.write("-----------------------\n")
        fh.write("annotation_class\tcount\tproportion\tpercentage\n")
        for row in class_rows:
            fh.write("\t".join(map(str, row)) + "\n")

        fh.write("\nSource coverage summary\n")
        fh.write("-----------------------\n")
        fh.write("source_field\tcount_yes\tproportion\tpercentage\n")
        for row in coverage_rows:
            fh.write("\t".join(map(str, row)) + "\n")

        fh.write("\nMismatch files\n")
        fh.write("--------------\n")
        for path in sorted(QC_DIR.glob("*mismatch*.tsv")):
            with path.open() as check_fh:
                n_lines = sum(1 for _ in check_fh)
            n_records = max(0, n_lines - 1)
            fh.write(f"{path}\t{n_records} records\n")

    print(f"Wrote QC report: {REPORT}")
    print(f"Failed checks: {failed_checks}")

    if failed_checks > 0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
PY

echo
echo "QC report:"
cat "$REPORT"