#!/usr/bin/env python3

from pathlib import Path
import argparse
import csv
from collections import defaultdict


ANNOTATION_CLASS_ORDER = [
    "annotated_swissprot_supported",
    "sequence_supported_trembl_cnidaria",
    "domain_supported_interpro",
    "orthology_supported_eggnog",
    "function_dark_but_signalp_secretory_candidate",
    "function_dark_no_current_annotation",
]


def parse_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def format_float(value, digits=4):
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def read_fasta_ids(path):
    ids = []
    lengths = {}
    current_id = None
    seq_chunks = []

    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = len("".join(seq_chunks))

                current_id = line[1:].strip().split()[0]
                ids.append(current_id)
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

        if current_id is not None:
            lengths[current_id] = len("".join(seq_chunks))

    return ids, lengths


def load_representative_lookup(path):
    rep_to_gene = {}
    rep_to_isoform_count = {}

    if not path.exists():
        return rep_to_gene, rep_to_isoform_count

    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            rep_id = row.get("representative_mrna_id", "")
            gene_id = row.get("gene_id", "")
            isoform_count = row.get("isoform_count", "")

            if rep_id:
                rep_to_gene[rep_id] = gene_id
                rep_to_isoform_count[rep_id] = isoform_count

    return rep_to_gene, rep_to_isoform_count


def load_top_hit(path, evalue_threshold):
    """
    DIAMOND / BLASTp top-hit table.

    Expected preferred outfmt 6 columns:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

    Legacy 12-column files are tolerated:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

    Only hits with e-value <= evalue_threshold are retained.
    The first retained hit per query is used.
    """
    hits = {}

    if not path.exists():
        return hits

    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 12:
                continue

            qseqid = parts[0]
            sseqid = parts[1]
            pident = parts[2]
            aln_length = parts[3]
            evalue = parts[10]
            bitscore = parts[11]
            qlen = parts[12] if len(parts) > 12 else ""
            slen = parts[13] if len(parts) > 13 else ""

            evalue_float = parse_float(evalue)

            if evalue_float is None:
                continue

            if evalue_float > evalue_threshold:
                continue

            if qseqid in hits:
                continue

            aln_length_float = parse_float(aln_length)
            qlen_float = parse_float(qlen)
            slen_float = parse_float(slen)

            query_coverage = None
            subject_coverage = None

            if aln_length_float is not None and qlen_float not in (None, 0):
                query_coverage = aln_length_float / qlen_float

            if aln_length_float is not None and slen_float not in (None, 0):
                subject_coverage = aln_length_float / slen_float

            hits[qseqid] = {
                "subject": sseqid,
                "pident": pident,
                "aln_length": aln_length,
                "evalue": evalue,
                "bitscore": bitscore,
                "qlen": qlen,
                "slen": slen,
                "query_coverage": format_float(query_coverage),
                "subject_coverage": format_float(subject_coverage),
                "evalue_pass": "yes",
            }

    return hits


def load_interproscan(path):
    """
    InterProScan TSV columns:
    1  protein accession
    2  sequence MD5 digest
    3  sequence length
    4  analysis
    5  signature accession
    6  signature description
    7  start
    8  stop
    9  score
    10 status
    11 date
    12 InterPro accession
    13 InterPro description
    14 GO terms
    15 pathways
    """
    data = defaultdict(
        lambda: {
            "interpro_hit": "no",
            "member_databases": set(),
            "signature_accessions": set(),
            "signature_descriptions": set(),
            "interpro_accessions": set(),
            "interpro_descriptions": set(),
            "go_terms_interpro": set(),
            "pathways_interpro": set(),
        }
    )

    if not path.exists():
        return data

    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 11:
                continue

            protein_id = parts[0]

            analysis = parts[3] if len(parts) > 3 else "-"
            signature_accession = parts[4] if len(parts) > 4 else "-"
            signature_description = parts[5] if len(parts) > 5 else "-"

            interpro_accession = parts[11] if len(parts) > 11 else "-"
            interpro_description = parts[12] if len(parts) > 12 else "-"
            go_terms = parts[13] if len(parts) > 13 else "-"
            pathways = parts[14] if len(parts) > 14 else "-"

            record = data[protein_id]
            record["interpro_hit"] = "yes"

            if analysis != "-":
                record["member_databases"].add(analysis)

            if signature_accession != "-":
                record["signature_accessions"].add(signature_accession)

            if signature_description != "-":
                record["signature_descriptions"].add(signature_description)

            if interpro_accession != "-":
                record["interpro_accessions"].add(interpro_accession)

            if interpro_description != "-":
                record["interpro_descriptions"].add(interpro_description)

            if go_terms != "-":
                for go_term in go_terms.split("|"):
                    if go_term:
                        record["go_terms_interpro"].add(go_term)

            if pathways != "-":
                for pathway in pathways.split("|"):
                    if pathway:
                        record["pathways_interpro"].add(pathway)

    return data


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

    if not path.exists():
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
                query_id = row.get("query", parts[0])
            else:
                row = {}
                query_id = parts[0]

                for i, name in enumerate(fallback_names):
                    row[name] = parts[i] if i < len(parts) else ""

            data[query_id] = {
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

    if not path.exists():
        return data

    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            raw = line.rstrip("\n")
            parts = raw.split()

            if not parts:
                continue

            protein_id = parts[0]
            lower = raw.lower()

            is_positive = (
                "signal peptide" in lower
                or "\tsp\t" in lower
                or " sp " in lower
                or " sp(" in lower
                or lower.endswith("\tsp")
                or lower.endswith(" sp")
            )

            data[protein_id] = {
                "signalp_result_present": "yes",
                "signalp_positive": "yes" if is_positive else "no",
                "signalp_raw": raw,
            }

    return data


def join_set(values):
    cleaned = sorted(v for v in values if v and v != "-")
    return "|".join(cleaned) if cleaned else ""


def add_homology_fields(row, prefix, hit):
    row[f"{prefix}_hit"] = "yes" if hit else "no"
    row[f"{prefix}_subject"] = hit.get("subject", "")
    row[f"{prefix}_pident"] = hit.get("pident", "")
    row[f"{prefix}_aln_length"] = hit.get("aln_length", "")
    row[f"{prefix}_evalue"] = hit.get("evalue", "")
    row[f"{prefix}_bitscore"] = hit.get("bitscore", "")
    row[f"{prefix}_qlen"] = hit.get("qlen", "")
    row[f"{prefix}_slen"] = hit.get("slen", "")
    row[f"{prefix}_query_coverage"] = hit.get("query_coverage", "")
    row[f"{prefix}_subject_coverage"] = hit.get("subject_coverage", "")
    row[f"{prefix}_evalue_pass"] = hit.get("evalue_pass", "no")


def choose_hierarchical_annotation(row):
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


def write_summary(summary_path, fasta_path, rows, input_paths, evalue_threshold):
    counts = defaultdict(int)

    for row in rows:
        counts[row["annotation_class"]] += 1

    total_rows = len(rows)

    with summary_path.open("w") as fh:
        fh.write("Master annotation test100 summary\n")
        fh.write("=================================\n\n")
        fh.write(f"Input FASTA: {fasta_path}\n")
        fh.write(f"Rows written: {total_rows}\n")
        fh.write(f"Homology e-value threshold: {evalue_threshold:g}\n\n")

        fh.write("Input files used:\n")
        for label, path in input_paths.items():
            status = "FOUND" if Path(path).exists() else "MISSING"
            fh.write(f"{label}\t{status}\t{path}\n")

        fh.write("\nAnnotation class counts and proportions:\n")
        fh.write("annotation_class\tcount\tproportion\tpercentage\n")

        for annotation_class in ANNOTATION_CLASS_ORDER:
            count = counts.get(annotation_class, 0)
            proportion = count / total_rows if total_rows else 0
            percentage = proportion * 100

            fh.write(
                f"{annotation_class}\t"
                f"{count}\t"
                f"{proportion:.4f}\t"
                f"{percentage:.2f}%\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Build master annotation table for the 100-protein representative test set."
    )

    parser.add_argument(
        "--fasta",
        default="02_annotation/eggnog/test_input/equina_representative_test100.no_stop.fa",
        help="Representative test FASTA used as the row universe.",
    )

    parser.add_argument(
        "--rep-lookup",
        default="01_qc/isoforms/representative_longest_protein_per_gene.tsv",
        help="Representative protein to gene lookup table.",
    )

    parser.add_argument(
        "--homology-evalue-threshold",
        type=float,
        default=1e-5,
        help="Maximum e-value allowed for DIAMOND/BLASTp hits used in classification.",
    )

    parser.add_argument(
        "--diamond-swissprot",
        default="02_annotation/diamond/filtered/equina_test100_vs_swissprot_all_representative.evalue_1e-5.top_hits.tsv",
        help="Filtered representative DIAMOND Swiss-Prot top hits.",
    )

    parser.add_argument(
        "--diamond-trembl-cnidaria",
        default="02_annotation/diamond/filtered/equina_test100_vs_trembl_cnidaria_selected_representative.evalue_1e-5.top_hits.tsv",
        help="Filtered representative DIAMOND Cnidaria-TrEMBL top hits.",
    )

    parser.add_argument(
        "--blastp-swissprot",
        default="02_annotation/blastp/filtered/equina_test100_vs_swissprot_all.evalue_1e-5.best_hit_representative.tsv",
        help="Filtered representative BLASTp Swiss-Prot best hits.",
    )

    parser.add_argument(
        "--blastp-trembl-cnidaria",
        default="02_annotation/blastp/filtered/equina_test100_vs_trembl_cnidaria_selected.evalue_1e-5.best_hit_representative.tsv",
        help="Filtered representative BLASTp Cnidaria-TrEMBL best hits.",
    )

    parser.add_argument(
        "--interproscan",
        default="02_annotation/interproscan/raw/equina_representative_test100.interproscan.tsv",
        help="Representative test InterProScan TSV.",
    )

    parser.add_argument(
        "--eggnog",
        default="02_annotation/eggnog/raw/equina_representative_test100.emapper.annotations",
        help="Representative test eggNOG-mapper annotations.",
    )

    parser.add_argument(
        "--signalp",
        default="02_annotation/signalp/summary/equina_representative_test100.signalp5_summary.tsv",
        help="Representative test SignalP 5 summary.",
    )

    parser.add_argument(
        "--out",
        default="02_annotation/master/test100/equina_representative_test100.master_annotation.tsv",
        help="Output master annotation table.",
    )

    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not fasta_path.exists():
        raise SystemExit(f"ERROR: FASTA not found: {fasta_path}")

    input_paths = {
        "representative_lookup": args.rep_lookup,
        "diamond_swissprot": args.diamond_swissprot,
        "diamond_trembl_cnidaria": args.diamond_trembl_cnidaria,
        "blastp_swissprot": args.blastp_swissprot,
        "blastp_trembl_cnidaria": args.blastp_trembl_cnidaria,
        "interproscan": args.interproscan,
        "eggnog": args.eggnog,
        "signalp": args.signalp,
    }

    protein_ids, protein_lengths = read_fasta_ids(fasta_path)

    rep_to_gene, rep_to_isoform_count = load_representative_lookup(Path(args.rep_lookup))

    diamond_swissprot = load_top_hit(
        Path(args.diamond_swissprot),
        args.homology_evalue_threshold,
    )
    diamond_trembl_cnidaria = load_top_hit(
        Path(args.diamond_trembl_cnidaria),
        args.homology_evalue_threshold,
    )

    blastp_swissprot = load_top_hit(
        Path(args.blastp_swissprot),
        args.homology_evalue_threshold,
    )
    blastp_trembl_cnidaria = load_top_hit(
        Path(args.blastp_trembl_cnidaria),
        args.homology_evalue_threshold,
    )

    interproscan = load_interproscan(Path(args.interproscan))
    eggnog = load_eggnog(Path(args.eggnog))
    signalp = load_signalp(Path(args.signalp))

    homology_suffixes = [
        "hit",
        "subject",
        "pident",
        "aln_length",
        "evalue",
        "bitscore",
        "qlen",
        "slen",
        "query_coverage",
        "subject_coverage",
        "evalue_pass",
    ]

    columns = [
        "protein_id",
        "gene_id",
        "protein_length",
        "isoform_count",
    ]

    for prefix in [
        "diamond_swissprot",
        "blastp_swissprot",
        "diamond_trembl_cnidaria",
        "blastp_trembl_cnidaria",
    ]:
        columns.extend([f"{prefix}_{suffix}" for suffix in homology_suffixes])

    columns.extend(
        [
            "interpro_hit",
            "interpro_member_databases",
            "interpro_signature_accessions",
            "interpro_accessions",
            "interpro_descriptions",
            "go_terms_interpro",
            "pathways_interpro",
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
            "signalp_result_present",
            "signalp_positive",
            "signalp_raw",
            "annotation_class",
        ]
    )

    rows = []

    for protein_id in protein_ids:
        row = {
            "protein_id": protein_id,
            "gene_id": rep_to_gene.get(protein_id, protein_id),
            "protein_length": protein_lengths.get(protein_id, ""),
            "isoform_count": rep_to_isoform_count.get(protein_id, "1"),
        }

        add_homology_fields(
            row,
            "diamond_swissprot",
            diamond_swissprot.get(protein_id, {}),
        )
        add_homology_fields(
            row,
            "blastp_swissprot",
            blastp_swissprot.get(protein_id, {}),
        )
        add_homology_fields(
            row,
            "diamond_trembl_cnidaria",
            diamond_trembl_cnidaria.get(protein_id, {}),
        )
        add_homology_fields(
            row,
            "blastp_trembl_cnidaria",
            blastp_trembl_cnidaria.get(protein_id, {}),
        )

        interpro_record = interproscan.get(protein_id, {})
        eggnog_record = eggnog.get(protein_id, {})
        signalp_record = signalp.get(protein_id, {})

        row.update(
            {
                "interpro_hit": interpro_record.get("interpro_hit", "no"),
                "interpro_member_databases": join_set(interpro_record.get("member_databases", set())),
                "interpro_signature_accessions": join_set(interpro_record.get("signature_accessions", set())),
                "interpro_accessions": join_set(interpro_record.get("interpro_accessions", set())),
                "interpro_descriptions": join_set(interpro_record.get("interpro_descriptions", set())),
                "go_terms_interpro": join_set(interpro_record.get("go_terms_interpro", set())),
                "pathways_interpro": join_set(interpro_record.get("pathways_interpro", set())),
                "eggnog_hit": eggnog_record.get("eggnog_hit", "no"),
                "eggnog_seed_ortholog": eggnog_record.get("eggnog_seed_ortholog", ""),
                "eggnog_evalue": eggnog_record.get("eggnog_evalue", ""),
                "eggnog_score": eggnog_record.get("eggnog_score", ""),
                "eggnog_ogs": eggnog_record.get("eggnog_ogs", ""),
                "eggnog_max_annot_lvl": eggnog_record.get("eggnog_max_annot_lvl", ""),
                "eggnog_cog_category": eggnog_record.get("eggnog_cog_category", ""),
                "eggnog_description": eggnog_record.get("eggnog_description", ""),
                "eggnog_preferred_name": eggnog_record.get("eggnog_preferred_name", ""),
                "go_terms_eggnog": eggnog_record.get("go_terms_eggnog", ""),
                "eggnog_ec": eggnog_record.get("eggnog_ec", ""),
                "eggnog_kegg_ko": eggnog_record.get("eggnog_kegg_ko", ""),
                "eggnog_kegg_pathway": eggnog_record.get("eggnog_kegg_pathway", ""),
                "eggnog_pfam": eggnog_record.get("eggnog_pfam", ""),
                "signalp_result_present": signalp_record.get("signalp_result_present", "no"),
                "signalp_positive": signalp_record.get("signalp_positive", "no"),
                "signalp_raw": signalp_record.get("signalp_raw", ""),
            }
        )

        row["annotation_class"] = choose_hierarchical_annotation(row)
        rows.append(row)

    with output_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary_path = output_path.with_suffix(".summary.txt")
    write_summary(
        summary_path,
        fasta_path,
        rows,
        input_paths,
        args.homology_evalue_threshold,
    )

    print(f"Wrote master table: {output_path}")
    print(f"Rows written: {len(rows)}")
    print(f"Wrote summary: {summary_path}")


if __name__ == "__main__":
    main()