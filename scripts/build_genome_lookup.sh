#!/bin/bash
#SBATCH --job-name=build_genome_lookup
#SBATCH --mem=32G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

cd /uoa/home/r02hw22/sharedscratch/project_dark_genes/

PROJECT_DIR="/uoa/scratch/users/r02hw22/project_dark_genes"

GFF3="${PROJECT_DIR}/00_raw/unpacked/equina.proteins.gff3"
CANDIDATES="${PROJECT_DIR}/03_dark_candidates/equina_dark_candidates.tsv"
OUTDIR="${PROJECT_DIR}/06_genome_lookup"

mkdir -p "${PROJECT_DIR}/scripts" "${OUTDIR}" "${PROJECT_DIR}/logs"

PYTHON_SCRIPT="${PROJECT_DIR}/scripts/build_genome_lookup.py"

cat > "${PYTHON_SCRIPT}" <<'PYTHON'
#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys
import urllib.parse
from collections import defaultdict, Counter


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


TARGET_CLASSES = {
    "function_dark_no_current_annotation",
    "function_dark_but_signalp_secretory_candidate",
}


TRANSCRIPT_TYPES = {
    "mrna",
    "transcript",
    "primary_transcript",
    "processed_transcript",
    "lncrna",
    "ncrna",
    "rrna",
    "trna",
    "snrna",
    "snorna",
    "mirna",
}


PROTEIN_TYPES = {
    "protein",
    "polypeptide",
    "peptide",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Parse GFF3 into gene-transcript-protein-contig lookup tables, "
            "then link dark candidate proteins back to genomic coordinates."
        )
    )

    parser.add_argument("--gff3", required=True, help="Input GFF3 file.")
    parser.add_argument("--candidates", required=True, help="Dark candidate TSV.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--candidate-id-column", default=None, help="Optional candidate protein ID column.")
    parser.add_argument("--candidate-class-column", default=None, help="Optional candidate class column.")

    return parser.parse_args()


def parse_attributes(attr_string):
    attrs = defaultdict(list)

    if attr_string in {"", "."}:
        return attrs

    for part in attr_string.strip().split(";"):
        if not part:
            continue

        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
        else:
            key, value = part, ""

        key = urllib.parse.unquote(key.strip())
        value = urllib.parse.unquote(value.strip())

        if not key:
            continue

        values = [v.strip() for v in value.split(",") if v.strip()]
        if not values:
            values = [""]

        attrs[key].extend(values)

    return attrs


def get_first_attr(attrs, keys):
    for key in keys:
        if key in attrs:
            for value in attrs[key]:
                if value not in {"", "."}:
                    return value
    return None


def get_all_attrs(attrs, keys):
    out = []
    for key in keys:
        for value in attrs.get(key, []):
            if value not in {"", "."}:
                out.append(value)
    return out


def alias_variants(value):
    raw = str(value or "").strip()

    if not raw:
        return set()

    variants = set()

    def add(v):
        v = str(v or "").strip()
        if not v:
            return

        variants.add(v)

        # common GFF3 prefixes
        for prefix in [
            "gene:",
            "transcript:",
            "rna:",
            "mRNA:",
            "mrna:",
            "protein:",
            "cds:",
            "CDS:",
            "ID=",
            "Name=",
            "Parent=",
        ]:
            if v.startswith(prefix):
                variants.add(v[len(prefix):])

        # versionless alias, e.g. abc.1 -> abc
        if "." in v:
            variants.add(v.rsplit(".", 1)[0])

    add(raw)

    # first whitespace-delimited token from FASTA/GFF-like strings
    add(raw.split()[0])

    # split common compound fields
    for token in re.split(r"[|,;]", raw):
        add(token)

    # parse embedded GFF-style key=value fragments
    for match in re.findall(
        r"(?:^|[\s;])(?:ID|Name|Parent|protein_id|transcript_id|gene_id)=([^;\s]+)",
        raw,
    ):
        add(match)

    return variants


def update_minmax(record, contig, start, end, strand):
    if contig and not record.get("contig"):
        record["contig"] = contig

    if strand and strand != "." and not record.get("strand"):
        record["strand"] = strand

    if start is not None:
        if record.get("start") is None or start < record["start"]:
            record["start"] = start

    if end is not None:
        if record.get("end") is None or end > record["end"]:
            record["end"] = end


def blank_gene(gene_id):
    return {
        "gene_id": gene_id,
        "aliases": set(alias_variants(gene_id)),
        "contig": None,
        "start": None,
        "end": None,
        "strand": None,
        "type": "synthetic_gene",
    }


def blank_transcript(transcript_id):
    return {
        "transcript_id": transcript_id,
        "gene_id": None,
        "aliases": set(alias_variants(transcript_id)),
        "protein_aliases": set(),
        "contig": None,
        "start": None,
        "end": None,
        "strand": None,
        "type": "synthetic_transcript",
        "exons": [],
        "cds": [],
    }


def interval_key(interval):
    return (interval["contig"], interval["start"], interval["end"], interval["strand"])


def sort_intervals(intervals):
    return sorted(
        intervals,
        key=lambda x: (
            str(x["contig"]),
            int(x["start"]),
            int(x["end"]),
            str(x["strand"]),
        ),
    )


def unique_intervals(intervals):
    seen = set()
    out = []

    for interval in sort_intervals(intervals):
        key = interval_key(interval)
        if key not in seen:
            out.append(interval)
            seen.add(key)

    return out


def interval_span(intervals):
    if not intervals:
        return None, None, None

    starts = [i["start"] for i in intervals]
    ends = [i["end"] for i in intervals]

    start = min(starts)
    end = max(ends)

    return start, end, end - start + 1


def interval_total_bp(intervals):
    return sum(i["end"] - i["start"] + 1 for i in unique_intervals(intervals))


def format_location(contig, start, end, strand):
    if not contig or start is None or end is None:
        return "NA"
    return f"{contig}:{start}-{end}({strand or '.'})"


def choose_primary(values, fallback="NA"):
    clean = sorted(v for v in values if v and v != ".")
    if clean:
        return clean[0]
    return fallback


def parse_gff3(gff3_path):
    genes = {}
    transcripts = {}
    protein_alias_to_transcripts = defaultdict(set)

    feature_counts = Counter()

    with open(gff3_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            contig, source, feature_type, start, end, score, strand, phase, attr_string = parts

            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue

            feature_type_raw = feature_type
            feature_type = feature_type.lower()
            feature_counts[feature_type] += 1

            attrs = parse_attributes(attr_string)

            if feature_type == "gene":
                gene_id = get_first_attr(attrs, ["ID", "gene_id", "Name", "locus_tag"])

                if not gene_id:
                    gene_id = f"{contig}:{start}-{end}:{strand}:gene"

                genes.setdefault(gene_id, blank_gene(gene_id))
                genes[gene_id]["type"] = feature_type_raw
                update_minmax(genes[gene_id], contig, start, end, strand)

                for key in ["ID", "Name", "gene_id", "locus_tag"]:
                    genes[gene_id]["aliases"].update(get_all_attrs(attrs, [key]))

                continue

            if feature_type in TRANSCRIPT_TYPES or "rna" in feature_type or "transcript" in feature_type:
                transcript_id = get_first_attr(attrs, ["ID", "transcript_id", "Name"])

                if not transcript_id:
                    transcript_id = f"{contig}:{start}-{end}:{strand}:transcript"

                transcripts.setdefault(transcript_id, blank_transcript(transcript_id))
                tx = transcripts[transcript_id]
                tx["type"] = feature_type_raw

                parent_gene = get_first_attr(attrs, ["Parent", "gene_id", "gene", "locus_tag"])
                if parent_gene:
                    tx["gene_id"] = parent_gene
                    genes.setdefault(parent_gene, blank_gene(parent_gene))
                    update_minmax(genes[parent_gene], contig, start, end, strand)

                update_minmax(tx, contig, start, end, strand)

                for key in ["ID", "Name", "transcript_id", "Parent"]:
                    tx["aliases"].update(get_all_attrs(attrs, [key]))

                for key in ["protein_id", "protein", "protein_accession", "Derives_from"]:
                    tx["protein_aliases"].update(get_all_attrs(attrs, [key]))

                continue

            if feature_type in PROTEIN_TYPES:
                protein_id = get_first_attr(attrs, ["ID", "protein_id", "Name", "protein", "protein_accession"])
                parent_transcripts = get_all_attrs(attrs, ["Parent", "Derives_from", "transcript_id"])

                if protein_id and not parent_transcripts:
                    parent_transcripts = get_all_attrs(attrs, ["mRNA", "mrna"])

                for transcript_id in parent_transcripts:
                    transcripts.setdefault(transcript_id, blank_transcript(transcript_id))
                    tx = transcripts[transcript_id]
                    update_minmax(tx, contig, start, end, strand)

                    if protein_id:
                        tx["protein_aliases"].add(protein_id)
                        tx["protein_aliases"].update(alias_variants(protein_id))
                        protein_alias_to_transcripts[protein_id].add(transcript_id)

                continue

            if feature_type in {"exon", "cds"}:
                parent_transcripts = get_all_attrs(attrs, ["Parent", "transcript_id"])

                if not parent_transcripts:
                    parent = get_first_attr(attrs, ["ID"])
                    if parent:
                        parent_transcripts = [parent]

                gene_id = get_first_attr(attrs, ["gene_id", "gene", "locus_tag"])

                protein_values = get_all_attrs(
                    attrs,
                    [
                        "protein_id",
                        "protein",
                        "protein_accession",
                        "Name",
                        "Derives_from",
                    ],
                )

                # CDS IDs can sometimes be the protein-like ID; keep them as weak aliases.
                if feature_type == "cds":
                    protein_values.extend(get_all_attrs(attrs, ["ID"]))

                for transcript_id in parent_transcripts:
                    transcripts.setdefault(transcript_id, blank_transcript(transcript_id))
                    tx = transcripts[transcript_id]

                    if gene_id and not tx["gene_id"]:
                        tx["gene_id"] = gene_id
                        genes.setdefault(gene_id, blank_gene(gene_id))

                    update_minmax(tx, contig, start, end, strand)

                    interval = {
                        "contig": contig,
                        "start": start,
                        "end": end,
                        "strand": strand,
                    }

                    if feature_type == "exon":
                        tx["exons"].append(interval)
                    else:
                        tx["cds"].append(interval)
                        for protein_value in protein_values:
                            tx["protein_aliases"].add(protein_value)
                            tx["protein_aliases"].update(alias_variants(protein_value))
                            protein_alias_to_transcripts[protein_value].add(transcript_id)

                    if tx["gene_id"]:
                        genes.setdefault(tx["gene_id"], blank_gene(tx["gene_id"]))
                        update_minmax(genes[tx["gene_id"]], contig, start, end, strand)

                continue

    # Final alias indexing after all transcript/protein aliases are available.
    for transcript_id, tx in transcripts.items():
        tx["aliases"].update(alias_variants(transcript_id))

        if tx["gene_id"]:
            genes.setdefault(tx["gene_id"], blank_gene(tx["gene_id"]))
            genes[tx["gene_id"]]["aliases"].update(alias_variants(tx["gene_id"]))

        # If no explicit protein alias exists, transcript ID is retained as a fallback.
        if not tx["protein_aliases"]:
            tx["protein_aliases"].add(transcript_id)

        all_aliases = set()
        for value in tx["protein_aliases"]:
            all_aliases.add(value)
            all_aliases.update(alias_variants(value))

        for value in tx["aliases"]:
            all_aliases.add(value)
            all_aliases.update(alias_variants(value))

        tx["all_match_aliases"] = all_aliases

        for alias in all_aliases:
            protein_alias_to_transcripts[alias].add(transcript_id)

    return genes, transcripts, protein_alias_to_transcripts, feature_counts


def make_lookup_rows(genes, transcripts):
    rows = []

    for transcript_id in sorted(transcripts):
        tx = transcripts[transcript_id]
        gene_id = tx.get("gene_id") or "NA"

        exons = unique_intervals(tx.get("exons", []))
        cds = unique_intervals(tx.get("cds", []))

        if tx.get("start") is not None and tx.get("end") is not None:
            transcript_start = tx["start"]
            transcript_end = tx["end"]
            transcript_span_bp = transcript_end - transcript_start + 1
        elif exons:
            transcript_start, transcript_end, transcript_span_bp = interval_span(exons)
        elif cds:
            transcript_start, transcript_end, transcript_span_bp = interval_span(cds)
        else:
            transcript_start, transcript_end, transcript_span_bp = None, None, None

        cds_start, cds_end, cds_span_bp = interval_span(cds)
        cds_total_bp = interval_total_bp(cds) if cds else None

        if exons:
            exon_count = len(exons)
            exon_count_basis = "exon_features"
        elif cds:
            exon_count = len(cds)
            exon_count_basis = "cds_features_used_as_fallback"
        else:
            exon_count = None
            exon_count_basis = "no_exon_or_cds_features"

        contig = tx.get("contig") or "NA"
        strand = tx.get("strand") or "."

        gene = genes.get(gene_id, blank_gene(gene_id))
        gene_start = gene.get("start")
        gene_end = gene.get("end")
        gene_span_bp = gene_end - gene_start + 1 if gene_start is not None and gene_end is not None else None

        primary_protein_id = choose_primary(tx.get("protein_aliases", set()), fallback=transcript_id)

        row = {
            "gene_id": gene_id,
            "transcript_id": transcript_id,
            "protein_id": primary_protein_id,
            "all_protein_aliases": ";".join(sorted(tx.get("protein_aliases", set()))),
            "all_match_aliases": ";".join(sorted(tx.get("all_match_aliases", set()))),
            "contig": contig,
            "strand": strand,
            "transcript_start": transcript_start if transcript_start is not None else "NA",
            "transcript_end": transcript_end if transcript_end is not None else "NA",
            "transcript_span_bp": transcript_span_bp if transcript_span_bp is not None else "NA",
            "cds_start": cds_start if cds_start is not None else "NA",
            "cds_end": cds_end if cds_end is not None else "NA",
            "cds_span_bp": cds_span_bp if cds_span_bp is not None else "NA",
            "cds_total_bp": cds_total_bp if cds_total_bp is not None else "NA",
            "exon_count": exon_count if exon_count is not None else "NA",
            "exon_count_basis": exon_count_basis,
            "cds_part_count": len(cds),
            "scaffold_location": format_location(contig, transcript_start, transcript_end, strand),
            "cds_location": format_location(contig, cds_start, cds_end, strand),
            "gene_start": gene_start if gene_start is not None else "NA",
            "gene_end": gene_end if gene_end is not None else "NA",
            "gene_span_bp": gene_span_bp if gene_span_bp is not None else "NA",
            "gene_location": format_location(contig, gene_start, gene_end, strand),
            "transcript_type": tx.get("type", "NA"),
            "gene_type": gene.get("type", "NA"),
        }

        rows.append(row)

    return rows


def write_tsv(path, rows, headers):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=headers, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def detect_column(headers, preferred, required_name):
    lower_to_original = {h.lower(): h for h in headers}

    for name in preferred:
        if name in lower_to_original:
            return lower_to_original[name]

    raise SystemExit(
        f"ERROR: could not detect {required_name}. "
        f"Available columns: {', '.join(headers)}"
    )


def read_candidates(path, explicit_id_col=None, explicit_class_col=None):
    with open(path, "r", newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames or []

        if not headers:
            raise SystemExit(f"ERROR: no header found in candidate TSV: {path}")

        id_col = explicit_id_col or detect_column(
            headers,
            [
                "protein_id",
                "candidate_id",
                "query_id",
                "seqid",
                "sequence_id",
                "transcript_id",
                "id",
            ],
            "candidate ID column",
        )

        class_col = explicit_class_col or detect_column(
            headers,
            [
                "dark_candidate_category",
                "annotation_class",
                "classification",
                "function_classification",
                "functional_classification",
                "class",
            ],
            "candidate class column",
        )

        rows = list(reader)

    return headers, rows, id_col, class_col


def build_lookup_index(lookup_rows):
    index = defaultdict(list)

    for row in lookup_rows:
        for field in ["protein_id", "transcript_id", "all_protein_aliases", "all_match_aliases"]:
            value = row.get(field, "")

            if not value or value == "NA":
                continue

            for token in str(value).split(";"):
                for alias in alias_variants(token):
                    index[alias].append(row)

    return index


def find_candidate_matches(candidate_row, candidate_id_col, lookup_index):
    search_values = []

    for col in [
        candidate_id_col,
        "protein_id",
        "source_fasta_id",
        "source_fasta_header",
        "transcript_id",
        "gene_id",
        "id",
    ]:
        if col in candidate_row and candidate_row[col] not in {"", ".", "NA", "none"}:
            search_values.append((col, candidate_row[col]))

    seen_lookup_keys = set()
    matches = []

    for col, value in search_values:
        for alias in alias_variants(value):
            if alias in lookup_index:
                for lookup_row in lookup_index[alias]:
                    key = (
                        lookup_row["gene_id"],
                        lookup_row["transcript_id"],
                        lookup_row["protein_id"],
                        lookup_row["contig"],
                        lookup_row["transcript_start"],
                        lookup_row["transcript_end"],
                    )

                    if key not in seen_lookup_keys:
                        seen_lookup_keys.add(key)
                        matches.append((col, alias, lookup_row))

    return matches


def link_candidates(candidate_path, lookup_rows, outdir, explicit_id_col=None, explicit_class_col=None):
    candidate_headers, candidate_rows, candidate_id_col, candidate_class_col = read_candidates(
        candidate_path,
        explicit_id_col=explicit_id_col,
        explicit_class_col=explicit_class_col,
    )

    lookup_index = build_lookup_index(lookup_rows)

    genome_cols = [
        "genome_match_status",
        "genome_matched_by_column",
        "genome_matched_by_alias",
        "gene_id",
        "transcript_id",
        "gff3_protein_id",
        "contig",
        "strand",
        "scaffold_location",
        "gene_location",
        "transcript_start",
        "transcript_end",
        "transcript_span_bp",
        "cds_location",
        "cds_start",
        "cds_end",
        "cds_span_bp",
        "cds_total_bp",
        "exon_count",
        "exon_count_basis",
        "cds_part_count",
    ]

    linked_headers = candidate_headers + [c for c in genome_cols if c not in candidate_headers]

    linked_rows = []
    unmatched_rows = []
    gene_summary = defaultdict(lambda: {
        "candidate_ids": set(),
        "classes": set(),
        "transcript_ids": set(),
        "protein_ids": set(),
        "contig": None,
        "strand": None,
        "gene_location": None,
        "scaffold_locations": set(),
        "exon_counts": [],
        "cds_spans": [],
        "transcript_spans": [],
    })

    match_counter = Counter()

    for candidate in candidate_rows:
        candidate_id = candidate.get(candidate_id_col, "")
        candidate_class = candidate.get(candidate_class_col, "")
        matches = find_candidate_matches(candidate, candidate_id_col, lookup_index)

        if not matches:
            out = dict(candidate)
            out.update({
                "genome_match_status": "unmatched",
                "genome_matched_by_column": "NA",
                "genome_matched_by_alias": "NA",
                "gene_id": "NA",
                "transcript_id": "NA",
                "gff3_protein_id": "NA",
                "contig": "NA",
                "strand": "NA",
                "scaffold_location": "NA",
                "gene_location": "NA",
                "transcript_start": "NA",
                "transcript_end": "NA",
                "transcript_span_bp": "NA",
                "cds_location": "NA",
                "cds_start": "NA",
                "cds_end": "NA",
                "cds_span_bp": "NA",
                "cds_total_bp": "NA",
                "exon_count": "NA",
                "exon_count_basis": "NA",
                "cds_part_count": "NA",
            })
            linked_rows.append(out)
            unmatched_rows.append(out)
            match_counter["unmatched"] += 1
            continue

        status = "matched_unique" if len(matches) == 1 else "matched_multiple"
        match_counter[status] += 1

        for matched_by_col, matched_by_alias, lookup in matches:
            out = dict(candidate)
            out.update({
                "genome_match_status": status,
                "genome_matched_by_column": matched_by_col,
                "genome_matched_by_alias": matched_by_alias,
                "gene_id": lookup["gene_id"],
                "transcript_id": lookup["transcript_id"],
                "gff3_protein_id": lookup["protein_id"],
                "contig": lookup["contig"],
                "strand": lookup["strand"],
                "scaffold_location": lookup["scaffold_location"],
                "gene_location": lookup["gene_location"],
                "transcript_start": lookup["transcript_start"],
                "transcript_end": lookup["transcript_end"],
                "transcript_span_bp": lookup["transcript_span_bp"],
                "cds_location": lookup["cds_location"],
                "cds_start": lookup["cds_start"],
                "cds_end": lookup["cds_end"],
                "cds_span_bp": lookup["cds_span_bp"],
                "cds_total_bp": lookup["cds_total_bp"],
                "exon_count": lookup["exon_count"],
                "exon_count_basis": lookup["exon_count_basis"],
                "cds_part_count": lookup["cds_part_count"],
            })
            linked_rows.append(out)

            gs = gene_summary[lookup["gene_id"]]
            gs["candidate_ids"].add(candidate_id)
            gs["classes"].add(candidate_class)
            gs["transcript_ids"].add(lookup["transcript_id"])
            gs["protein_ids"].add(lookup["protein_id"])
            gs["contig"] = lookup["contig"]
            gs["strand"] = lookup["strand"]
            gs["gene_location"] = lookup["gene_location"]
            gs["scaffold_locations"].add(lookup["scaffold_location"])

            for key, target in [
                ("exon_count", "exon_counts"),
                ("cds_span_bp", "cds_spans"),
                ("transcript_span_bp", "transcript_spans"),
            ]:
                try:
                    value = int(lookup[key])
                    gs[target].append(value)
                except Exception:
                    pass

    linked_path = os.path.join(outdir, "equina_dark_candidates.genome_linked.tsv")
    unmatched_path = os.path.join(outdir, "equina_dark_candidates.genome_unmatched.tsv")
    gene_summary_path = os.path.join(outdir, "equina_dark_candidates.genome_summary_by_gene.tsv")

    write_tsv(linked_path, linked_rows, linked_headers)
    write_tsv(unmatched_path, unmatched_rows, linked_headers)

    gene_summary_rows = []

    for gene_id in sorted(gene_summary):
        gs = gene_summary[gene_id]

        gene_summary_rows.append({
            "gene_id": gene_id,
            "contig": gs["contig"],
            "strand": gs["strand"],
            "gene_location": gs["gene_location"],
            "n_dark_candidate_ids": len(gs["candidate_ids"]),
            "n_dark_candidate_transcripts": len(gs["transcript_ids"]),
            "n_dark_candidate_proteins": len(gs["protein_ids"]),
            "candidate_ids": ";".join(sorted(gs["candidate_ids"])),
            "annotation_classes": ";".join(sorted(gs["classes"])),
            "transcript_ids": ";".join(sorted(gs["transcript_ids"])),
            "protein_ids": ";".join(sorted(gs["protein_ids"])),
            "scaffold_locations": ";".join(sorted(gs["scaffold_locations"])),
            "min_exon_count": min(gs["exon_counts"]) if gs["exon_counts"] else "NA",
            "max_exon_count": max(gs["exon_counts"]) if gs["exon_counts"] else "NA",
            "min_cds_span_bp": min(gs["cds_spans"]) if gs["cds_spans"] else "NA",
            "max_cds_span_bp": max(gs["cds_spans"]) if gs["cds_spans"] else "NA",
            "min_transcript_span_bp": min(gs["transcript_spans"]) if gs["transcript_spans"] else "NA",
            "max_transcript_span_bp": max(gs["transcript_spans"]) if gs["transcript_spans"] else "NA",
        })

    gene_summary_headers = [
        "gene_id",
        "contig",
        "strand",
        "gene_location",
        "n_dark_candidate_ids",
        "n_dark_candidate_transcripts",
        "n_dark_candidate_proteins",
        "candidate_ids",
        "annotation_classes",
        "transcript_ids",
        "protein_ids",
        "scaffold_locations",
        "min_exon_count",
        "max_exon_count",
        "min_cds_span_bp",
        "max_cds_span_bp",
        "min_transcript_span_bp",
        "max_transcript_span_bp",
    ]

    write_tsv(gene_summary_path, gene_summary_rows, gene_summary_headers)

    return {
        "candidate_id_col": candidate_id_col,
        "candidate_class_col": candidate_class_col,
        "linked_path": linked_path,
        "unmatched_path": unmatched_path,
        "gene_summary_path": gene_summary_path,
        "n_candidate_rows": len(candidate_rows),
        "n_linked_rows": len(linked_rows),
        "n_unmatched": len(unmatched_rows),
        "match_counter": match_counter,
    }


def scaffold_summary(linked_path, outdir):
    rows = []

    with open(linked_path, "r", newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        by_contig = defaultdict(lambda: {
            "candidate_ids": set(),
            "genes": set(),
            "transcripts": set(),
            "classes": set(),
            "starts": [],
            "ends": [],
        })

        for row in reader:
            if row.get("genome_match_status") == "unmatched":
                continue

            contig = row.get("contig", "NA")
            if contig == "NA":
                continue

            candidate_id = row.get("protein_id") or row.get("candidate_id") or row.get("id") or ""
            by_contig[contig]["candidate_ids"].add(candidate_id)
            by_contig[contig]["genes"].add(row.get("gene_id", "NA"))
            by_contig[contig]["transcripts"].add(row.get("transcript_id", "NA"))

            cls = row.get("dark_candidate_category") or row.get("annotation_class") or row.get("classification") or "NA"
            by_contig[contig]["classes"].add(cls)

            try:
                by_contig[contig]["starts"].append(int(row["transcript_start"]))
                by_contig[contig]["ends"].append(int(row["transcript_end"]))
            except Exception:
                pass

    for contig in sorted(by_contig):
        item = by_contig[contig]

        if item["starts"] and item["ends"]:
            span = f"{contig}:{min(item['starts'])}-{max(item['ends'])}"
        else:
            span = "NA"

        rows.append({
            "contig": contig,
            "n_dark_candidate_ids": len(item["candidate_ids"]),
            "n_genes": len(item["genes"]),
            "n_transcripts": len(item["transcripts"]),
            "annotation_classes": ";".join(sorted(item["classes"])),
            "candidate_region_span": span,
        })

    path = os.path.join(outdir, "equina_dark_candidates.summary_by_scaffold.tsv")

    headers = [
        "contig",
        "n_dark_candidate_ids",
        "n_genes",
        "n_transcripts",
        "annotation_classes",
        "candidate_region_span",
    ]

    write_tsv(path, rows, headers)

    return path


def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    if not os.path.exists(args.gff3):
        raise SystemExit(f"ERROR: GFF3 not found: {args.gff3}")

    if not os.path.exists(args.candidates):
        raise SystemExit(f"ERROR: candidate TSV not found: {args.candidates}")

    genes, transcripts, protein_alias_to_transcripts, feature_counts = parse_gff3(args.gff3)
    lookup_rows = make_lookup_rows(genes, transcripts)

    lookup_headers = [
        "gene_id",
        "transcript_id",
        "protein_id",
        "all_protein_aliases",
        "all_match_aliases",
        "contig",
        "strand",
        "transcript_start",
        "transcript_end",
        "transcript_span_bp",
        "cds_start",
        "cds_end",
        "cds_span_bp",
        "cds_total_bp",
        "exon_count",
        "exon_count_basis",
        "cds_part_count",
        "scaffold_location",
        "cds_location",
        "gene_start",
        "gene_end",
        "gene_span_bp",
        "gene_location",
        "transcript_type",
        "gene_type",
    ]

    lookup_path = os.path.join(
        args.outdir,
        "equina_gff3_gene_transcript_protein_contig_lookup.tsv",
    )

    write_tsv(lookup_path, lookup_rows, lookup_headers)

    link_info = link_candidates(
        args.candidates,
        lookup_rows,
        args.outdir,
        explicit_id_col=args.candidate_id_column,
        explicit_class_col=args.candidate_class_column,
    )

    scaffold_summary_path = scaffold_summary(link_info["linked_path"], args.outdir)

    summary_path = os.path.join(args.outdir, "equina_genome_lookup.summary.txt")

    with open(summary_path, "w", encoding="utf-8") as summary:
        summary.write("Genome-aware dark candidate lookup summary\n")
        summary.write("==========================================\n")
        summary.write(f"GFF3: {args.gff3}\n")
        summary.write(f"Candidate TSV: {args.candidates}\n")
        summary.write(f"Output directory: {args.outdir}\n")
        summary.write("\nDetected candidate columns:\n")
        summary.write(f"- Candidate ID column: {link_info['candidate_id_col']}\n")
        summary.write(f"- Candidate class column: {link_info['candidate_class_col']}\n")
        summary.write("\nGFF3 parsed counts:\n")
        summary.write(f"- Genes parsed: {len(genes)}\n")
        summary.write(f"- Transcripts parsed: {len(transcripts)}\n")
        summary.write(f"- Lookup rows written: {len(lookup_rows)}\n")
        summary.write(f"- Indexed match aliases: {len(protein_alias_to_transcripts)}\n")
        summary.write("\nCandidate linking:\n")
        summary.write(f"- Candidate rows read: {link_info['n_candidate_rows']}\n")
        summary.write(f"- Linked output rows: {link_info['n_linked_rows']}\n")
        summary.write(f"- Unmatched candidates: {link_info['n_unmatched']}\n")

        for key in sorted(link_info["match_counter"]):
            summary.write(f"- {key}: {link_info['match_counter'][key]}\n")

        summary.write("\nMost common GFF3 feature types:\n")
        for feature_type, count in feature_counts.most_common(25):
            summary.write(f"- {feature_type}: {count}\n")

        summary.write("\nOutputs:\n")
        summary.write(f"- GFF3 lookup table: {lookup_path}\n")
        summary.write(f"- Genome-linked dark candidates: {link_info['linked_path']}\n")
        summary.write(f"- Genome-unmatched dark candidates: {link_info['unmatched_path']}\n")
        summary.write(f"- Gene-level dark candidate summary: {link_info['gene_summary_path']}\n")
        summary.write(f"- Scaffold-level dark candidate summary: {scaffold_summary_path}\n")
        summary.write(f"- Summary: {summary_path}\n")

    print("Done.")
    print(f"Summary: {summary_path}")
    print(f"GFF3 lookup table: {lookup_path}")
    print(f"Genome-linked candidates: {link_info['linked_path']}")
    print(f"Unmatched candidates: {link_info['unmatched_path']}")
    print(f"Gene summary: {link_info['gene_summary_path']}")
    print(f"Scaffold summary: {scaffold_summary_path}")


if __name__ == "__main__":
    main()
PYTHON

chmod +x "${PYTHON_SCRIPT}"

echo "Running genome-aware lookup build..."
echo "GFF3: ${GFF3}"
echo "Candidates: ${CANDIDATES}"
echo "Output directory: ${OUTDIR}"

python3 "${PYTHON_SCRIPT}" \
    --gff3 "${GFF3}" \
    --candidates "${CANDIDATES}" \
    --outdir "${OUTDIR}"

echo
echo "Summary:"
cat "${OUTDIR}/equina_genome_lookup.summary.txt"

echo
echo "Output files:"
ls -lh "${OUTDIR}"
