#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys
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


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Add BUSCO duplication context, near-identical dark-candidate clustering, "
            "duplication interpretation, and priority tiers."
        )
    )

    p.add_argument("--linked", required=True, help="Genome-linked dark candidate TSV.")
    p.add_argument("--lookup", required=True, help="GFF3 gene-transcript-protein-contig lookup TSV.")
    p.add_argument("--dark-fasta", required=True, help="Dark candidate FASTA.")
    p.add_argument("--busco-full-table", required=True, help="BUSCO full_table.tsv or NA.")
    p.add_argument("--cluster-file", required=True, help="MMseqs/CD-HIT/exact cluster output file.")
    p.add_argument("--cluster-method", required=True, help="mmseqs2, cd-hit, or exact_python_fallback.")
    p.add_argument("--outdir", required=True, help="Output directory.")
    p.add_argument("--min-dup-buscos-per-scaffold", type=int, default=3)

    return p.parse_args()


def open_tsv(path):
    return open(path, "r", newline="", encoding="utf-8", errors="replace")


def write_tsv(path, rows, headers):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=headers, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def alias_variants(value):
    raw = str(value or "").strip()

    if not raw or raw in {"NA", "none", "."}:
        return set()

    variants = set()

    def add(v):
        v = str(v or "").strip()
        if not v or v in {"NA", "none", "."}:
            return

        variants.add(v)
        variants.add(v.replace(" ", "_"))

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
            "protein_id=",
            "transcript_id=",
            "gene_id=",
        ]:
            if v.startswith(prefix):
                variants.add(v[len(prefix):])

        if "." in v:
            variants.add(v.rsplit(".", 1)[0])

    add(raw)
    add(raw.split()[0])

    for token in re.split(r"[|,;]", raw):
        add(token)

    for match in re.findall(
        r"(?:^|[\s;])(?:ID|Name|Parent|protein_id|transcript_id|gene_id)=([^;\s]+)",
        raw,
    ):
        add(match)

    return variants


def read_table(path):
    with open_tsv(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames or []
        rows = list(reader)
    return headers, rows


def build_lookup_index(lookup_rows):
    index = defaultdict(list)
    contigs = set()

    for row in lookup_rows:
        contig = row.get("contig")
        if contig and contig not in {"NA", ".", ""}:
            contigs.add(contig)

        for field in [
            "gene_id",
            "transcript_id",
            "protein_id",
            "all_protein_aliases",
            "all_match_aliases",
        ]:
            value = row.get(field, "")
            for token in str(value).split(";"):
                for alias in alias_variants(token):
                    index[alias].append(row)

    return index, contigs


def best_lookup_match(value, lookup_index):
    for alias in alias_variants(value):
        if alias in lookup_index:
            return alias, lookup_index[alias][0]
    return "NA", None


def to_int(value):
    try:
        if value in {"", ".", "NA", None}:
            return None
        return int(float(value))
    except Exception:
        return None


def parse_busco_full_table(path, lookup_index, contigs):
    if path == "NA" or not os.path.exists(path):
        return [], []

    raw_rows = []

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        header = None

        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith("#"):
                candidate = line.lstrip("#").strip()
                cols = candidate.split("\t")

                lower = [c.lower().replace(" ", "_") for c in cols]
                if "busco_id" in lower or "busco" in lower:
                    header = lower

                continue

            parts = line.split("\t")

            if header and len(header) <= len(parts):
                row = {header[i]: parts[i] for i in range(len(header))}
            else:
                # BUSCO v5 commonly uses:
                # BUSCO id, Status, Sequence, Gene Start, Gene End, Score, Length, ...
                row = {
                    "busco_id": parts[0] if len(parts) > 0 else "NA",
                    "status": parts[1] if len(parts) > 1 else "NA",
                    "sequence": parts[2] if len(parts) > 2 else "NA",
                    "gene_start": parts[3] if len(parts) > 3 else "NA",
                    "gene_end": parts[4] if len(parts) > 4 else "NA",
                    "score": parts[5] if len(parts) > 5 else "NA",
                    "length": parts[6] if len(parts) > 6 else "NA",
                }

            raw_rows.append(row)

    loci = []
    duplicated_loci = []

    for row in raw_rows:
        lower_keys = {k.lower().replace(" ", "_"): k for k in row}

        def get_any(*names):
            for name in names:
                key = lower_keys.get(name.lower().replace(" ", "_"))
                if key is not None:
                    return row.get(key, "NA")
            return "NA"

        busco_id = get_any("busco_id", "busco", "busco id")
        status = get_any("status")
        sequence = get_any("sequence", "seq_name", "scaffold", "contig", "target")
        gene = get_any("gene", "gene_id", "prediction", "transcript", "protein")
        start = to_int(get_any("start", "gene_start", "from"))
        end = to_int(get_any("end", "gene_end", "to"))

        mapped_contig = "NA"
        mapped_start = "NA"
        mapped_end = "NA"
        mapped_gene_id = "NA"
        mapped_transcript_id = "NA"
        mapped_protein_id = "NA"
        mapped_by = "unmapped"

        if sequence in contigs:
            mapped_contig = sequence
            mapped_start = start if start is not None else "NA"
            mapped_end = end if end is not None else "NA"
            mapped_by = "busco_sequence_is_contig"

        else:
            match_alias, match = best_lookup_match(sequence, lookup_index)

            if match is None and gene not in {"NA", ".", ""}:
                match_alias, match = best_lookup_match(gene, lookup_index)

            if match is not None:
                mapped_contig = match.get("contig", "NA")
                mapped_start = match.get("transcript_start", "NA")
                mapped_end = match.get("transcript_end", "NA")
                mapped_gene_id = match.get("gene_id", "NA")
                mapped_transcript_id = match.get("transcript_id", "NA")
                mapped_protein_id = match.get("protein_id", "NA")
                mapped_by = f"lookup_alias:{match_alias}"

        out = {
            "busco_id": busco_id,
            "busco_status": status,
            "busco_sequence": sequence,
            "busco_gene_or_prediction": gene,
            "busco_start": start if start is not None else "NA",
            "busco_end": end if end is not None else "NA",
            "mapped_contig": mapped_contig,
            "mapped_start": mapped_start,
            "mapped_end": mapped_end,
            "mapped_gene_id": mapped_gene_id,
            "mapped_transcript_id": mapped_transcript_id,
            "mapped_protein_id": mapped_protein_id,
            "mapping_method": mapped_by,
        }

        loci.append(out)

        if "duplicated" in str(status).lower():
            duplicated_loci.append(out)

    return loci, duplicated_loci


def summarise_busco_scaffolds(duplicated_loci, min_dup_buscos):
    by_scaffold = defaultdict(lambda: {
        "duplicated_busco_ids": set(),
        "duplicated_busco_records": 0,
        "mapped_loci": [],
    })

    for row in duplicated_loci:
        contig = row.get("mapped_contig", "NA")
        if contig in {"NA", "", "."}:
            continue

        by_scaffold[contig]["duplicated_busco_ids"].add(row.get("busco_id", "NA"))
        by_scaffold[contig]["duplicated_busco_records"] += 1

        start = to_int(row.get("mapped_start"))
        end = to_int(row.get("mapped_end"))

        if start is not None and end is not None:
            by_scaffold[contig]["mapped_loci"].append((start, end))

    rows = []
    high_dup_scaffolds = set()

    for contig in sorted(by_scaffold):
        item = by_scaffold[contig]
        starts = [x[0] for x in item["mapped_loci"]]
        ends = [x[1] for x in item["mapped_loci"]]

        n_unique = len(item["duplicated_busco_ids"])
        is_high = n_unique >= min_dup_buscos

        if is_high:
            high_dup_scaffolds.add(contig)

        rows.append({
            "contig": contig,
            "duplicated_BUSCO_unique_ids": n_unique,
            "duplicated_BUSCO_records": item["duplicated_busco_records"],
            "duplicated_BUSCO_ids": ";".join(sorted(item["duplicated_busco_ids"])),
            "duplicated_BUSCO_region_span": (
                f"{contig}:{min(starts)}-{max(ends)}" if starts and ends else "NA"
            ),
            "on_BUSCO_duplication_rich_scaffold": "yes" if is_high else "no",
            "BUSCO_duplication_scaffold_threshold": min_dup_buscos,
        })

    return rows, high_dup_scaffolds


def read_fasta(path):
    records = {}
    current_id = None
    chunks = []

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = "".join(chunks)

                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())

    if current_id is not None:
        records[current_id] = "".join(chunks)

    return records


def write_exact_clusters_from_fasta(fasta_path, out_path):
    records = read_fasta(fasta_path)
    by_seq = defaultdict(list)

    for seq_id, seq in records.items():
        by_seq[seq].append(seq_id)

    with open(out_path, "w", encoding="utf-8") as handle:
        for i, seq in enumerate(sorted(by_seq), start=1):
            cluster_id = f"exact_cluster_{i:06d}"
            for member in sorted(by_seq[seq]):
                handle.write(f"{cluster_id}\t{member}\n")

    return out_path


def parse_mmseqs_cluster_tsv(path):
    clusters = defaultdict(set)

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue

            rep, member = parts[0], parts[1]
            cluster_id = rep
            clusters[cluster_id].add(rep)
            clusters[cluster_id].add(member)

    return clusters


def parse_exact_cluster_tsv(path):
    clusters = defaultdict(set)

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue

            cluster_id, member = parts[0], parts[1]
            clusters[cluster_id].add(member)

    return clusters


def parse_cdhit_clstr(path):
    clusters = defaultdict(set)
    current_cluster = None

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")

            if line.startswith(">Cluster"):
                cluster_num = line.split()[-1]
                current_cluster = f"cdhit_cluster_{int(cluster_num) + 1:06d}"
                continue

            if current_cluster is None:
                continue

            match = re.search(r">([^\.]+(?:\.[^\.]+)*?)\.\.\.", line)
            if match:
                seq_id = match.group(1)
                clusters[current_cluster].add(seq_id)

    return clusters


def load_clusters(cluster_method, cluster_file, dark_fasta, outdir):
    if cluster_method == "exact_python_fallback":
        if not os.path.exists(cluster_file) or os.path.getsize(cluster_file) == 0:
            write_exact_clusters_from_fasta(dark_fasta, cluster_file)
        clusters = parse_exact_cluster_tsv(cluster_file)

    elif cluster_method == "mmseqs2":
        clusters = parse_mmseqs_cluster_tsv(cluster_file)

    elif cluster_method == "cd-hit":
        clusters = parse_cdhit_clstr(cluster_file)

    else:
        clusters = {}

    # Ensure singleton sequences are present if the clustering tool omitted them.
    fasta_records = read_fasta(dark_fasta)
    clustered_members = set()
    for members in clusters.values():
        clustered_members.update(members)

    singleton_i = 1
    for seq_id in fasta_records:
        if seq_id not in clustered_members:
            clusters[f"singleton_{singleton_i:06d}"].add(seq_id)
            singleton_i += 1

    cluster_rows = []
    member_to_cluster = {}

    for i, cluster_id in enumerate(sorted(clusters), start=1):
        members = sorted(clusters[cluster_id])
        stable_cluster_id = f"dark_cluster_{i:06d}"
        representative = members[0]
        size = len(members)

        for member in members:
            member_to_cluster[member] = {
                "cluster_id": stable_cluster_id,
                "representative_id": representative,
                "cluster_size": size,
                "cluster_members": ";".join(members),
            }

            cluster_rows.append({
                "dark_cluster_id": stable_cluster_id,
                "representative_id": representative,
                "candidate_id": member,
                "cluster_size": size,
                "cluster_members": ";".join(members),
                "cluster_method": cluster_method,
            })

    return cluster_rows, member_to_cluster


def candidate_id_from_row(row):
    for col in [
        "protein_id",
        "candidate_id",
        "source_fasta_id",
        "transcript_id",
        "id",
        "query_id",
        "seqid",
    ]:
        value = row.get(col)
        if value and value not in {"NA", ".", "none", ""}:
            return value

    return "NA"


def match_cluster_for_candidate(candidate_id, member_to_cluster):
    for alias in alias_variants(candidate_id):
        if alias in member_to_cluster:
            return member_to_cluster[alias]

    return None


def get_candidate_class(row):
    for col in [
        "dark_candidate_category",
        "annotation_class",
        "classification",
        "function_classification",
        "functional_classification",
        "class",
    ]:
        value = row.get(col)
        if value and value not in {"NA", ".", ""}:
            return value
    return "NA"


def summarise_cluster_context(linked_rows, member_to_cluster, high_dup_scaffolds):
    by_cluster = defaultdict(lambda: {
        "candidate_ids": set(),
        "genes": set(),
        "transcripts": set(),
        "contigs": set(),
        "locations": set(),
        "high_busco_contigs": set(),
        "unmatched_count": 0,
    })

    for row in linked_rows:
        cid = candidate_id_from_row(row)
        cluster = match_cluster_for_candidate(cid, member_to_cluster)

        if cluster is None:
            continue

        cluster_id = cluster["cluster_id"]
        item = by_cluster[cluster_id]
        item["candidate_ids"].add(cid)

        if row.get("genome_match_status") == "unmatched":
            item["unmatched_count"] += 1
            continue

        for field, target in [
            ("gene_id", "genes"),
            ("transcript_id", "transcripts"),
            ("contig", "contigs"),
            ("scaffold_location", "locations"),
        ]:
            value = row.get(field)
            if value and value not in {"NA", ".", ""}:
                item[target].add(value)

        contig = row.get("contig")
        if contig in high_dup_scaffolds:
            item["high_busco_contigs"].add(contig)

    return by_cluster


def interpret_duplication(row, cluster_info, cluster_context, high_dup_scaffolds):
    status = row.get("genome_match_status", "NA")
    contig = row.get("contig", "NA")

    if cluster_info is None:
        cluster_size = 1
        cluster_id = "not_clustered"
    else:
        cluster_size = int(cluster_info["cluster_size"])
        cluster_id = cluster_info["cluster_id"]

    on_high_busco = contig in high_dup_scaffolds
    ctx = cluster_context.get(cluster_id, {})

    n_genes = len(ctx.get("genes", set()))
    n_contigs = len(ctx.get("contigs", set()))
    n_locations = len(ctx.get("locations", set()))
    high_busco_cluster_contigs = len(ctx.get("high_busco_contigs", set()))

    reasons = []

    if status == "unmatched":
        reasons.append("no_genomic_match")

    if on_high_busco:
        reasons.append("candidate_on_BUSCO_duplication_rich_scaffold")

    if cluster_size >= 5:
        reasons.append("large_near_identical_dark_cluster")
    elif cluster_size > 1:
        reasons.append("near_identical_dark_cluster")

    if cluster_size <= 1:
        interpretation = "single_copy_or_no_near_identical_dark_duplicate"

    elif n_genes <= 1 and status != "unmatched":
        interpretation = "likely_annotation_or_isoform_redundancy"

    elif high_busco_cluster_contigs > 0:
        interpretation = "possible_assembly_redundancy_or_haplotig_duplication"

    elif n_genes > 1 and n_locations > 1 and not on_high_busco:
        interpretation = "possible_biological_gene_family_expansion"

    elif n_contigs > 1 and not on_high_busco:
        interpretation = "possible_biological_duplication_requires_synteny_check"

    else:
        interpretation = "ambiguous_duplication_requires_manual_review"

    if interpretation != "single_copy_or_no_near_identical_dark_duplicate":
        reasons.append(interpretation)

    return {
        "cluster_size": cluster_size,
        "cluster_id": cluster_id,
        "n_cluster_genes": n_genes,
        "n_cluster_contigs": n_contigs,
        "n_cluster_locations": n_locations,
        "on_BUSCO_duplication_rich_scaffold": "yes" if on_high_busco else "no",
        "cluster_has_BUSCO_duplication_rich_scaffold": "yes" if high_busco_cluster_contigs > 0 else "no",
        "duplication_interpretation": interpretation,
        "manual_review_reason": ";".join(reasons) if reasons else "none",
    }


def score_candidate(row, dup):
    score = 0
    candidate_class = get_candidate_class(row)

    if candidate_class == "function_dark_but_signalp_secretory_candidate":
        score += 3

    if row.get("genome_match_status") in {"matched_unique", "matched_multiple"}:
        score += 2

    if row.get("sequence_status") == "sequence_found_in_external_fasta":
        score += 2
    elif row.get("source_fasta_id") not in {"NA", "none", "", "."}:
        score += 1

    if row.get("cds_span_bp") not in {"NA", "", "."}:
        score += 1

    if row.get("transcript_span_bp") not in {"NA", "", "."}:
        score += 1

    try:
        exon_count = int(row.get("exon_count", "0"))
        if exon_count >= 2:
            score += 1
    except Exception:
        pass

    if dup["cluster_size"] == 1:
        score += 1
    elif dup["cluster_size"] >= 5:
        score -= 3
    elif dup["cluster_size"] > 1:
        score -= 1

    if dup["on_BUSCO_duplication_rich_scaffold"] == "yes":
        score -= 3

    if dup["duplication_interpretation"] == "likely_annotation_or_isoform_redundancy":
        score -= 3
    elif dup["duplication_interpretation"] == "possible_assembly_redundancy_or_haplotig_duplication":
        score -= 3
    elif dup["duplication_interpretation"] == "ambiguous_duplication_requires_manual_review":
        score -= 2
    elif dup["duplication_interpretation"] == "possible_biological_gene_family_expansion":
        score += 1

    if row.get("genome_match_status") == "unmatched":
        score -= 3

    if score >= 5:
        tier = "high_priority"
    elif score >= 2:
        tier = "medium_priority"
    else:
        tier = "low_priority_or_manual_review"

    return score, tier


def scaffold_counts_from_summary(scaffold_rows):
    by_contig = {}

    for row in scaffold_rows:
        by_contig[row["contig"]] = row

    return by_contig


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    lookup_headers, lookup_rows = read_table(args.lookup)
    linked_headers, linked_rows = read_table(args.linked)

    lookup_index, contigs = build_lookup_index(lookup_rows)

    all_busco_loci, duplicated_busco_loci = parse_busco_full_table(
        args.busco_full_table,
        lookup_index,
        contigs,
    )

    busco_loci_path = os.path.join(args.outdir, "equina_busco_duplicate_loci.tsv")
    busco_loci_headers = [
        "busco_id",
        "busco_status",
        "busco_sequence",
        "busco_gene_or_prediction",
        "busco_start",
        "busco_end",
        "mapped_contig",
        "mapped_start",
        "mapped_end",
        "mapped_gene_id",
        "mapped_transcript_id",
        "mapped_protein_id",
        "mapping_method",
    ]

    write_tsv(busco_loci_path, duplicated_busco_loci, busco_loci_headers)

    scaffold_rows, high_dup_scaffolds = summarise_busco_scaffolds(
        duplicated_busco_loci,
        args.min_dup_buscos_per_scaffold,
    )

    scaffold_path = os.path.join(args.outdir, "equina_scaffold_busco_duplication_summary.tsv")
    scaffold_headers = [
        "contig",
        "duplicated_BUSCO_unique_ids",
        "duplicated_BUSCO_records",
        "duplicated_BUSCO_ids",
        "duplicated_BUSCO_region_span",
        "on_BUSCO_duplication_rich_scaffold",
        "BUSCO_duplication_scaffold_threshold",
    ]
    write_tsv(scaffold_path, scaffold_rows, scaffold_headers)

    cluster_rows, member_to_cluster = load_clusters(
        args.cluster_method,
        args.cluster_file,
        args.dark_fasta,
        args.outdir,
    )

    cluster_path = os.path.join(args.outdir, "equina_dark_candidate_protein_clusters.tsv")
    cluster_headers = [
        "dark_cluster_id",
        "representative_id",
        "candidate_id",
        "cluster_size",
        "cluster_members",
        "cluster_method",
    ]
    write_tsv(cluster_path, cluster_rows, cluster_headers)

    cluster_context = summarise_cluster_context(
        linked_rows,
        member_to_cluster,
        high_dup_scaffolds,
    )

    scaffold_summary_by_contig = scaffold_counts_from_summary(scaffold_rows)

    added_headers = [
        "BUSCO_duplicate_scaffold_unique_ids",
        "BUSCO_duplicate_scaffold_records",
        "on_BUSCO_duplication_rich_scaffold",
        "dark_cluster_id",
        "dark_cluster_representative_id",
        "dark_cluster_size",
        "dark_cluster_members",
        "n_cluster_genes",
        "n_cluster_contigs",
        "n_cluster_locations",
        "cluster_has_BUSCO_duplication_rich_scaffold",
        "duplication_interpretation",
        "priority_score",
        "priority_tier",
        "manual_review_reason",
    ]

    out_headers = linked_headers + [h for h in added_headers if h not in linked_headers]

    context_rows = []
    prioritised_rows = []

    for row in linked_rows:
        out = dict(row)
        candidate_id = candidate_id_from_row(row)
        cluster_info = match_cluster_for_candidate(candidate_id, member_to_cluster)

        if cluster_info is None:
            cluster_info = {
                "cluster_id": "not_clustered",
                "representative_id": candidate_id,
                "cluster_size": 1,
                "cluster_members": candidate_id,
            }

        dup = interpret_duplication(
            row,
            cluster_info,
            cluster_context,
            high_dup_scaffolds,
        )

        score, tier = score_candidate(row, dup)

        contig = row.get("contig", "NA")
        scaffold_info = scaffold_summary_by_contig.get(contig, {})

        out["BUSCO_duplicate_scaffold_unique_ids"] = scaffold_info.get("duplicated_BUSCO_unique_ids", 0)
        out["BUSCO_duplicate_scaffold_records"] = scaffold_info.get("duplicated_BUSCO_records", 0)
        out["on_BUSCO_duplication_rich_scaffold"] = dup["on_BUSCO_duplication_rich_scaffold"]
        out["dark_cluster_id"] = cluster_info["cluster_id"]
        out["dark_cluster_representative_id"] = cluster_info["representative_id"]
        out["dark_cluster_size"] = cluster_info["cluster_size"]
        out["dark_cluster_members"] = cluster_info["cluster_members"]
        out["n_cluster_genes"] = dup["n_cluster_genes"]
        out["n_cluster_contigs"] = dup["n_cluster_contigs"]
        out["n_cluster_locations"] = dup["n_cluster_locations"]
        out["cluster_has_BUSCO_duplication_rich_scaffold"] = dup["cluster_has_BUSCO_duplication_rich_scaffold"]
        out["duplication_interpretation"] = dup["duplication_interpretation"]
        out["priority_score"] = score
        out["priority_tier"] = tier
        out["manual_review_reason"] = dup["manual_review_reason"]

        context_rows.append(out)

    def priority_sort_key(row):
        tier_rank = {
            "high_priority": 0,
            "medium_priority": 1,
            "low_priority_or_manual_review": 2,
        }
        try:
            score = int(row.get("priority_score", 0))
        except Exception:
            score = 0

        return (
            tier_rank.get(row.get("priority_tier"), 9),
            -score,
            row.get("contig", ""),
            row.get("transcript_start", "NA"),
            candidate_id_from_row(row),
        )

    prioritised_rows = sorted(context_rows, key=priority_sort_key)

    context_path = os.path.join(args.outdir, "equina_dark_candidates.with_busco_duplication_context.tsv")
    prioritised_path = os.path.join(args.outdir, "equina_dark_candidates.prioritised.tsv")

    write_tsv(context_path, context_rows, out_headers)
    write_tsv(prioritised_path, prioritised_rows, out_headers)

    interpretation_counts = Counter(row["duplication_interpretation"] for row in context_rows)
    tier_counts = Counter(row["priority_tier"] for row in context_rows)
    cluster_size_counts = Counter()

    seen_clusters = {}
    for row in cluster_rows:
        seen_clusters[row["dark_cluster_id"]] = int(row["cluster_size"])

    for size in seen_clusters.values():
        if size == 1:
            cluster_size_counts["singleton_clusters"] += 1
        elif size < 5:
            cluster_size_counts["small_near_identical_clusters_2_to_4"] += 1
        else:
            cluster_size_counts["large_near_identical_clusters_5_plus"] += 1

    summary_path = os.path.join(args.outdir, "equina_duplication_context.summary.txt")

    with open(summary_path, "w", encoding="utf-8") as summary:
        summary.write("Dark candidate BUSCO/duplication context summary\n")
        summary.write("================================================\n")
        summary.write(f"Genome-linked candidate TSV: {args.linked}\n")
        summary.write(f"GFF3 lookup TSV: {args.lookup}\n")
        summary.write(f"Dark FASTA: {args.dark_fasta}\n")
        summary.write(f"BUSCO full table: {args.busco_full_table}\n")
        summary.write(f"Cluster method: {args.cluster_method}\n")
        summary.write(f"Cluster file: {args.cluster_file}\n")
        summary.write(f"BUSCO duplicated-scaffold threshold: {args.min_dup_buscos_per_scaffold} unique duplicated BUSCO IDs\n")
        summary.write("\nBUSCO duplication context:\n")
        summary.write(f"- Total BUSCO rows parsed: {len(all_busco_loci)}\n")
        summary.write(f"- Duplicated BUSCO loci parsed: {len(duplicated_busco_loci)}\n")
        summary.write(f"- BUSCO-duplication-rich scaffolds: {len(high_dup_scaffolds)}\n")
        summary.write("\nProtein clustering:\n")
        summary.write(f"- Cluster rows written: {len(cluster_rows)}\n")
        summary.write(f"- Unique clusters: {len(seen_clusters)}\n")
        for key in sorted(cluster_size_counts):
            summary.write(f"- {key}: {cluster_size_counts[key]}\n")
        summary.write("\nDuplication interpretation counts:\n")
        for key, value in interpretation_counts.most_common():
            summary.write(f"- {key}: {value}\n")
        summary.write("\nPriority tier counts:\n")
        for key, value in tier_counts.most_common():
            summary.write(f"- {key}: {value}\n")
        summary.write("\nOutputs:\n")
        summary.write(f"- Duplicated BUSCO loci: {busco_loci_path}\n")
        summary.write(f"- Scaffold BUSCO duplication summary: {scaffold_path}\n")
        summary.write(f"- Dark candidate protein clusters: {cluster_path}\n")
        summary.write(f"- Candidates with BUSCO/duplication context: {context_path}\n")
        summary.write(f"- Prioritised candidates: {prioritised_path}\n")
        summary.write(f"- Summary: {summary_path}\n")

    print("Done.")
    print(f"Summary: {summary_path}")
    print(f"Duplicated BUSCO loci: {busco_loci_path}")
    print(f"Scaffold BUSCO duplication summary: {scaffold_path}")
    print(f"Dark candidate protein clusters: {cluster_path}")
    print(f"Candidates with BUSCO/duplication context: {context_path}")
    print(f"Prioritised candidates: {prioritised_path}")


if __name__ == "__main__":
    main()
