#!/usr/bin/env python3
"""Add DMR/methylation interval overlap context to dark candidates.

Inputs:
  --candidates : candidate TSV with contig/start/end columns
  --dmrs       : BED-like or TSV DMR interval file

Outputs:
  <prefix>.dmr_overlap.tsv
  <prefix>.dmr_overlap.summary.txt

Candidate and generic TSV coordinates are treated as 1-based closed intervals.
BED coordinates can be interpreted as 0-based half-open with --dmr-format bed.
"""

import argparse
import bisect
import csv
import os
import sys
from collections import Counter, defaultdict

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    limit = sys.maxsize
    while True:
        limit //= 10
        try:
            csv.field_size_limit(limit)
            break
        except OverflowError:
            continue


def parse_args():
    p = argparse.ArgumentParser(description="Annotate dark candidates with DMR interval overlap context.")
    p.add_argument("--candidates", required=True)
    p.add_argument("--dmrs", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--prefix", default="equina_dark_candidates")
    p.add_argument("--candidate-id-column", default="protein_id")
    p.add_argument("--contig-column", default="contig")
    p.add_argument("--start-column", default="transcript_start")
    p.add_argument("--end-column", default="transcript_end")
    p.add_argument("--dmr-format", choices=["auto", "bed", "tsv"], default="auto")
    p.add_argument("--dmr-contig-column", default=None)
    p.add_argument("--dmr-start-column", default=None)
    p.add_argument("--dmr-end-column", default=None)
    p.add_argument("--dmr-id-column", default=None)
    p.add_argument("--dmr-effect-column", default=None, help="Optional column such as hyper/hypo/delta methylation.")
    p.add_argument("--min-overlap-bp", type=int, default=1)
    return p.parse_args()


def clean(x):
    x = str(x or "").strip()
    return "" if x in {"", ".", "NA", "none", "None", "nan"} else x


def to_int(x):
    x = clean(x)
    if not x:
        return None
    try:
        return int(float(x.replace(",", "")))
    except Exception:
        return None


def detect_col(headers, explicit, candidates, label, required=True):
    if explicit:
        if explicit not in headers:
            raise SystemExit(f"ERROR: requested {label} column {explicit!r} not found")
        return explicit
    lower = {h.lower(): h for h in headers}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    if required:
        raise SystemExit(f"ERROR: could not detect {label}. Available columns: {', '.join(headers)}")
    return None


def looks_like_header(parts):
    if len(parts) < 3:
        return False
    return to_int(parts[1]) is None or to_int(parts[2]) is None


def sign_effect(value):
    value = clean(value)
    if not value:
        return "NA"
    try:
        numeric = float(value)
    except Exception:
        return value
    if numeric > 0:
        return f"hyper_or_positive_meth_diff:{value}"
    if numeric < 0:
        return f"hypo_or_negative_meth_diff:{value}"
    return f"zero_meth_diff:{value}"


def read_dmrs(path, args):
    by_contig = defaultdict(list)
    parse_counts = Counter()
    detected = {}

    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        first_data = None
        for raw in handle:
            if raw.strip() and not raw.startswith("#"):
                first_data = raw.rstrip("\n")
                break
        if first_data is None:
            return by_contig, {}, parse_counts, detected

    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        if args.dmr_format == "bed" or (args.dmr_format == "auto" and not looks_like_header(first_data.split("\t"))):
            detected.update({"dmr_format_used": "bed"})
            for i, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith("#"):
                    parse_counts["comment_or_blank"] += 1
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    parse_counts["bad_bed_line"] += 1
                    continue
                contig = parts[0]
                start0 = to_int(parts[1])
                end0 = to_int(parts[2])
                if start0 is None or end0 is None:
                    parse_counts["bad_coordinates"] += 1
                    continue
                dmr_id = parts[3] if len(parts) > 3 and clean(parts[3]) else f"dmr_{i}"
                effect = parts[4] if len(parts) > 4 and clean(parts[4]) else "NA"
                start = start0 + 1
                end = end0
                if end < start:
                    parse_counts["bad_interval"] += 1
                    continue
                by_contig[contig].append({"contig": contig, "start": start, "end": end, "dmr_id": dmr_id, "effect": effect})
                parse_counts["parsed_bed"] += 1
        else:
            detected.update({"dmr_format_used": "tsv"})
            reader = csv.DictReader(handle, delimiter="\t")
            headers = reader.fieldnames or []
            contig_col = detect_col(
                headers,
                args.dmr_contig_column,
                ["contig", "scaffold", "seqid", "seqnames", "chrom", "chromosome"],
                "DMR contig",
            )
            start_col = detect_col(headers, args.dmr_start_column, ["start", "dmr_start", "begin"], "DMR start")
            end_col = detect_col(headers, args.dmr_end_column, ["end", "dmr_end", "stop"], "DMR end")
            id_col = detect_col(
                headers,
                args.dmr_id_column,
                ["dmr_id", "id", "ID", "Name", "name", "region", "Parent"],
                "DMR ID",
                required=False,
            )
            effect_col = detect_col(
                headers,
                args.dmr_effect_column,
                ["direction", "effect", "status", "meth.diff", "meth_diff", "methylation_difference", "delta_methylation", "logfc"],
                "DMR effect",
                required=False,
            )
            detected.update({
                "dmr_contig_column": contig_col,
                "dmr_start_column": start_col,
                "dmr_end_column": end_col,
                "dmr_id_column": id_col or "NA",
                "dmr_effect_column": effect_col or "NA",
            })
            for i, row in enumerate(reader, start=1):
                contig = clean(row.get(contig_col))
                start = to_int(row.get(start_col))
                end = to_int(row.get(end_col))
                if not contig or start is None or end is None:
                    parse_counts["bad_coordinates"] += 1
                    continue
                if end < start:
                    start, end = end, start
                dmr_id = clean(row.get(id_col)) if id_col else f"dmr_{i}"
                raw_effect = clean(row.get(effect_col)) if effect_col else "NA"
                effect = sign_effect(raw_effect)
                by_contig[contig].append({"contig": contig, "start": start, "end": end, "dmr_id": dmr_id or f"dmr_{i}", "effect": effect or "NA"})
                parse_counts["parsed_tsv"] += 1

    for contig in by_contig:
        by_contig[contig].sort(key=lambda x: (x["start"], x["end"]))
    starts = {contig: [x["start"] for x in rows] for contig, rows in by_contig.items()}
    return by_contig, starts, parse_counts, detected


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def find_overlaps(contig, start, end, dmrs, starts, min_bp):
    rows = dmrs.get(contig, [])
    st = starts.get(contig, [])
    if not rows:
        return []
    stop = bisect.bisect_right(st, end)
    hits = []
    for dmr in rows[:stop]:
        if dmr["end"] < start:
            continue
        os = max(start, dmr["start"])
        oe = min(end, dmr["end"])
        bp = oe - os + 1
        if bp >= min_bp:
            hit = dict(dmr)
            hit["overlap_bp"] = bp
            hits.append(hit)
    return hits


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    dmrs, starts, parse_counts, detected = read_dmrs(args.dmrs, args)

    with open(args.candidates, "r", newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames or []
        rows = list(reader)

    for col in [args.candidate_id_column, args.contig_column, args.start_column, args.end_column]:
        if col not in headers:
            raise SystemExit(f"ERROR: missing column {col!r} in candidate TSV")

    extra = ["dmr_overlap_status", "dmr_overlap_count", "dmr_overlap_bp_unique", "dmr_overlap_fraction", "dmr_overlap_ids", "dmr_overlap_effects"]
    out_headers = headers + [x for x in extra if x not in headers]
    out_rows = []
    status_counts = Counter()
    effect_counts = Counter()

    for row in rows:
        contig = clean(row.get(args.contig_column))
        start = to_int(row.get(args.start_column))
        end = to_int(row.get(args.end_column))
        out = dict(row)
        if not contig or start is None or end is None:
            status = "no_candidate_coordinates"
            out.update({"dmr_overlap_status": status, "dmr_overlap_count": 0, "dmr_overlap_bp_unique": "NA", "dmr_overlap_fraction": "NA", "dmr_overlap_ids": "NA", "dmr_overlap_effects": "NA"})
            status_counts[status] += 1
            out_rows.append(out)
            continue
        if end < start:
            start, end = end, start
        interval_bp = end - start + 1
        hits = find_overlaps(contig, start, end, dmrs, starts, args.min_overlap_bp)
        clipped = merge_intervals([(max(start, h["start"]), min(end, h["end"])) for h in hits])
        unique_bp = sum(e - s + 1 for s, e in clipped)
        effects = Counter(h["effect"] for h in hits if clean(h.get("effect")))
        for effect, count in effects.items():
            effect_counts[effect] += count
        status = "dmr_overlap" if hits else "no_dmr_overlap"
        status_counts[status] += 1
        out.update({
            "dmr_overlap_status": status,
            "dmr_overlap_count": len(hits),
            "dmr_overlap_bp_unique": unique_bp,
            "dmr_overlap_fraction": f"{unique_bp / interval_bp:.6f}" if interval_bp > 0 else "NA",
            "dmr_overlap_ids": ";".join(h["dmr_id"] for h in hits[:50]) if hits else "NA",
            "dmr_overlap_effects": ";".join(f"{k}:{v}" for k, v in effects.most_common()) if effects else "NA",
        })
        out_rows.append(out)

    out_tsv = os.path.join(args.outdir, f"{args.prefix}.dmr_overlap.tsv")
    summary_txt = os.path.join(args.outdir, f"{args.prefix}.dmr_overlap.summary.txt")

    with open(out_tsv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=out_headers, extrasaction="ignore")
        writer.writeheader()
        for row in out_rows:
            writer.writerow(row)

    with open(summary_txt, "w", encoding="utf-8") as summary:
        summary.write("Dark candidate DMR-overlap summary\n")
        summary.write("==================================\n")
        summary.write(f"Candidate TSV: {args.candidates}\n")
        summary.write(f"DMR intervals: {args.dmrs}\n")
        summary.write(f"Output directory: {args.outdir}\n")
        summary.write(f"Candidate rows read: {len(rows)}\n")
        summary.write("\nDetected DMR columns/settings:\n")
        for key, value in detected.items():
            summary.write(f"- {key}: {value}\n")
        summary.write("\nDMR parse counts:\n")
        for key, value in parse_counts.most_common():
            summary.write(f"- {key}: {value}\n")
        summary.write(f"- DMR contigs indexed: {len(dmrs)}\n")
        summary.write("\nDMR-overlap status counts:\n")
        for key, value in status_counts.most_common():
            summary.write(f"- {key}: {value}\n")
        summary.write("\nDMR effect counts among overlaps:\n")
        for key, value in effect_counts.most_common(25):
            summary.write(f"- {key}: {value}\n")
        summary.write("\nOutputs:\n")
        summary.write(f"- Candidate DMR-overlap table: {out_tsv}\n")
        summary.write(f"- Summary: {summary_txt}\n")

    print("Done.")
    print(f"Summary: {summary_txt}")
    print(f"DMR-overlap table: {out_tsv}")


if __name__ == "__main__":
    main()
