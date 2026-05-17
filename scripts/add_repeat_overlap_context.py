#!/usr/bin/env python3
"""Add repeat/TE overlap context to dark candidates.

Inputs:
  --candidates : TSV with contig and genomic coordinate columns
  --repeats    : GFF3 or BED repeat annotation

Outputs:
  <prefix>.repeat_overlap.tsv
  <prefix>.repeat_overlap.summary.txt

Candidate and GFF3 coordinates are treated as 1-based closed intervals.
BED intervals are converted from 0-based half-open to 1-based closed.
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
    p = argparse.ArgumentParser(description="Annotate dark candidates with repeat-overlap context.")
    p.add_argument("--candidates", required=True)
    p.add_argument("--repeats", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--prefix", default="equina_dark_candidates")
    p.add_argument("--candidate-id-column", default="protein_id")
    p.add_argument("--contig-column", default="contig")
    p.add_argument("--start-column", default="transcript_start")
    p.add_argument("--end-column", default="transcript_end")
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


def parse_attrs(attr):
    out = {}
    for part in str(attr or "").split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
        elif " " in part:
            k, v = part.split(" ", 1)
        else:
            continue
        out[k.strip()] = v.strip().strip('"')
    return out


def get_repeat_label(attrs, fallback):
    for key in ["Name", "Target", "ID", "Class", "classification", "family", "repeat_class"]:
        value = clean(attrs.get(key, ""))
        if value:
            return value.replace("#", ":")
    return fallback


def parse_repeat_line(line):
    if not line.strip() or line.startswith("#"):
        return None
    parts = line.rstrip("\n").split("\t")
    if len(parts) == 9:
        contig, source, ftype, start, end, score, strand, phase, attrs_raw = parts
        start_i, end_i = to_int(start), to_int(end)
        if start_i is None or end_i is None:
            return None
        if end_i < start_i:
            start_i, end_i = end_i, start_i
        attrs = parse_attrs(attrs_raw)
        label = get_repeat_label(attrs, ftype)
        return {
            "contig": contig,
            "start": start_i,
            "end": end_i,
            "label": label,
            "feature_type": ftype,
            "source": source,
        }
    if len(parts) >= 3:
        start0, end0 = to_int(parts[1]), to_int(parts[2])
        if start0 is None or end0 is None:
            return None
        label = parts[3] if len(parts) >= 4 else "repeat_interval"
        return {
            "contig": parts[0],
            "start": start0 + 1,
            "end": end0,
            "label": label,
            "feature_type": "BED_repeat_interval",
            "source": "BED",
        }
    return None


def read_repeats(path):
    by_contig = defaultdict(list)
    parsed = 0
    skipped = 0
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            repeat = parse_repeat_line(line)
            if repeat is None:
                skipped += 1
                continue
            by_contig[repeat["contig"]].append(repeat)
            parsed += 1
    for contig in by_contig:
        by_contig[contig].sort(key=lambda r: (r["start"], r["end"]))
    starts = {contig: [r["start"] for r in rows] for contig, rows in by_contig.items()}
    return by_contig, starts, parsed, skipped


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


def find_overlaps(contig, start, end, repeats, starts, min_bp):
    rows = repeats.get(contig, [])
    st = starts.get(contig, [])
    if not rows:
        return []
    stop = bisect.bisect_right(st, end)
    hits = []
    for rep in rows[:stop]:
        if rep["end"] < start:
            continue
        os = max(start, rep["start"])
        oe = min(end, rep["end"])
        bp = oe - os + 1
        if bp >= min_bp:
            hit = dict(rep)
            hit["overlap_bp"] = bp
            hits.append(hit)
    return hits


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    repeats, starts, n_repeats, n_skipped = read_repeats(args.repeats)

    with open(args.candidates, "r", newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames or []
        rows = list(reader)

    for col in [args.candidate_id_column, args.contig_column, args.start_column, args.end_column]:
        if col not in headers:
            raise SystemExit(f"ERROR: missing column {col!r} in candidate TSV")

    extra = [
        "repeat_overlap_status",
        "repeat_overlap_count",
        "repeat_overlap_bp_unique",
        "repeat_overlap_fraction",
        "repeat_overlap_labels",
    ]
    out_headers = headers + [x for x in extra if x not in headers]
    out_rows = []
    status_counts = Counter()
    label_counts = Counter()

    for row in rows:
        contig = clean(row.get(args.contig_column))
        start = to_int(row.get(args.start_column))
        end = to_int(row.get(args.end_column))
        out = dict(row)
        if not contig or start is None or end is None:
            status = "no_candidate_coordinates"
            out.update({
                "repeat_overlap_status": status,
                "repeat_overlap_count": 0,
                "repeat_overlap_bp_unique": "NA",
                "repeat_overlap_fraction": "NA",
                "repeat_overlap_labels": "NA",
            })
            status_counts[status] += 1
            out_rows.append(out)
            continue
        if end < start:
            start, end = end, start
        interval_bp = end - start + 1
        hits = find_overlaps(contig, start, end, repeats, starts, args.min_overlap_bp)
        clipped = merge_intervals([(max(start, h["start"]), min(end, h["end"])) for h in hits])
        unique_bp = sum(e - s + 1 for s, e in clipped)
        labels = Counter(h["label"] for h in hits)
        for label, count in labels.items():
            label_counts[label] += count
        status = "repeat_overlap" if hits else "no_repeat_overlap"
        status_counts[status] += 1
        out.update({
            "repeat_overlap_status": status,
            "repeat_overlap_count": len(hits),
            "repeat_overlap_bp_unique": unique_bp,
            "repeat_overlap_fraction": f"{unique_bp / interval_bp:.6f}" if interval_bp > 0 else "NA",
            "repeat_overlap_labels": ";".join(f"{k}:{v}" for k, v in labels.most_common()) if labels else "NA",
        })
        out_rows.append(out)

    out_tsv = os.path.join(args.outdir, f"{args.prefix}.repeat_overlap.tsv")
    summary_txt = os.path.join(args.outdir, f"{args.prefix}.repeat_overlap.summary.txt")

    with open(out_tsv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=out_headers, extrasaction="ignore")
        writer.writeheader()
        for row in out_rows:
            writer.writerow(row)

    with open(summary_txt, "w", encoding="utf-8") as summary:
        summary.write("Dark candidate repeat-overlap summary\n")
        summary.write("=====================================\n")
        summary.write(f"Candidate TSV: {args.candidates}\n")
        summary.write(f"Repeat annotation: {args.repeats}\n")
        summary.write(f"Output directory: {args.outdir}\n")
        summary.write(f"Candidate rows read: {len(rows)}\n")
        summary.write(f"Repeat records parsed: {n_repeats}\n")
        summary.write(f"Repeat records skipped: {n_skipped}\n")
        summary.write("\nRepeat-overlap status counts:\n")
        for key, value in status_counts.most_common():
            summary.write(f"- {key}: {value}\n")
        summary.write("\nTop overlapping repeat labels:\n")
        for key, value in label_counts.most_common(25):
            summary.write(f"- {key}: {value}\n")
        summary.write("\nOutputs:\n")
        summary.write(f"- Candidate repeat-overlap table: {out_tsv}\n")
        summary.write(f"- Summary: {summary_txt}\n")

    print("Done.")
    print(f"Summary: {summary_txt}")
    print(f"Repeat-overlap table: {out_tsv}")


if __name__ == "__main__":
    main()
