#!/usr/bin/env python3

import argparse
import csv
import os
import sys

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
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    if args.input in {"", "NA"} or not os.path.exists(args.input) or os.path.getsize(args.input) == 0:
        print(f"SKIP\t{args.input}\tmissing")
        return

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    long_rows = 0
    short_rows = 0
    total_rows = 0

    with open(args.input, "r", newline="", encoding="utf-8", errors="replace") as fin, \
            open(args.output, "w", newline="", encoding="utf-8") as fout:

        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(
            fout,
            delimiter="\t",
            quotechar='"',
            quoting=csv.QUOTE_MINIMAL,
            lineterminator="\n",
        )

        try:
            header = next(reader)
        except StopIteration:
            print(f"SKIP\t{args.input}\tempty")
            return

        ncol = len(header)
        writer.writerow(header)

        for row in reader:
            total_rows += 1

            if len(row) < ncol:
                row = row + [""] * (ncol - len(row))
                short_rows += 1

            elif len(row) > ncol:
                # Preserve early columns and collapse any excess fields into the final column.
                # This handles annotation text fields that contain extra literal tab characters.
                row = row[:ncol - 1] + [" | ".join(row[ncol - 1:])]
                long_rows += 1

            writer.writerow(row)

    print(
        f"OK\t{args.input}\t{args.output}\trows={total_rows}\tlong_rows_fixed={long_rows}\tshort_rows_fixed={short_rows}"
    )


if __name__ == "__main__":
    main()
