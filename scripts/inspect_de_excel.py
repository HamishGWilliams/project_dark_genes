#!/usr/bin/env python3
"""Inspect an Excel workbook containing differential-expression results.

This utility is intentionally lightweight: it reports sheet names, dimensions,
column names, and the first few rows so the DE integration step can be configured
without opening the workbook interactively on the cluster.
"""

import argparse
import os
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Inspect sheets and columns in a DE Excel workbook.")
    p.add_argument("--excel", required=True, help="Path to .xlsx/.xls workbook")
    p.add_argument("--out", default=None, help="Optional text report path")
    p.add_argument("--head", type=int, default=3, help="Rows to preview per sheet")
    return p.parse_args()


def main():
    args = parse_args()
    if not os.path.exists(args.excel):
        raise SystemExit(f"ERROR: workbook not found: {args.excel}")

    xls = pd.ExcelFile(args.excel)
    lines = []
    lines.append("Differential-expression workbook inspection")
    lines.append("===========================================")
    lines.append(f"Workbook: {args.excel}")
    lines.append(f"Sheets: {len(xls.sheet_names)}")
    lines.append("")

    for sheet in xls.sheet_names:
        df = pd.read_excel(args.excel, sheet_name=sheet, nrows=max(args.head, 1))
        lines.append(f"## Sheet: {sheet}")
        lines.append(f"Preview rows read: {len(df)}")
        lines.append(f"Columns ({len(df.columns)}):")
        for col in df.columns:
            lines.append(f"- {col}")
        if args.head > 0:
            lines.append("Preview:")
            lines.append(df.head(args.head).to_string(index=False))
        lines.append("")

    report = "\n".join(lines)
    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as handle:
            handle.write(report + "\n")
        print(f"Wrote report: {args.out}")
    else:
        print(report)


if __name__ == "__main__":
    main()
