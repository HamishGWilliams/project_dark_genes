#!/usr/bin/env python3
"""Inspect an .xlsx workbook containing differential-expression results.

This version uses only the Python standard library. It can run on the cluster even
when pandas/openpyxl are unavailable. It reports sheet names, detected columns,
and a small row preview for each sheet.

Supported input: .xlsx. Legacy .xls files are not supported without external
libraries.
"""

import argparse
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "pkgrel": "http://schemas.openxmlformats.org/package/2006/relationships",
}


def parse_args():
    p = argparse.ArgumentParser(description="Inspect sheets and columns in a DE .xlsx workbook without pandas.")
    p.add_argument("--excel", required=True, help="Path to .xlsx workbook")
    p.add_argument("--out", default=None, help="Optional text report path")
    p.add_argument("--head", type=int, default=3, help="Data rows to preview per sheet")
    return p.parse_args()


def col_index_from_cell_ref(ref):
    letters = re.sub(r"[^A-Z]", "", ref.upper())
    idx = 0
    for char in letters:
        idx = idx * 26 + (ord(char) - ord("A") + 1)
    return idx - 1


def read_xml(zf, path):
    with zf.open(path) as handle:
        return ET.parse(handle).getroot()


def read_shared_strings(zf):
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []
    root = read_xml(zf, "xl/sharedStrings.xml")
    strings = []
    for si in root.findall("main:si", NS):
        texts = []
        for t in si.findall(".//main:t", NS):
            texts.append(t.text or "")
        strings.append("".join(texts))
    return strings


def sheet_paths(zf):
    workbook = read_xml(zf, "xl/workbook.xml")
    rels = read_xml(zf, "xl/_rels/workbook.xml.rels")

    rel_map = {}
    for rel in rels.findall("pkgrel:Relationship", NS):
        rid = rel.attrib.get("Id")
        target = rel.attrib.get("Target")
        if rid and target:
            if not target.startswith("xl/"):
                target = "xl/" + target.lstrip("/")
            rel_map[rid] = target

    sheets = []
    for sheet in workbook.findall("main:sheets/main:sheet", NS):
        name = sheet.attrib.get("name", "unknown_sheet")
        rid = sheet.attrib.get("{%s}id" % NS["rel"])
        path = rel_map.get(rid)
        if path:
            sheets.append((name, path))
    return sheets


def cell_value(cell, shared_strings):
    cell_type = cell.attrib.get("t")
    value_node = cell.find("main:v", NS)
    inline_node = cell.find("main:is/main:t", NS)

    if cell_type == "inlineStr":
        return inline_node.text if inline_node is not None and inline_node.text is not None else ""

    if value_node is None or value_node.text is None:
        return ""

    raw = value_node.text
    if cell_type == "s":
        try:
            return shared_strings[int(raw)]
        except Exception:
            return raw
    return raw


def read_sheet_preview(zf, sheet_path, shared_strings, data_rows=3):
    root = read_xml(zf, sheet_path)
    rows = []
    max_cols = 0

    for row in root.findall("main:sheetData/main:row", NS):
        values = []
        for cell in row.findall("main:c", NS):
            ref = cell.attrib.get("r", "")
            col_idx = col_index_from_cell_ref(ref) if ref else len(values)
            while len(values) <= col_idx:
                values.append("")
            values[col_idx] = cell_value(cell, shared_strings)
        max_cols = max(max_cols, len(values))
        rows.append(values)
        if len(rows) >= data_rows + 1:
            break

    for row in rows:
        while len(row) < max_cols:
            row.append("")

    header = rows[0] if rows else []
    preview = rows[1:] if len(rows) > 1 else []
    return header, preview


def format_preview(header, preview):
    if not header and not preview:
        return "<empty sheet>"
    rows = [header] + preview
    widths = []
    ncols = max(len(r) for r in rows) if rows else 0
    for i in range(ncols):
        widths.append(min(40, max(len(str(r[i])) if i < len(r) else 0 for r in rows)))
    lines = []
    for row in rows:
        cells = []
        for i in range(ncols):
            value = str(row[i]) if i < len(row) else ""
            if len(value) > 40:
                value = value[:37] + "..."
            cells.append(value.ljust(widths[i]))
        lines.append(" | ".join(cells).rstrip())
    return "\n".join(lines)


def main():
    args = parse_args()
    if not os.path.exists(args.excel):
        raise SystemExit(f"ERROR: workbook not found: {args.excel}")
    if not args.excel.lower().endswith(".xlsx"):
        raise SystemExit("ERROR: this dependency-free inspector supports .xlsx files only. Convert .xls to .xlsx first.")

    lines = []
    with zipfile.ZipFile(args.excel) as zf:
        shared_strings = read_shared_strings(zf)
        sheets = sheet_paths(zf)
        lines.append("Differential-expression workbook inspection")
        lines.append("===========================================")
        lines.append(f"Workbook: {args.excel}")
        lines.append(f"Sheets: {len(sheets)}")
        lines.append("")

        for name, path in sheets:
            header, preview = read_sheet_preview(zf, path, shared_strings, data_rows=max(args.head, 0))
            lines.append(f"## Sheet: {name}")
            lines.append(f"Path: {path}")
            lines.append(f"Detected columns ({len(header)}):")
            for col in header:
                lines.append(f"- {col}")
            if args.head > 0:
                lines.append("Preview:")
                lines.append(format_preview(header, preview))
            lines.append("")

    report = "\n".join(lines)
    if args.out:
        outdir = os.path.dirname(args.out)
        if outdir:
            os.makedirs(outdir, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as handle:
            handle.write(report + "\n")
        print(f"Wrote report: {args.out}")
    else:
        print(report)


if __name__ == "__main__":
    main()
