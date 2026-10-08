#!/usr/bin/env python3
"""
Convert the supplementary-tables workbook into JSON files for the
interactive viewer in docs/.

Usage:
    python scripts/build_supp_tables.py path/to/supplementary_tables.xlsx [docs]

Writes:
    docs/data/manifest.json   list of tables (id, title, group, size)
    docs/data/S<n>.json       one file per table (columns, rows, notes, glossary)
    docs/data/search.json     small cross-table index of text values
    docs/data/supplementary_tables.xlsx   copy of the workbook for download

Requires: openpyxl  (pip install openpyxl)

The workbook mixes layouts: one- and two-row headers, section rows such as
"Forward MR" or "Validation samples:", repeated header blocks, notes and
glossaries under the data. TABLES below records the layout of each sheet;
everything else is detected automatically. If you add or reorder sheets,
update TABLES and GROUPS and re-run.
"""
import datetime as dt
import json
import math
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import openpyxl

# Sidebar groups, in display order. Table IDs follow the sheet names
# ("Table S5 ..." -> "S5"). If two sheets carry the same number, the
# second and later ones get a letter suffix in sheet order (S19a, S19b).
GROUPS = [
    ("cohorts", "ENIGMA cohorts", ["S1", "S2", "S3", "S4", "S9"]),
    ("gwas", "BAG Han GWAS", ["S5", "S6", "S7", "S8", "S10", "S11", "S12", "S13"]),
    ("factor", "BAG factor", ["S14", "S15", "S16", "S17", "S18", "S19a", "S19b", "S20", "S21"]),
    ("mr", "Mendelian randomisation", ["S22", "S23", "S24", "S25"]),
    ("pgs", "Polygenic scores", ["S26", "S27", "S28", "S29", "S30", "S33"]),
    ("phewas", "PheWAS", ["S31", "S32"]),
]

# Per-table layout. Row numbers are 1-based Excel rows.
#   header:   row holding the column names
#   group:    optional row above it with spanning group labels
#   start_col: first used column (1 = A)
#   first_section: label for rows before the first section row
#   filldown: column indexes (0-based, after start_col) to fill downwards
#   notes_from: first row of notes when they lack a "Note." prefix
#   heatmap:  shade numeric cells ("seq" = sequential colour scale)
#   title:    override the title written in the sheet
TABLES = {
    "S1": dict(group=3, header=4, first_section="Discovery cohorts"),
    "S2": dict(group=3, header=4),
    "S3": dict(header=3, first_section="Discovery cohorts"),
    "S4": dict(group=3, header=4, first_section="Discovery cohorts"),
    "S5": dict(group=3, header=4),
    "S6": dict(header=3, notes_from=15),
    "S9": dict(group=3, header=4),
    "S14": dict(header=3, heatmap="seq"),
    "S26": dict(header=4, start_col=2),
    "S29": dict(group=3, header=4, filldown=[0]),
    "S32": dict(group=2, header=3),
}
DEFAULT = dict(header=3)

P_NAME = re.compile(
    r"^(p|pval|p-value|pvalue|p\.uncorrected|p\.fdr\.corrected|adjp|fdr|fdr_bh|"
    r"top p|top fdr|hetpva|q_pval|intercept_pval|p_difference|.*_p|"
    r"eqtlmapminp|eqtlmapminq|p\s*\(.*\))$", re.I)
NUM = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")
NA_VALUES = {"NA", "N/A", "na", "NaN", "–", "-", ""}


def clean_text(v):
    if v is None:
        return None
    if isinstance(v, (dt.datetime, dt.date)):
        return v.strftime("%Y-%m-%d")
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return None
        return int(v) if v.is_integer() and abs(v) < 1e15 else v
    if isinstance(v, (int, bool)):
        return int(v)
    s = str(v).replace("\xa0", " ").strip()
    if s == "":
        return None
    if NUM.match(s.replace(",", "")) and not (s.startswith("0") and len(s) > 1 and s[1].isdigit()):
        # Numeric text. Commas only count as thousand separators, e.g. "670,739".
        if "," in s and not re.match(r"^\d{1,3}(,\d{3})+$", s):
            return s
        f = float(s.replace(",", ""))
        return int(f) if f.is_integer() and "e" not in s.lower() and "." not in s else f
    return s


def header_name(v):
    if v is None:
        return ""
    return re.sub(r"\s+", " ", str(v).replace("\xa0", " ")).strip()


def nonempty(row):
    return [c for c in row if c is not None]


def build_table(ws, tid):
    cfg = {**DEFAULT, **TABLES.get(tid, {})}
    start = cfg.get("start_col", 1) - 1
    raw = [list(r) for r in ws.iter_rows(values_only=True)]
    raw = [[clean_text(c) for c in r[start:]] for r in raw]
    width = max((len(r) for r in raw), default=0)
    raw = [r + [None] * (width - len(r)) for r in raw]

    h = cfg["header"] - 1
    g = cfg.get("group")
    g = g - 1 if g else None

    title = cfg.get("title")
    if not title:
        for r in raw[:h]:
            ne = nonempty(r)
            if ne and str(ne[0]).lower().startswith("table s"):
                title = str(ne[0]).strip()
                break
    title = re.sub(r"\s+", " ", title or f"Table {tid}")

    sub = [header_name(c) for c in raw[h]]
    grp = [""] * width
    if g is not None:
        # Spread each group label over its merged range, or else rightwards
        # until the next label or a blank spacer column.
        merged = {}
        for m in ws.merged_cells.ranges:
            if m.min_row == g + 1:
                for c in range(m.min_col, m.max_col + 1):
                    merged[c - 1 - start] = header_name(ws.cell(m.min_row, m.min_col).value)
        cur = ""
        for i in range(width):
            label = header_name(raw[g][i])
            if label:
                cur = label
            elif i in merged:
                cur = merged[i]
            elif not sub[i] or (merged and i not in merged):
                cur = ""
            grp[i] = cur if (cur and (label or i in merged or not merged)) else ""
            if label and not sub[i]:
                # The group label is the column's own name (e.g. "Cohort").
                sub[i], grp[i] = label, ""
    header_tokens = {str(x).lower() for x in sub + grp if x}
    n_named = sum(1 for x in sub if x)

    rows, sections, notes, glossary = [], [], [], []
    section = cfg.get("first_section")
    mode = "data"
    prev_blank = True
    notes_from = cfg.get("notes_from")
    for rn, r in enumerate(raw[h + 1:], start=h + 2):
        ne = nonempty(r)
        if notes_from and rn >= notes_from:
            mode = "notes"
        if not ne:
            prev_blank = True
            continue
        first = str(ne[0])
        if mode == "data" and first.lower().startswith(("note", "annotations")):
            mode = "notes"
        if mode == "notes":
            if first.lower() == "annotations":
                prev_blank = False
                continue
            if len(ne) == 1:
                m = re.match(r"^([A-Za-z_][\w.]*)\s*=\s*(.+)$", first, re.S)
                if m:                       # "Term = definition"
                    glossary.append([m.group(1), m.group(2).strip()])
                elif not isinstance(ne[0], (int, float)) and len(first.strip(" .")) > 4:  # skip stray numbers / bare "Note."
                    notes.append(first)
            else:
                glossary.append([str(ne[0]), str(ne[-1])])
            continue
        # Repeated header block (e.g. before "Validation samples").
        texts = [str(x).lower() for x in ne]
        if len(ne) >= 2 and sum(t in header_tokens for t in texts) >= 0.5 * len(ne):
            prev_blank = False
            continue
        only_first = r[0] is not None and len(ne) == 1
        if only_first and (n_named >= 4 or prev_blank):
            section = first.rstrip(":").strip()
            if section.lower() == "validation samples":
                section = "Validation samples"
            prev_blank = False
            continue
        rows.append(r)
        sections.append(section)
        prev_blank = False

    for col in cfg.get("filldown", []):
        last = None
        for r in rows:
            if r[col] is None:
                r[col] = last
            else:
                last = r[col]

    # Drop columns empty in both header and data.
    keep = [i for i in range(width)
            if sub[i] or grp[i] or any(r[i] is not None for r in rows)]
    has_sections = any(s for s in sections) and len(set(sections)) > 1 or (
        any(sections) and cfg.get("first_section") is None)

    columns = []
    if has_sections:
        columns.append(dict(name="Section", group="", type="text", kind="section"))
    for i in keep:
        vals = [r[i] for r in rows if r[i] is not None and not (isinstance(r[i], str) and r[i] in NA_VALUES)]
        numeric = bool(vals) and sum(isinstance(v, (int, float)) for v in vals) >= 0.9 * len(vals)
        name = sub[i] or grp[i] or f"Column {i + 1}"
        kind = ""
        strs = [v for v in vals if isinstance(v, str)]
        if numeric and P_NAME.match(name):
            kind = "p"
        elif strs and sum(bool(re.match(r"^rs\d+$", v)) for v in strs) > 0.8 * len(strs):
            kind = "snp"
        elif strs and sum(bool(re.match(r"^ENSG\d+$", v)) for v in strs) > 0.8 * len(strs):
            kind = "ensg"
        elif name.upper() == "PMID":
            kind = "pmid"
        elif strs and sum(bool(re.match(r"^(https?://|www\.)\S+$", v)) for v in strs) > 0.5 * len(strs):
            kind = "url"
        elif name.lower() in ("symbol", "hugo", "nearest gene", "mappedgene", "reportedgene", "topgene"):
            kind = "gene"
        if numeric and kind == "" and name.lower() in ("pmid", "entrezid", "chr", "chromosome"):
            kind = "id"
        columns.append(dict(name=name, group=grp[i] if sub[i] else "",
                            type="num" if numeric else "text", kind=kind))

    out_rows = []
    for r, s in zip(rows, sections):
        vals = [r[i] for i in keep]
        out_rows.append(([s] if has_sections else []) + vals)

    heat = cfg.get("heatmap")
    return dict(
        id=tid, sheet=ws.title, title=title,
        caption=re.sub(r"^Table\s*S\d+[a-z]?\.\s*", "", title),
        columns=columns, rows=out_rows, notes=notes, glossary=glossary,
        heatmap=heat,
    )


def search_tokens(table):
    """Distinct short text values per table, for cross-table search."""
    seen = defaultdict(int)
    for row in table["rows"]:
        for col, v in zip(table["columns"], row):
            if not isinstance(v, str) or col["kind"] in ("url",) or len(v) > 90:
                continue
            parts = [v]
            # Gene/SNP lists such as "MAPT, CRHR1", "ENSG1:ENSG2" or "MIR6500 - C1orf185".
            split = [s for s in re.split(r"\s*[,:;]\s*|\s+-\s+", v) if s]
            if len(split) > 1 and all(" " not in s for s in split):
                parts = split
            for p in set(s.strip() for s in parts):
                if p and p not in NA_VALUES and not NUM.match(p):
                    seen[p] += 1
    return seen


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    xlsx = Path(sys.argv[1])
    docs = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("docs")
    out = docs / "data"
    out.mkdir(parents=True, exist_ok=True)

    wb = openpyxl.load_workbook(xlsx, data_only=True)  # cached values, not formulas
    sheets = []                                   # (id, worksheet) in sheet order
    nums = [re.match(r"Table\s*S(\d+)", ws.title) for ws in wb.worksheets]
    counts = defaultdict(int)
    for m in nums:
        if m:
            counts[m.group(1)] += 1
    seen = defaultdict(int)
    for ws, m in zip(wb.worksheets, nums):
        if not m:
            print(f"skipping sheet without a table number: {ws.title}")
            continue
        n = m.group(1)
        tid = f"S{n}" if counts[n] == 1 else f"S{n}{'abcdefgh'[seen[n]]}"
        seen[n] += 1
        if counts[n] > 1:
            print(f"warning: more than one sheet is numbered S{n}; {ws.title!r} -> {tid}")
        sheets.append((tid, ws))

    group_of = {i: key for key, _, ids in GROUPS for i in ids}
    manifest, index = [], defaultdict(dict)
    for tid, ws in sheets:
        t = build_table(ws, tid)
        if tid not in group_of:
            print(f"warning: {tid} is not in GROUPS; listed under 'Other'")
        t["group"] = group_of.get(tid, "other")
        with open(out / f"{t['id']}.json", "w", encoding="utf-8") as f:
            json.dump(t, f, ensure_ascii=False, separators=(",", ":"))
        manifest.append(dict(id=t["id"], title=t["title"], caption=t["caption"],
                             group=t["group"], rows=len(t["rows"]), cols=len(t["columns"])))
        for token, n in search_tokens(t).items():
            index[token][t["id"]] = n
        print(f"{t['id']:>4}  {len(t['rows']):>5} rows × {len(t['columns']):>2} cols  {t['caption'][:70]}")

    groups = [dict(key=k, label=l) for k, l, _ in GROUPS] + [dict(key="other", label="Other")]
    with open(out / "manifest.json", "w", encoding="utf-8") as f:
        json.dump(dict(groups=groups, tables=manifest,
                       built=dt.date.today().isoformat(), source=xlsx.name),
                  f, ensure_ascii=False, indent=1)
    with open(out / "search.json", "w", encoding="utf-8") as f:
        json.dump(index, f, ensure_ascii=False, separators=(",", ":"))
    shutil.copyfile(xlsx, out / "supplementary_tables.xlsx")


if __name__ == "__main__":
    main()
