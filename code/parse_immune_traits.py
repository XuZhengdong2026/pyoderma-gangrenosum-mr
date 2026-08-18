# -*- coding: utf-8 -*-
"""Parse immune_traits.tsv (R-style c(...) vectors) into id/trait CSV."""
import csv
import re

SRC = r"F:\gwas_data\immune_traits.tsv"
OUT = r"F:\gwas_data\immune_traits_parsed.csv"

text = open(SRC, encoding="utf-8").read()

blocks = []
cur = None
for ln in text.splitlines():
    ln = ln.strip()
    if not ln:
        continue
    if ln.startswith("c("):
        if cur is not None:
            blocks.append(cur)
        cur = ln
    elif cur is not None:
        cur += " " + ln
if cur is not None:
    blocks.append(cur)


def parse(block):
    inner = block[block.index("(") + 1: block.rindex(")")]
    return re.findall(r'"((?:[^"\\]|\\.)*)"', inner)


ids = parse(blocks[0]) if len(blocks) > 0 else []
names = parse(blocks[1]) if len(blocks) > 1 else []
print("blocks:", len(blocks), "ids:", len(ids), "names:", len(names))

rows = list(zip(ids, names))
with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["id", "trait"])
    w.writerows(rows)

for kw in ["neutrophil", "monocyte", "regulatory", "CD4", "CD8",
           "dendritic", "myeloid", "B cell", "NK", "basophil", "eosinophil"]:
    hits = [r for r in rows if kw.lower() in r[1].lower()]
    print(kw, len(hits))
    for h in hits[:8]:
        print("   ", h)
