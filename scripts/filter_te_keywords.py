"""Drop TE-related gene models by keyword in their functional annotation.

Used by workflow 06, step 9. Removes gene blocks whose attributes (Note=
descriptions, Dbxref) contain any configured keyword, case-insensitively,
e.g. 'transpos', 'gypsy', 'copia', 'reverse transcriptase'. The log lists
how many genes each keyword removed.
"""

import os
import sys

sys.path.insert(0, os.path.join(snakemake.scriptdir, "lib"))  # noqa: F821
from genome_utils import (count_genes, drop_genes, gene_blocks,  # noqa: E402
                          gene_id_of_block, load_gff, write_gff)

sample = snakemake.wildcards.sample  # noqa: F821
hap = snakemake.wildcards.hap  # noqa: F821
keywords = snakemake.params.keywords  # noqa: F821
if isinstance(keywords, str):
    keywords = [k.strip() for k in keywords.split(",") if k.strip()]
kws = [k.lower() for k in keywords]

header, rows = load_gff(snakemake.input.gff)  # noqa: F821
blocks = gene_blocks(rows)

te_genes, hits = set(), {}
for b in blocks:
    blob = " ".join(p[8] for p in b).lower()
    matched = [k for k in kws if k in blob]
    gid = gene_id_of_block(b)
    if matched and gid:
        te_genes.add(gid)
        for k in matched:
            hits[k] = hits.get(k, 0) + 1

kept_rows = drop_genes(rows, te_genes)
write_gff(snakemake.output.gff, header, kept_rows)  # noqa: F821
with open(snakemake.output.ids, "w") as out:  # noqa: F821
    out.write("\n".join(sorted(te_genes)) + "\n" if te_genes else "")

lines = [f"[{sample}\t{hap}] keyword TE filter: {len(blocks):,} -> "
         f"{count_genes(kept_rows):,} genes ({len(te_genes):,} removed)"]
for k in sorted(hits, key=lambda x: -hits[x]):
    lines.append(f"    {k:<32}{hits[k]:>8,} genes")
with open(snakemake.log[0], "w") as lg:  # noqa: F821
    lg.write("\n".join(lines) + "\n")
print("  " + lines[0])
