"""Drop TE-related gene models that carry a TE-associated Pfam domain.

Used by workflow 06, step 8. Reads a list of Pfam IDs (one per line) and
removes every gene block whose GFF3 attributes mention any of them (the
IDs were written by merge_interproscan_gff.py). If no list is configured,
the step passes the GFF3 through unchanged and says so in the log.
"""

import os
import sys

sys.path.insert(0, os.path.join(snakemake.scriptdir, "lib"))  # noqa: F821
from genome_utils import (count_genes, drop_genes, gene_blocks,  # noqa: E402
                          gene_id_of_block, load_gff, write_gff)

sample = snakemake.wildcards.sample  # noqa: F821
hap = snakemake.wildcards.hap  # noqa: F821
pfam_list = snakemake.params.pfam_list  # noqa: F821

header, rows = load_gff(snakemake.input.gff)  # noqa: F821
blocks = gene_blocks(rows)

if not pfam_list or not os.path.exists(pfam_list):
    msg = f"[{sample}\t{hap}] Pfam TE filter SKIPPED (te_pfam_list not set or missing)"
    write_gff(snakemake.output.gff, header, rows)  # noqa: F821
    open(snakemake.output.ids, "w").close()  # noqa: F821
else:
    with open(pfam_list) as fh:
        pfam_ids = {line.strip().split()[0] for line in fh
                    if line.strip() and not line.startswith("#")}
    te_genes = set()
    for b in blocks:
        blob = " ".join(p[8] for p in b)
        if any(pid in blob for pid in pfam_ids):
            gid = gene_id_of_block(b)
            if gid:
                te_genes.add(gid)
    kept_rows = drop_genes(rows, te_genes)
    write_gff(snakemake.output.gff, header, kept_rows)  # noqa: F821
    with open(snakemake.output.ids, "w") as out:  # noqa: F821
        out.write("\n".join(sorted(te_genes)) + "\n" if te_genes else "")
    msg = (f"[{sample}\t{hap}] Pfam TE filter ({len(pfam_ids):,} IDs): "
           f"{len(blocks):,} -> {count_genes(kept_rows):,} genes "
           f"({len(te_genes):,} removed)")

with open(snakemake.log[0], "w") as lg:  # noqa: F821
    lg.write(msg + "\n")
print("  " + msg)
