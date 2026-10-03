"""Sort a normalised GFF3 gene-by-gene into the 'all' track order.

Used by workflow 06, step 1. Gene blocks (gene + its children) are kept
intact and ordered Chr01..ChrNN, then Chr0N_UA bins, then unplaced unitigs,
and by start coordinate within each sequence.
"""

import os
import sys

sys.path.insert(0, os.path.join(snakemake.scriptdir, "lib"))  # noqa: F821
from genome_utils import gene_blocks, load_gff, seqid_sort_key, write_gff  # noqa: E402

header, rows = load_gff(snakemake.input.gff)  # noqa: F821
blocks = gene_blocks(rows)
blocks.sort(key=lambda b: (seqid_sort_key(b[0][0]), int(b[0][3])))
write_gff(snakemake.output.gff, header, [p for b in blocks for p in b])  # noqa: F821
print(f"  [{snakemake.wildcards.sample}\t{snakemake.wildcards.hap}\tall] "  # noqa: F821
      f"{len(blocks):,} genes, sorted Chr01..ChrNN -> UA -> utg")
