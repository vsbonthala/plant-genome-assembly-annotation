"""Drop gene models whose proteins are all shorter than ``min_aa``.

Used by workflow 06, step 7. A gene is kept if ANY of its transcripts
encodes a protein of at least ``min_aa`` residues (terminal '*' ignored).
mRNA -> gene links come from the GFF3 Parent attribute, and IDs are matched
exactly (never by substring).
"""

import os
import sys

sys.path.insert(0, os.path.join(snakemake.scriptdir, "lib"))  # noqa: F821
from genome_utils import (count_genes, drop_genes, gene_blocks,  # noqa: E402
                          gene_id_of_block, load_gff, mrna_to_gene, read_fasta,
                          write_gff)

sample = snakemake.wildcards.sample  # noqa: F821
hap = snakemake.wildcards.hap  # noqa: F821
min_aa = int(snakemake.params.min_aa)  # noqa: F821

prot = read_fasta(snakemake.input.pep)  # noqa: F821
keep_mrna = {n for n, s in prot.items() if len(s.rstrip("*")) >= min_aa}
with open(snakemake.output.ids, "w") as out:  # noqa: F821
    out.write("\n".join(sorted(keep_mrna)) + "\n")

header, rows = load_gff(snakemake.input.gff)  # noqa: F821
m2g = mrna_to_gene(rows)
keep_genes = {m2g[m] for m in keep_mrna if m in m2g}
all_genes = {gene_id_of_block(b) for b in gene_blocks(rows)}
drop = all_genes - keep_genes
kept_rows = drop_genes(rows, drop)
write_gff(snakemake.output.gff, header, kept_rows)  # noqa: F821

msg = (f"[{sample}\t{hap}] < {min_aa} aa filter: {len(all_genes):,} -> "
       f"{count_genes(kept_rows):,} genes ({len(drop):,} removed); "
       f"proteins kept {len(keep_mrna):,}/{len(prot):,}")
with open(snakemake.log[0], "w") as lg:  # noqa: F821
    lg.write(msg + "\n")
print("  " + msg)
