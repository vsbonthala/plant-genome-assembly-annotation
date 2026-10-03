"""Write the final non-TE protein and transcript sets.

Used by workflow 06, step 10. Keeps exactly the mRNAs that survive all
filters in the final GFF3, and writes their proteins, transcripts and an
ID list.
"""

import os
import sys

sys.path.insert(0, os.path.join(snakemake.scriptdir, "lib"))  # noqa: F821
from genome_utils import gff_attr, load_gff, read_fasta, write_fasta  # noqa: E402

inp, outp = snakemake.input, snakemake.output  # noqa: F821

_, rows = load_gff(inp.gff)
keep_mrna = {gff_attr(p[8], "ID") for p in rows if p[2] in ("mRNA", "transcript")}
keep_mrna.discard(None)
with open(outp.ids, "w") as out:
    out.write("\n".join(sorted(keep_mrna)) + "\n")

n_pep = write_fasta(outp.pep, read_fasta(inp.pep), keep=keep_mrna)
n_trn = write_fasta(outp.trn, read_fasta(inp.trn), keep=keep_mrna)
print(f"  [{snakemake.wildcards.sample}\t{snakemake.wildcards.hap}] final non-TE: "  # noqa: F821
      f"{len(keep_mrna):,} mRNAs -> {n_pep:,} proteins, {n_trn:,} transcripts")
