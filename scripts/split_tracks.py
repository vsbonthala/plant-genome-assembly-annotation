"""Split the final gene models into three complementary tracks.

Used by workflow 06, step 11.

    chrom : Chr01..ChrNN            placed chromosomes
    utg   : Chr0N_UA + unitigs      unanchored bins + unplaced sequence
    all   : everything (all = chrom + utg)

Splitting happens after annotation and filtering, so BLAST/MMseqs2 and
InterProScan (the most expensive steps) run once on 'all' rather than once
per track.
"""

import os
import sys

sys.path.insert(0, os.path.join(snakemake.scriptdir, "lib"))  # noqa: F821
from genome_utils import (gene_blocks, gff_attr, load_gff, read_fasta,  # noqa: E402
                          seq_track, seqid_sort_key, write_fasta, write_gff)

sample = snakemake.wildcards.sample  # noqa: F821
hap = snakemake.wildcards.hap  # noqa: F821
track = snakemake.wildcards.track  # noqa: F821
inp, outp = snakemake.input, snakemake.output  # noqa: F821

header, rows = load_gff(inp.gff)
blocks = gene_blocks(rows)
keep = blocks if track == "all" else [b for b in blocks if seq_track(b[0][0]) == track]
keep.sort(key=lambda b: (seqid_sort_key(b[0][0]), int(b[0][3])))
kept_rows = [p for b in keep for p in b]
write_gff(outp.gff, header, kept_rows)

keep_mrna = {gff_attr(p[8], "ID") for p in kept_rows if p[2] in ("mRNA", "transcript")}
keep_mrna.discard(None)
n_pep = write_fasta(outp.pep, read_fasta(inp.pep), keep=keep_mrna)
n_trn = write_fasta(outp.trn, read_fasta(inp.trn), keep=keep_mrna)

seqs = sorted({b[0][0] for b in keep}, key=seqid_sort_key)
msg = (f"[{sample}\t{hap}\t{track}] {len(keep):,} genes on {len(seqs)} sequences, "
       f"{n_pep:,} proteins, {n_trn:,} transcripts")
with open(snakemake.log[0], "w") as lg:  # noqa: F821
    lg.write(msg + "\n")
    lg.write("sequence order: " + " ".join(seqs[:20]) + (" ..." if len(seqs) > 20 else "") + "\n")
print("  " + msg)
