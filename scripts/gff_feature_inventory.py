"""Inventory of the raw Helixer GFF3 before any processing.

Used by workflow 06, step 0. Reports feature types, genes per track
(chromosomes vs unplaced), CDS phase values, strands, and whether UTRs,
explicit introns and CDS phases are present - the things that decide
whether gffread will translate the models correctly.
"""

import os
import sys
from collections import Counter

sys.path.insert(0, os.path.join(snakemake.scriptdir, "lib"))  # noqa: F821
from genome_utils import load_gff, seq_track  # noqa: E402

sample = snakemake.wildcards.sample  # noqa: F821
hap = snakemake.wildcards.hap  # noqa: F821
gff_path = snakemake.input.gff  # noqa: F821

_, rows = load_gff(gff_path)
types = Counter(p[2] for p in rows)
phases = Counter(p[7] for p in rows if p[2] == "CDS")
strands = Counter(p[6] for p in rows)
tracks = Counter(seq_track(p[0]) for p in rows if p[2] == "gene")

with open(snakemake.output.inv, "w") as out:  # noqa: F821
    out.write(f"# sample:    {sample}\n# haplotype: {hap}\n# source:    {gff_path}\n\n")
    out.write("feature types:\n")
    for t, n in types.most_common():
        out.write(f"  {t:<20}{n:>10,}\n")
    out.write("\ngenes per track:\n")
    for t in ("chrom", "utg"):
        out.write(f"  {t:<20}{tracks.get(t, 0):>10,}\n")
    out.write("\nCDS phase values ('.' means phase NOT set):\n")
    for ph, n in phases.most_common():
        out.write(f"  {ph:<20}{n:>10,}\n")
    out.write("\nstrands:\n")
    for s, n in strands.most_common():
        out.write(f"  {s:<20}{n:>10,}\n")
    out.write("\nnotes:\n")
    out.write("  UTR features present : "
              f"{'yes' if any('utr' in t.lower() for t in types) else 'NO'}\n")
    out.write("  explicit introns     : "
              f"{'yes' if 'intron' in types else 'NO (inferred from exon gaps - GFF3-legal)'}\n")
    out.write("  CDS phase set        : "
              f"{'NO - gffread may mis-translate' if set(phases) <= {'.'} else 'yes'}\n")
