"""Gene counts at each filtering stage for one haplotype.

Used by workflow 06, step 12a. Reports the number of genes after Helixer,
after the minimum-length filter, after the Pfam TE filter and after the
keyword TE filter, each as a percentage of the raw count.
"""

import os
import sys

sys.path.insert(0, os.path.join(snakemake.scriptdir, "lib"))  # noqa: F821
from genome_utils import count_genes, load_gff  # noqa: E402

sample = snakemake.wildcards.sample  # noqa: F821
hap = snakemake.wildcards.hap  # noqa: F821
inp = snakemake.input  # noqa: F821
min_aa = snakemake.params.min_aa  # noqa: F821

stages = [("raw_helixer", inp.raw),
          (f"after_min{min_aa}aa", inp.aa),
          ("after_pfam_TE", inp.pfam),
          ("after_keyword_TE", inp.final)]
counts = [(name, count_genes(load_gff(path)[1])) for name, path in stages]
start = counts[0][1] or 1
with open(snakemake.output.tsv, "w") as out:  # noqa: F821
    out.write("sample\thaplotype\tstage\tn_genes\tpct_of_raw\n")
    for name, n in counts:
        out.write(f"{sample}\t{hap}\t{name}\t{n}\t{100 * n / start:.1f}\n")
