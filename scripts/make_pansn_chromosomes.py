"""Write a chromosome-only FASTA in PanSN-spec naming for pangenome graphs.

Used by workflow 05 (finishing). Keeps only the placed chromosomes
(Chr01..ChrNN; no Chr0N_UA bins or unplaced unitigs, which would create
spurious graph components) and renames them to PanSN-spec
(https://github.com/pangenome/PanSN-spec):

    [sample]#[haplotype_id]#[contig]      e.g.  SampleA#1#Chr01

PanSN requires an integer haplotype ID, taken from the config
(``pansn_hap_ids``). Compression and indexing for pggb happen in the next
rule (bgzip + samtools faidx).
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(snakemake.scriptdir, "lib"))  # noqa: F821
from genome_utils import read_fasta, write_fasta  # noqa: E402

sample = snakemake.wildcards.sample  # noqa: F821
hap = snakemake.wildcards.hap  # noqa: F821
hap_id = snakemake.params.hap_id  # noqa: F821
delim = snakemake.params.delim  # noqa: F821
n_chrom = int(snakemake.params.n_chrom)  # noqa: F821

prefix = f"{sample}{delim}{hap_id}{delim}"
chrom_re = re.compile(r"^Chr(\d+)$")

records = {k: v for k, v in read_fasta(snakemake.input.sorted_fa).items()  # noqa: F821
           if chrom_re.match(k)}
order = sorted(records, key=lambda k: int(chrom_re.match(k).group(1)))
renamed = {f"{prefix}{k}": records[k] for k in order}
write_fasta(snakemake.output.pansn, renamed)  # noqa: F821

total_bp = sum(len(v) for v in records.values())
msgs = [f"[{sample}] {hap}: PanSN prefix '{prefix}'",
        f"  chromosomes written: {len(records)} ({', '.join(order)})",
        f"  total bp: {total_bp:,}"]
if len(records) != n_chrom:
    msgs.append(f"  *** WARNING: expected {n_chrom} chromosomes, found "
                f"{len(records)} - pggb needs complete sets ***")
with open(snakemake.log[0], "w") as lg:  # noqa: F821
    lg.write("\n".join(msgs) + "\n")
print("\n".join(msgs))
