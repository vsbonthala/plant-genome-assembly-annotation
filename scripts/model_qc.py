"""Quality metrics for the final gene models of one haplotype and track.

Used by workflow 06, step 12b. Reports gene, mRNA, exon and CDS counts,
exons per mRNA, the mono-exonic fraction, protein length statistics,
internal stop codons, start-with-M and end-with-stop fractions, and domain
coverage (proteins with a Pfam or any InterPro signature). Low Pfam
coverage is expected on the 'utg' track (fragmentary sequence) but would
point to a model problem on 'chrom'.
"""

import os
import statistics as st
import sys

sys.path.insert(0, os.path.join(snakemake.scriptdir, "lib"))  # noqa: F821
from genome_utils import count_genes, gff_attr, load_gff, parse_ipscan, read_fasta  # noqa: E402

w = snakemake.wildcards  # noqa: F821
inp = snakemake.input  # noqa: F821

_, rows = load_gff(inp.gff)
n_gene = count_genes(rows)
n_mrna = sum(1 for p in rows if p[2] in ("mRNA", "transcript"))
n_exon = sum(1 for p in rows if p[2] == "exon")
n_cds = sum(1 for p in rows if p[2] == "CDS")
n_seqs = len({p[0] for p in rows})

exons_per = {}
for p in rows:
    if p[2] == "exon":
        par = gff_attr(p[8], "Parent")
        if par:
            exons_per[par] = exons_per.get(par, 0) + 1
nan = float("nan")
mono_frac = (sum(1 for v in exons_per.values() if v == 1) / len(exons_per)) if exons_per else nan

prot = read_fasta(inp.pep)
lens = [len(s.rstrip("*")) for s in prot.values()]
n_prot = len(prot)
internal_stop = sum(1 for s in prot.values() if "*" in s.rstrip("*"))
starts_m = sum(1 for s in prot.values() if s.startswith("M"))
ends_stop = sum(1 for s in prot.values() if s.endswith("*"))

pfam_all, ipr_all = parse_ipscan(inp.ipscan)
pfam_here = len(pfam_all & set(prot))
ipr_here = len(ipr_all & set(prot))


def frac(a, b):
    return a / b if b else nan


metrics = [
    ("sample", w.sample), ("haplotype", w.hap), ("track", w.track),
    ("n_sequences", n_seqs), ("n_genes", n_gene), ("n_mrna", n_mrna),
    ("n_exons", n_exon), ("n_cds", n_cds),
    ("mean_exons_per_mrna", f"{st.mean(exons_per.values()) if exons_per else nan:.3f}"),
    ("mono_exonic_frac", f"{mono_frac:.4f}"),
    ("n_proteins", n_prot),
    ("n_proteins_with_pfam", pfam_here),
    ("frac_proteins_with_pfam", f"{frac(pfam_here, n_prot):.4f}"),
    ("n_proteins_with_interpro", ipr_here),
    ("frac_proteins_with_interpro", f"{frac(ipr_here, n_prot):.4f}"),
    ("mean_prot_len_aa", f"{st.mean(lens) if lens else nan:.1f}"),
    ("median_prot_len_aa", f"{st.median(lens) if lens else nan:.1f}"),
    ("internal_stop_n", internal_stop),
    ("internal_stop_frac", f"{frac(internal_stop, n_prot):.4f}"),
    ("starts_with_M_frac", f"{frac(starts_m, n_prot):.4f}"),
    ("ends_with_stop_frac", f"{frac(ends_stop, n_prot):.4f}"),
]
with open(snakemake.output.qc, "w") as out:  # noqa: F821
    out.write("metric\tvalue\n")
    for k, v in metrics:
        out.write(f"{k}\t{v}\n")
