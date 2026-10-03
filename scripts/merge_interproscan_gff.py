"""Merge InterProScan results into the GFF3.

Used by workflow 06, step 6. Adds ``Dbxref`` (e.g. Pfam:PF03151,
InterPro:IPR004853) and ``Ontology_term`` (GO IDs) to each mRNA, and rolls
them up to the parent gene. ``ipr_analyses`` limits which analyses are
written; Pfam must be included because the TE filter reads Pfam IDs from
the GFF3.
"""

import os
import sys

sys.path.insert(0, os.path.join(snakemake.scriptdir, "lib"))  # noqa: F821
from genome_utils import (add_attr, gff_attr, load_gff, mrna_to_gene,  # noqa: E402
                          parse_ipscan_annotations, write_gff)

sample = snakemake.wildcards.sample  # noqa: F821
hap = snakemake.wildcards.hap  # noqa: F821
analyses = snakemake.params.ipr_analyses  # noqa: F821

incl = {a.strip() for a in analyses.split(",") if a.strip()} or None
xref, go = parse_ipscan_annotations(snakemake.input.tsv, include_analyses=incl)  # noqa: F821

header, rows = load_gff(snakemake.input.gff)  # noqa: F821
m2g = mrna_to_gene(rows)

gx, gg = {}, {}
for m, vals in xref.items():
    if m in m2g:
        gx.setdefault(m2g[m], set()).update(vals)
for m, vals in go.items():
    if m in m2g:
        gg.setdefault(m2g[m], set()).update(vals)

out_rows, n_mrna, n_gene = [], 0, 0
for p in rows:
    fid = gff_attr(p[8], "ID")
    if p[2] in ("mRNA", "transcript") and fid in xref:
        a = add_attr(add_attr(p[8], "Dbxref", xref[fid]), "Ontology_term", go.get(fid, set()))
        p = p[:8] + [a]
        n_mrna += 1
    elif p[2] == "gene" and fid in gx:
        a = add_attr(add_attr(p[8], "Dbxref", gx[fid]), "Ontology_term", gg.get(fid, set()))
        p = p[:8] + [a]
        n_gene += 1
    out_rows.append(p)
write_gff(snakemake.output.gff, header, out_rows)  # noqa: F821

n_pfam = sum(1 for v in xref.values() if any(x.startswith("Pfam:") for x in v))
msgs = [f"[{sample}\t{hap}] InterProScan merged into GFF3",
        f"  analyses written            : {sorted(incl) if incl else 'all'}",
        f"  proteins with any signature : {len(xref):,}",
        f"  proteins with a Pfam domain : {n_pfam:,}",
        f"  mRNA rows annotated         : {n_mrna:,}",
        f"  gene rows annotated         : {n_gene:,}"]
if n_pfam == 0:
    msgs.append("  *** WARNING: no Pfam IDs written - the Pfam TE filter will "
                "have nothing to match ***")
with open(snakemake.log[0], "w") as lg:  # noqa: F821
    lg.write("\n".join(msgs) + "\n")
print("\n".join(msgs))
