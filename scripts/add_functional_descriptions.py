"""Attach homology-derived functional descriptions to genes and sequences.

Used by workflow 06, step 4. For each protein's best hit, looks up the
reference description in the UniProt/Swiss-Prot FASTA (trimming the
OS=/OX=/GN=/PE=/SV= tail) and:

  * adds ``Note=<description>`` to the mRNA and its parent gene in the GFF3
    (';' and '=' are neutralised, as they are GFF3 delimiters);
  * appends the description to protein and transcript FASTA headers.

Proteins without a hit become 'hypothetical protein'. Warns if most hit IDs
cannot be found in the reference FASTA, which means the search database was
built from a different FASTA.
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(snakemake.scriptdir, "lib"))  # noqa: F821
from genome_utils import gff_attr, id_variants, load_gff, mrna_to_gene, write_gff  # noqa: E402

sample = snakemake.wildcards.sample  # noqa: F821
hap = snakemake.wildcards.hap  # noqa: F821
uniprot_fa = snakemake.params.uniprot_fa  # noqa: F821
inp, outp = snakemake.input, snakemake.output  # noqa: F821

# 1. reference ID -> description, registered under every ID form
desc, n_entries = {}, 0
with open(uniprot_fa) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        parts = line[1:].rstrip("\n").split(None, 1)
        if not parts:
            continue
        d = parts[1] if len(parts) > 1 else ""
        d = re.split(r'\s+(?:OS|OX|GN|PE|SV)=', d)[0].strip() or "hypothetical protein"
        for v in id_variants(parts[0]):
            desc.setdefault(v, d)
        n_entries += 1
if n_entries == 0:
    raise ValueError(f"no FASTA headers parsed from {uniprot_fa}")

# 2. query mRNA -> description via its best hit
q2d, n_hit, n_miss = {}, 0, 0
with open(inp.best) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        d = desc.get(f[1])
        if d is None:
            n_miss += 1
            d = "hypothetical protein"
        else:
            n_hit += 1
        q2d[f[0]] = d

# 3. GFF3: Note= on mRNA and gene rows
header, rows = load_gff(inp.gff)
m2g = mrna_to_gene(rows)
g2d = {m2g[m]: d for m, d in q2d.items() if m in m2g}
out_rows = []
for p in rows:
    fid = gff_attr(p[8], "ID")
    d = q2d.get(fid) if p[2] in ("mRNA", "transcript") else g2d.get(fid) if p[2] == "gene" else None
    if d and "Note=" not in p[8]:
        safe = d.replace(";", ",").replace("=", "-")
        p = p[:8] + [p[8].rstrip(";") + f";Note={safe}"]
    out_rows.append(p)
write_gff(outp.gff, header, out_rows)

# 4. FASTA headers
for src, dst in ((inp.pep, outp.pep), (inp.trn, outp.trn)):
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                fo.write(f">{name} {q2d.get(name, 'hypothetical protein')}\n")
            else:
                fo.write(line)

msgs = [f"[{sample}\t{hap}] reference entries parsed: {n_entries:,}",
        f"  queries with a hit          : {len(q2d):,}",
        f"  target id resolved          : {n_hit:,}",
        f"  target id NOT in reference  : {n_miss:,}"]
if n_miss > 0.5 * (n_hit + n_miss):
    msgs.append("  *** WARNING: majority of target ids did not resolve - "
                "is the search database built from THIS FASTA? ***")
with open(snakemake.log[0], "w") as lg:  # noqa: F821
    lg.write("\n".join(msgs) + "\n")
print("\n".join(msgs))
