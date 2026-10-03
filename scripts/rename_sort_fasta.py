"""Normalise scaffold sequence IDs and sort them.

Used by workflow 05 (finishing). Renames RagTag scaffold headers so every
haplotype carries identical chromosome names:

    <ref>_Chr1_RagTag     ->  Chr01      (reference prefix stripped)
    <ref>_Chr1UA_RagTag   ->  Chr01_UA   (reference-unanchored bin)
    utg000123l            ->  utg000123l (unplaced unitig, kept verbatim)

Output order: Chr01..ChrNN, then Chr0N_UA bins, then unplaced unitigs,
wrapped at 60 columns. Writes an old->new ID map with lengths as an audit
trail, fails on duplicate IDs instead of silently overwriting, and warns if
the chromosome count differs from the expected number.

Note: because IDs repeat across haplotypes, do not concatenate haplotype
FASTAs without re-prefixing (see make_pansn_chromosomes.py).
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(snakemake.scriptdir, "lib"))  # noqa: F821
from genome_utils import write_fasta  # noqa: E402

sample = snakemake.wildcards.sample  # noqa: F821
hap = snakemake.wildcards.hap  # noqa: F821
n_chrom = int(snakemake.params.n_chrom)  # noqa: F821


def new_id(old_header):
    """Map one RagTag header to the normalised ID. Always returns a value."""
    m = re.search(r"[Cc]hr0*(\d+)\s*_?UA", old_header)       # UA bins first
    if m:
        return f"Chr{int(m.group(1)):02d}_UA"
    m = re.search(r"[Cc]hr0*(\d+)", old_header)              # placed chromosome
    if m:
        return f"Chr{int(m.group(1)):02d}"
    return re.sub(r"_RagTag$", "", old_header)               # unplaced unitig


def sort_key(k):
    m = re.fullmatch(r"Chr(\d+)(_UA)?", k)
    if m:
        return (1 if m.group(2) else 0, int(m.group(1)), k)
    return (2, 0, k)


seqs, id_pairs, n_records = {}, [], 0
current_id, current_seq = None, []
with open(snakemake.input.fa) as fin:  # noqa: F821
    for line in fin:
        if line.startswith(">"):
            if current_id is not None:
                seqs[current_id] = "".join(current_seq)
            old_header = line[1:].strip().split()[0]
            current_seq = []
            n_records += 1
            current_id = new_id(old_header)
            if current_id in seqs:                     # fail loud, never overwrite
                raise ValueError(f"[{sample} {hap}] duplicate ID '{current_id}' "
                                 f"generated from header '{old_header}'")
            id_pairs.append((old_header, current_id))
        else:
            current_seq.append(line.strip())
    if current_id is not None:
        seqs[current_id] = "".join(current_seq)

write_fasta(snakemake.output.sorted_fa, seqs,  # noqa: F821
            order=sorted(seqs, key=sort_key))

with open(snakemake.output.id_map, "w") as mout:  # noqa: F821
    mout.write("old_header\tnew_id\tlength_bp\n")
    for old, new in id_pairs:
        mout.write(f"{old}\t{new}\t{len(seqs[new])}\n")

chroms = sorted(k for k in seqs if re.fullmatch(r"Chr\d+", k))
ua = [k for k in seqs if k.endswith("_UA")]
other = [k for k in seqs if k not in chroms and k not in ua]
msgs = [
    f"[{sample}] {hap}: {n_records} records in, {len(seqs)} written",
    f"  chromosomes : {len(chroms)}  ({', '.join(chroms) if chroms else 'none'})",
    f"  UA bins     : {len(ua)}",
    f"  unplaced    : {len(other)}",
]
if n_records != len(seqs):
    msgs.append(f"  *** WARNING: {n_records - len(seqs)} records lost ***")
if len(chroms) != n_chrom:
    msgs.append(f"  *** WARNING: expected {n_chrom} chromosomes, found {len(chroms)} ***")
with open(snakemake.log[0], "w") as log_out:  # noqa: F821
    log_out.write("\n".join(msgs) + "\n")
print("\n".join(msgs))
