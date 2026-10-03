"""Keep the single best homology hit per protein.

Used by workflow 06, step 3. Takes tabular (BLAST -outfmt 6 layout) hits
from MMseqs2, DIAMOND or BLASTP and keeps, per query, the hit with the
highest bitscore (column 12). Several hits are requested upstream because
``--max-seqs 1`` / ``-max_target_seqs 1`` does not reliably return the top
hit.
"""

sample = snakemake.wildcards.sample  # noqa: F821
hap = snakemake.wildcards.hap  # noqa: F821
aligner = snakemake.params.aligner  # noqa: F821

best, n_raw = {}, 0
with open(snakemake.input.raw) as fh:  # noqa: F821
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        n_raw += 1
        try:
            bits = float(f[11])
        except ValueError:
            continue
        if f[0] not in best or bits > best[f[0]][0]:
            best[f[0]] = (bits, line)

with open(snakemake.output.best, "w") as out:  # noqa: F821
    for q in sorted(best):
        out.write(best[q][1])

msg = (f"[{sample}\t{hap}] {aligner}: {n_raw:,} hits -> "
       f"{len(best):,} queries with a best hit")
with open(snakemake.log[0], "w") as lg:  # noqa: F821
    lg.write(msg + "\n")
print("  " + msg)
