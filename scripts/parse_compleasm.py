"""Fold a compleasm summary.txt into one standard BUSCO line.

Used by the assembly workflows (01, 03, 04) and workflow 06. compleasm reports S/D/F/I/M/N on separate
lines; this writes a one-row table with the familiar string

    C:99.76%[S:3.04%,D:96.72%],F:0.12%,M:0.12%,n:1614

where C = S + D. The 'I' (incomplete) field is deliberately excluded.
"""

import os
import sys

sys.path.insert(0, os.path.join(snakemake.scriptdir, "lib"))  # noqa: F821
from genome_utils import busco_one_line, parse_compleasm_summary  # noqa: E402

w = snakemake.wildcards  # noqa: F821
hap = getattr(w, "hap", "NA")
track = getattr(w, "track", "NA")
kind = snakemake.params.kind  # noqa: F821

lineage, pct, cnt, n = parse_compleasm_summary(snakemake.input.summary)  # noqa: F821
line = busco_one_line(pct, n)
with open(snakemake.output.tsv, "w") as out:  # noqa: F821
    out.write("sample\thaplotype\ttrack\ttype\tlineage\tCOMPLEASM\tS_n\tD_n\tF_n\tM_n\tN\n")
    out.write(f"{w.sample}\t{hap}\t{track}\t{kind}\t{lineage or 'NA'}\t{line}\t"
              f"{cnt.get('S', 'NA')}\t{cnt.get('D', 'NA')}\t{cnt.get('F', 'NA')}\t"
              f"{cnt.get('M', 'NA')}\t{n if n is not None else 'NA'}\n")
print(f"  [{w.sample}\t{hap}\t{track}\t{kind}] {line}")
