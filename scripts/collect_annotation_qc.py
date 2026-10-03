"""Collect per-haplotype annotation QC into sample-level tables and flags.

Used by workflow 06, step 13. Writes, for one sample:

  * <sample>_annotation_qc.tsv      haplotype x track metrics + BUSCO line
  * <sample>_compleasm.tsv          all compleasm rows
  * <sample>_filtering_summary.tsv  gene counts per filtering stage
  * <sample>_annotation_qc_flags.txt  [ok]/[FLAG] checks:
      - gene-count spread across haplotypes on 'chrom' above a threshold
      - internal-stop fraction above a threshold (CDS phase/strand problems)
      - Pfam domain coverage below 40% on 'chrom'
      - more than 25% of genes on unplaced sequence (scaffolding problems)
"""

sample = snakemake.wildcards.sample  # noqa: F821
p = snakemake.params  # noqa: F821
haps, tracks = list(p.haps), list(p.tracks)
ann_dir = p.ann_dir
hap_spread_max = float(p.hap_spread_max)
stop_frac_max = float(p.stop_frac_max)
outp = snakemake.output  # noqa: F821


def path(hap, suffix):
    return f"{ann_dir}/{hap}/{sample}_{hap}.{suffix}"


per = {}
for hap in haps:
    for track in tracks:
        d = {}
        with open(path(hap, f"{track}.model_qc.tsv")) as fh:
            for line in fh:
                k, _, v = line.rstrip("\n").partition("\t")
                if k != "metric":
                    d[k] = v
        per[(hap, track)] = d

busco = {}
for hap in haps:
    for track in tracks:
        try:
            with open(path(hap, f"{track}.compleasm.tsv")) as fh:
                rows_ = [line.rstrip("\n").split("\t") for line in fh]
            busco[(hap, track)] = dict(zip(rows_[0], rows_[1]))
        except (OSError, IndexError):
            busco[(hap, track)] = {}

cols = ["n_sequences", "n_genes", "n_mrna", "n_exons", "n_cds",
        "mean_exons_per_mrna", "mono_exonic_frac", "n_proteins",
        "n_proteins_with_pfam", "frac_proteins_with_pfam",
        "n_proteins_with_interpro", "frac_proteins_with_interpro",
        "mean_prot_len_aa", "median_prot_len_aa",
        "internal_stop_n", "internal_stop_frac",
        "starts_with_M_frac", "ends_with_stop_frac"]

with open(outp.tsv, "w") as out:
    out.write("sample\thaplotype\ttrack\t" + "\t".join(cols) + "\tCOMPLEASM\n")
    for hap in haps:
        for track in tracks:
            out.write(f"{sample}\t{hap}\t{track}\t"
                      + "\t".join(per[(hap, track)].get(c, "NA") for c in cols)
                      + "\t" + busco[(hap, track)].get("COMPLEASM", "NA") + "\n")


def concat_tables(dest, suffixes, header):
    with open(dest, "w") as out:
        out.write(header)
        for hap in haps:
            for suffix in suffixes:
                try:
                    with open(path(hap, suffix)) as fh:
                        for i, line in enumerate(fh):
                            if i:
                                out.write(line)
                except OSError:
                    pass


concat_tables(outp.cmp, [f"{t}.compleasm.tsv" for t in tracks],
              "sample\thaplotype\ttrack\ttype\tlineage\tCOMPLEASM\tS_n\tD_n\tF_n\tM_n\tN\n")
concat_tables(outp.filt, ["filtering_summary.tsv"],
              "sample\thaplotype\tstage\tn_genes\tpct_of_raw\n")

# ---------------- flags ----------------
msgs = []
for track in tracks:
    counts = []
    for hap in haps:
        try:
            counts.append(int(per[(hap, track)]["n_genes"]))
        except (KeyError, ValueError):
            pass
    if len(counts) == len(haps) and max(counts) > 0:
        spread = (max(counts) - min(counts)) / max(counts)
        line = f"{track}: gene count across haplotypes {min(counts):,}-{max(counts):,} ({spread:.1%})"
        if track == "chrom" and spread > hap_spread_max:
            msgs.append(f"[FLAG] {line} > {hap_spread_max:.0%} - check assembly "
                        "balance / prediction consistency")
        else:
            msgs.append(f"[ok]   {line}")

for hap in haps:
    for track in tracks:
        try:
            isf = float(per[(hap, track)]["internal_stop_frac"])
        except (KeyError, ValueError):
            continue
        if isf > stop_frac_max:
            msgs.append(f"[FLAG] {hap} {track}: internal-stop fraction {isf:.1%} - "
                        "CDS phase or strand handling likely wrong")

for hap in haps:
    try:
        fp = float(per[(hap, "chrom")]["frac_proteins_with_pfam"])
        n = int(per[(hap, "chrom")]["n_proteins_with_pfam"])
        line = f"{hap} chrom: {fp:.1%} of proteins carry a Pfam domain ({n:,})"
        msgs.append(f"[FLAG] {line} - unusually low; check models/InterProScan"
                    if fp < 0.40 else f"[ok]   {line}")
    except (KeyError, ValueError):
        pass

for hap in haps:
    try:
        g_chrom = int(per[(hap, "chrom")]["n_genes"])
        g_utg = int(per[(hap, "utg")]["n_genes"])
    except (KeyError, ValueError):
        continue
    tot = g_chrom + g_utg
    if tot:
        f_utg = g_utg / tot
        line = f"{hap}: {f_utg:.1%} of genes on unplaced sequence ({g_utg:,}/{tot:,})"
        msgs.append(f"[FLAG] {line} - unusually high; check scaffolding"
                    if f_utg > 0.25 else f"[ok]   {line}")

with open(outp.flags, "w") as out:
    out.write(f"Annotation QC - {sample}\n" + "=" * 60 + "\n")
    out.write("\n".join(f"{sample}\t{m}" for m in msgs) + "\n"
              if msgs else f"{sample}\tno checks could be computed\n")
print(open(outp.flags).read())
