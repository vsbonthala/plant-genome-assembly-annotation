"""Shared helpers for the finishing and functional-annotation workflows.

Pure Python standard library: FASTA/GFF3 parsing and writing, sequence-ID
ordering, gene-block handling, and parsers for InterProScan and compleasm
output. Imported by the scripts in ../ via ``sys.path``.
"""

import os
import re
from collections import defaultdict

# -----------------------------------------------------------------------------
# FASTA
# -----------------------------------------------------------------------------

def read_fasta(path):
    """Return an insertion-ordered dict {first header token: sequence}."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(path, seqs, keep=None, wrap=60, order=None):
    """Write sequences wrapped at ``wrap`` columns; return the number written.

    ``keep``  - optional set of names to keep.
    ``order`` - optional iterable of names giving the output order.
    """
    names = order if order is not None else seqs.keys()
    n = 0
    with open(path, "w") as out:
        for name in names:
            if keep is not None and name not in keep:
                continue
            s = seqs[name]
            out.write(f">{name}\n")
            for i in range(0, len(s), wrap):
                out.write(s[i:i + wrap] + "\n")
            n += 1
    return n


# -----------------------------------------------------------------------------
# Sequence-ID conventions: Chr01..ChrNN, Chr0N_UA, unplaced unitigs
# -----------------------------------------------------------------------------

def seq_track(seqid):
    """'chrom' for placed chromosomes (ChrNN), otherwise 'utg'."""
    return "chrom" if re.fullmatch(r"Chr\d+", seqid) else "utg"


def seqid_sort_key(seqid):
    """Chr01..ChrNN first, then Chr0N_UA bins, then unplaced unitigs."""
    m = re.fullmatch(r"Chr(\d+)", seqid)
    if m:
        return (0, int(m.group(1)), "")
    m = re.fullmatch(r"Chr(\d+)_UA", seqid)
    if m:
        return (1, int(m.group(1)), "")
    return (2, 0, seqid)


# -----------------------------------------------------------------------------
# GFF3
# -----------------------------------------------------------------------------

def gff_attr(attrs, key):
    m = re.search(rf'(?:^|;)\s*{key}=([^;]+)', attrs)
    return m.group(1).strip() if m else None


def load_gff(path):
    header, rows = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
            elif line.strip():
                p = line.rstrip("\n").split("\t")
                if len(p) == 9:
                    rows.append(p)
    return header, rows


def write_gff(path, header, rows):
    with open(path, "w") as out:
        for h in header:
            out.write(h + "\n")
        for p in rows:
            out.write("\t".join(p) + "\n")


def gene_blocks(rows):
    """Group rows into per-gene blocks, preserving within-gene feature order."""
    blocks, current = [], None
    for p in rows:
        if p[2] == "gene":
            if current:
                blocks.append(current)
            current = [p]
        elif current is not None:
            current.append(p)
        else:
            blocks.append([p])
    if current:
        blocks.append(current)
    return blocks


def gene_id_of_block(block):
    return gff_attr(block[0][8], "ID")


def mrna_to_gene(rows):
    """mRNA/transcript ID -> parent gene ID, taken from the GFF3 itself."""
    m2g = {}
    for p in rows:
        if p[2] in ("mRNA", "transcript"):
            mid = gff_attr(p[8], "ID")
            par = gff_attr(p[8], "Parent")
            if mid and par:
                m2g[mid] = par.split(",")[0]
    return m2g


def drop_genes(rows, bad_gene_ids):
    """Remove whole gene blocks whose gene ID is in ``bad_gene_ids``."""
    return [p for b in gene_blocks(rows)
            if gene_id_of_block(b) not in bad_gene_ids
            for p in b]


def count_genes(rows):
    return sum(1 for p in rows if p[2] == "gene")


def add_attr(attrs, key, values):
    """Append comma-separated values to a GFF3 attribute, merging if present."""
    if not values:
        return attrs
    m = re.search(rf'(?:^|;)\s*{key}=([^;]*)', attrs)
    if m:
        existing = {v for v in m.group(1).split(",") if v}
        merged = ",".join(sorted(existing | set(values)))
        return attrs[:m.start(1)] + merged + attrs[m.end(1):]
    return attrs.rstrip(";") + f";{key}=" + ",".join(sorted(values))


# -----------------------------------------------------------------------------
# Homology-search reference IDs
# -----------------------------------------------------------------------------

def id_variants(tok):
    """Every ID form a homology-search tool might report for one reference
    FASTA header token (sp__ACC__NAME, sp|ACC|NAME, ACC, NAME), so the
    description lookup keeps working whichever aligner or FASTA is used."""
    out = {tok}
    norm = tok.replace("__", "|")
    out.add(norm)
    m = re.match(r'^(sp|tr)\|([^|]+)\|(.+)$', norm)
    if m:
        db, acc, name = m.groups()
        out.update({acc, name, f"{db}|{acc}|{name}", f"{db}__{acc}__{name}"})
    return out


# -----------------------------------------------------------------------------
# InterProScan TSV
# -----------------------------------------------------------------------------

def parse_ipscan(path):
    """InterProScan TSV -> (proteins with a Pfam hit, proteins with any IPR)."""
    pfam, ipr = set(), set()
    if not os.path.exists(path):
        return pfam, ipr
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            acc, analysis = f[0], f[3]
            if analysis.strip().lower().startswith("pfam"):
                pfam.add(acc)
            if len(f) > 11 and f[11].strip().startswith("IPR"):
                ipr.add(acc)
    return pfam, ipr


def parse_ipscan_annotations(path, include_analyses=None):
    """InterProScan TSV -> (mRNA -> Dbxref set, mRNA -> GO set).

    Columns used: 1 protein acc, 4 analysis, 5 signature acc,
    12 InterPro acc, 14 GO terms (pipe-separated). Pathways (col 15) are
    deliberately not written to the GFF3 to keep it compact.
    """
    xref, go = defaultdict(set), defaultdict(set)
    if not os.path.exists(path):
        return xref, go
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            acc, analysis, sig = f[0], f[3].strip(), f[4].strip()
            if include_analyses and analysis not in include_analyses:
                continue
            if sig and sig != "-":
                xref[acc].add(f"{analysis}:{sig}")
            if len(f) > 11 and f[11].strip().startswith("IPR"):
                xref[acc].add(f"InterPro:{f[11].strip()}")
            if len(f) > 13 and f[13].strip() not in ("", "-"):
                for g in f[13].split("|"):
                    g = g.strip()
                    if g.startswith("GO:"):
                        go[acc].add(g)
    return xref, go


# -----------------------------------------------------------------------------
# compleasm (BUSCO) summary
# -----------------------------------------------------------------------------

def parse_compleasm_summary(path):
    """compleasm summary.txt -> (lineage, pct dict, count dict, N)."""
    lineage, pct, cnt, n = None, {}, {}, None
    if not os.path.exists(path):
        return lineage, pct, cnt, n
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("## lineage:"):
                lineage = line.split(":", 1)[1].strip()
                continue
            m = re.match(r'^([SDFIM]):\s*([\d.]+)%\s*,\s*(\d+)\s*$', line)
            if m:
                pct[m.group(1)] = float(m.group(2))
                cnt[m.group(1)] = int(m.group(3))
                continue
            m = re.match(r'^N:\s*(\d+)\s*$', line)
            if m:
                n = int(m.group(1))
    return lineage, pct, cnt, n


def busco_one_line(pct, n):
    """Standard single-line BUSCO string; C = S + D. 'I' is excluded."""
    if not pct or n is None:
        return "NA"
    c = pct.get("S", 0.0) + pct.get("D", 0.0)
    return (f"C:{c:.2f}%[S:{pct.get('S', 0.0):.2f}%,D:{pct.get('D', 0.0):.2f}%],"
            f"F:{pct.get('F', 0.0):.2f}%,M:{pct.get('M', 0.0):.2f}%,n:{n}")
