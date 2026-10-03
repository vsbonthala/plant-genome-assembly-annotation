# 05 · Finishing and gene prediction

Per-haplotype finishing of reference-guided scaffolds (consistent sequence
names, ordering, indexing, PanSN sets for pangenome graphs) and ab initio
gene prediction with Helixer.

## Steps

| Rule | What it does | Tool / script |
|---|---|---|
| `rename_sort_fasta` | Renames sequences so every haplotype uses identical names, sorts them, writes an old→new ID map, and stops on duplicate IDs | `scripts/rename_sort_fasta.py` |
| `index_and_stats` | FASTA index and contiguity statistics | samtools, BBMap |
| `make_pansn_chromosomes` | Placed chromosomes only, renamed `sample#hap#ChrNN` ([PanSN-spec](https://github.com/pangenome/PanSN-spec)) | `scripts/make_pansn_chromosomes.py` |
| `compress_pansn` | bgzip + faidx, as pggb requires | htslib, samtools |
| `helixer` | Ab initio gene prediction, `land_plant` model (GPU) | [Helixer](https://github.com/weberlab-hhu/Helixer) |

### Sequence naming

```
<ref>_Chr1_RagTag     ->  Chr01        placed chromosome (reference prefix removed)
<ref>_Chr1UA_RagTag   ->  Chr01_UA     sequence assigned to Chr01 but not anchored
utg000123l            ->  utg000123l   unplaced unitig, unchanged
```

Order: Chr01…ChrNN, then the UA bins, then unplaced unitigs. Identical
names across haplotypes allow direct haplotype-vs-haplotype and
haplotype-vs-reference comparison. Because names repeat across haplotypes,
use the PanSN files (not the per-haplotype FASTAs) when combining genomes.

### Building a pangenome graph from the PanSN files

```bash
# one chromosome from every sample and haplotype
for f in results/05_finishing_helixer/*/*/*.chroms.pansn.fa.gz; do
  samtools faidx "$f" $(cut -f1 "$f.fai" | grep '#Chr01$')
done > chr01.pan.fa
bgzip chr01.pan.fa && samtools faidx chr01.pan.fa.gz
pggb -i chr01.pan.fa.gz -n <number of haplotypes> -o pggb_chr01
```

## Input

Reference-guided scaffolds per sample and haplotype (for example RagTag
output restricted to nuclear sequence), found through the `scaffolds`
pattern, e.g. `data/scaffolds/{sample}/{hap}/ragtag.scaffold.nuclear.fa`.

## Configuration — `config/05_finishing_helixer.yaml`

| Key | Meaning |
|---|---|
| `haplotypes`, `pansn_hap_ids` | Haplotype labels and their integer PanSN IDs |
| `n_chrom` | Expected chromosome number (sanity check, logged) |
| `scaffolds` | Path pattern to the scaffolds |
| `helixer.setup`, `helixer.cmd` | How to activate and call Helixer |
| `helixer.lineage`, `helixer.threads` | Helixer model and threads |

## Output — `results/05_finishing_helixer/<sample>/<hap>/`

| File | Content |
|---|---|
| `<sample>_<hap>.fa` (+ `.fai`) | Finished haplotype |
| `<sample>_<hap>.id_map.tsv` | Old → new IDs with lengths |
| `<sample>_<hap>.stats.txt` | Contiguity statistics |
| `<sample>_<hap>.chroms.pansn.fa.gz` (+ `.fai`) | PanSN chromosome set for pggb |
| `<sample>_<hap>.Helixer.gff3` | Gene models |

## Requirements

Helixer needs a GPU and is installed from its container or a Python
virtualenv; set `helixer.setup` to activate it.
