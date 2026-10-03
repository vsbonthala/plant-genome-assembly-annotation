# 08 · Pan-genome TE annotation with panEDTA

Consistent TE annotation across all genomes of a pan-genome with panEDTA,
part of [EDTA](https://github.com/oushujun/EDTA).

## Why panEDTA

Independent EDTA runs (workflow 07) produce a separate library for each
genome, so the same TE family can be named and classified differently from
genome to genome. panEDTA builds one pan-genome TE library from all genomes
and re-annotates every genome with it, so TE annotations can be compared
directly across genomes and haplotypes.

## Steps

| Rule | What it does | Tool |
|---|---|---|
| `genome_list` | Writes `genome_list.txt` with the absolute path of every genome | Python |
| `panedta` | Annotates each genome with EDTA, combines the libraries into a pan-genome TE library (using the species CDS to remove gene sequences), and re-annotates every genome with it | `panEDTA.sh` |

## Input

* All genome or haplotype FASTAs, e.g. the finished haplotypes from
  workflow 05.
* Coding sequences of the species, e.g. `<sample>_<hap>.cds.fa` from
  workflow 06 or a reference annotation.

## Configuration — `config/08_pan_repeats_panedta.yaml`

| Key | Meaning |
|---|---|
| `genomes` | Name → FASTA path for every genome in the pan-genome |
| `cds` | CDS FASTA (`panEDTA.sh -c`) |
| `panedta.cmd` | Path to `panEDTA.sh` if it is not on the `PATH` |
| `panedta.extra` | Further options, e.g. a curated library; see `panEDTA.sh -h` |
| `panedta.threads`, `panedta.mem_mb` | Resources |

## Output — `results/08_pan_repeats_panedta/`

| File | Content |
|---|---|
| `genome_list.txt` | Genomes passed to panEDTA |
| `panEDTA.log` | Full panEDTA log |
| `panEDTA.done` | Marker written when panEDTA finished successfully |
| other files | Pan-genome TE library and per-genome annotations, as written by panEDTA |

## Run

```bash
snakemake -s workflows/08_pan_repeats_panedta/Snakefile \
          --configfile config/08_pan_repeats_panedta.yaml --use-conda --cores 10
```

The command run is the same as calling panEDTA directly:

```bash
panEDTA.sh -g genome_list.txt -c cds.fasta -t 10
```
