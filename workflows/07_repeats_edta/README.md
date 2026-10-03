# 07 · Repeat annotation with EDTA

Builds a de novo transposable-element library for each genome and annotates
TEs genome-wide with [EDTA](https://github.com/oushujun/EDTA).

## Steps

| Rule | What it does | Tool |
|---|---|---|
| `edta` | Identifies LTR, TIR, Helitron and non-LTR elements, builds a filtered TE library (optionally merged with a curated library), and annotates the whole genome | EDTA |

Each genome is an independent job, so on a cluster many genomes run in
parallel; this replaces a Slurm job array.

## Configuration — `config/07_repeats_edta.yaml`

| Key | Meaning |
|---|---|
| `genomes` | Name → FASTA path for each genome or haplotype |
| `edta.species` | EDTA `--species` (`others`, `Rice` or `Maize`) |
| `edta.curated_lib` | Optional curated TE library, e.g. a plant subset of RepBase (`--curatedlib`) |
| `edta.sensitive`, `edta.anno` | `--sensitive` (RepeatModeler for remaining TEs) and `--anno` (whole-genome annotation) |
| `edta.threads`, `edta.mem_mb` | Resources |

## Output — `results/07_repeats_edta/<genome>/`

| File | Content |
|---|---|
| `<genome>.fa.mod.EDTA.TElib.fa` | Non-redundant TE library |
| `<genome>.fa.mod.EDTA.TEanno.gff3` | Whole-genome TE annotation |
| other `*.EDTA.*` files | Summaries and intermediate results from EDTA |

## Run

```bash
snakemake -s workflows/07_repeats_edta/Snakefile \
          --configfile config/07_repeats_edta.yaml --use-conda --cores 48
```

For annotations that are comparable across many genomes, use workflow 08
(panEDTA).
