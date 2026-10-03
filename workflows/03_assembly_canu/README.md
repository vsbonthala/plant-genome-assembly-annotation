# 03 · Assembly with Canu

De novo assembly of noisy or accurate long reads (PacBio CLR, PacBio HiFi or
ONT) with Canu's correction, trimming and assembly stages.

## Steps

| Rule | What it does | Tool |
|---|---|---|
| `filter_reads` | Removes reads shorter than `min_read_length` | SeqKit |
| `canu` | Read correction, trimming and assembly → `<sample>.contigs.fasta` | Canu |
| `compleasm_genome`, `parse_compleasm_genome` | Gene-space completeness, one BUSCO line | compleasm |
| `assembly_stats` | Contiguity statistics | BBMap `stats.sh` |

## Configuration — `config/03_assembly_canu.yaml`

| Key | Meaning |
|---|---|
| `reads` | Path pattern to each sample's reads |
| `min_read_length` | Minimum read length kept (0 = keep all) |
| `canu.genome_size` | `genomeSize`, e.g. `850m` |
| `canu.read_type` | `-pacbio`, `-pacbio-hifi` or `-nanopore` |
| `canu.max_memory`, `meryl_memory`, `meryl_threads` | Canu resource limits (GB / threads) |
| `canu.use_grid`, `canu.extra` | Let Canu submit its own cluster jobs (see below) |
| `threads.canu` | `maxThreads` |

## Running Canu on a cluster

By default the whole Canu run is one Snakemake job (`useGrid=false`), which
Snakemake can submit to Slurm like any other step. Alternatively, Canu can
split itself into many cluster jobs: set `use_grid: true` and pass your
scheduler options, e.g.

```yaml
canu:
  use_grid: true
  extra: 'gridOptions="--account=<acct> --partition=<part> --time=21-00:00:00"'
```

## Output — `results/03_assembly_canu/<sample>/`

| Path | Content |
|---|---|
| `canu/<sample>.contigs.fasta` | Canu contigs |
| `canu/` | All Canu intermediate files and reports |
| `eval/<sample>.canu.compleasm.tsv`, `eval/<sample>.canu.stats.txt` | Evaluation |

## Run

```bash
snakemake -s workflows/03_assembly_canu/Snakefile \
          --configfile config/03_assembly_canu.yaml --use-conda --cores 20
```
