# 04 · Assembly with Flye

De novo repeat-graph assembly of PacBio or ONT reads with Flye, including
its built-in polishing iterations.

## Steps

| Rule | What it does | Tool |
|---|---|---|
| `filter_reads` | Removes reads shorter than `min_read_length` | SeqKit |
| `flye` | Repeat-graph assembly and polishing → `assembly.fasta` | Flye |
| `compleasm_genome`, `parse_compleasm_genome` | Gene-space completeness, one BUSCO line | compleasm |
| `assembly_stats` | Contiguity statistics | BBMap `stats.sh` |

## Configuration — `config/04_assembly_flye.yaml`

| Key | Meaning |
|---|---|
| `reads` | Path pattern to each sample's reads |
| `min_read_length` | Minimum read length kept (0 = keep all) |
| `flye.read_type` | `pacbio-raw`, `pacbio-corr`, `pacbio-hifi`, `nano-raw`, `nano-corr` or `nano-hq` |
| `flye.genome_size` | `--genome-size`, e.g. `850m` |
| `flye.iterations` | Number of polishing iterations |
| `flye.extra` | Other Flye options |
| `threads.flye`, `mem_mb.flye` | Resources |

## Output — `results/04_assembly_flye/<sample>/`

| Path | Content |
|---|---|
| `flye/assembly.fasta` | Flye assembly |
| `flye/assembly_info.txt` | Per-contig length, coverage, circularity and repeat status |
| `eval/<sample>.flye.compleasm.tsv`, `eval/<sample>.flye.stats.txt` | Evaluation |

## Run

```bash
snakemake -s workflows/04_assembly_flye/Snakefile \
          --configfile config/04_assembly_flye.yaml --use-conda --cores 50
```

Because workflows 01, 03 and 04 share the same evaluation module, their
`eval/` outputs can be compared directly when choosing an assembler.
