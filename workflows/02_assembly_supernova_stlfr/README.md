# 02 · Assembly with Supernova (stLFR linked reads)

De novo assembly of BGI stLFR linked reads with Supernova, followed by
haplotig purging to remove allelic duplication from heterozygous genomes,
and evaluation.

## Steps

| Rule | What it does | Tool |
|---|---|---|
| `stlfr2supernova` | Converts stLFR barcodes to 10x format and runs Supernova, driven by a per-sample profile file | [stlfr2supernova_pipeline](https://github.com/BGI-Qingdao/stlfr2supernova_pipeline), Supernova |
| `supernova_mkoutput` | Writes the assembly graph as one pseudohaplotype FASTA | Supernova |
| `map_reads` | Maps the cleaned short reads back to the draft | minimap2, samtools |
| `purge_hist` | Read-depth histogram for choosing the purging cut-offs | purge_haplotigs |
| `purge_cov` | Flags suspect contigs from the low / mid / high depth cut-offs | purge_haplotigs |
| `purge_purge` | Removes haplotigs and artefactual contigs → `curated.fasta` | purge_haplotigs |
| `min_length_filter` | Drops contigs shorter than `min_contig_len` | SeqKit |
| `purge_stats` | Contiguity before and after the length filter | BBMap `stats.sh` |
| `busco`, `busco_plot` | Completeness of both versions and a comparison plot | BUSCO |

## Input

stLFR reads described in one profile file per sample (read files and
primers), as required by the stlfr2supernova pipeline.

## Configuration — `config/02_assembly_supernova_stlfr.yaml`

| Key | Meaning |
|---|---|
| `stlfr2supernova.run_script` | Path to the pipeline's `run.sh` |
| `stlfr2supernova.profile` | Profile file per sample, e.g. `config/stlfr_profiles/{sample}.profile` |
| `supernova.asm_dir`, `clean_r1`, `clean_r2` | Paths the pipeline writes, relative to the sample folder |
| `supernova.minsize`, `supernova.style` | `mkoutput` options; `pseudohap` gives one FASTA for purging |
| `purge_haplotigs.default` / `.<sample>` | Low, mid and high depth cut-offs |
| `min_contig_len` | Minimum contig length kept after purging |
| `busco.lineage` | BUSCO lineage |

## Choosing the purging cut-offs

The cut-offs depend on each dataset's read-depth distribution:

```bash
snakemake -s workflows/02_assembly_supernova_stlfr/Snakefile \
          --configfile config/02_assembly_supernova_stlfr.yaml --use-conda --cores 80 \
          --until purge_hist
# inspect results/02_assembly_supernova_stlfr/<sample>/purge/aligned.bam.histogram.png
# set purge_haplotigs in the config, then run the full workflow
```

## Output — `results/02_assembly_supernova_stlfr/<sample>/`

| Path | Content |
|---|---|
| `<sample>_output.fasta` | Supernova pseudohaplotype assembly |
| `purge/curated.fasta` | Purged assembly |
| `purge/curated_min_len.fasta` | Purged assembly without short contigs |
| `purge/*.stats.txt` | Contiguity statistics of both |
| `busco/busco_figure.png` | BUSCO comparison plot |

## Requirements

Supernova and the stlfr2supernova pipeline are not on Conda; install them
separately and make `supernova` available on the `PATH`.
