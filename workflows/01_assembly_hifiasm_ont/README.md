# 01 · Assembly with Hifiasm (ONT reads)

Polyploid-aware de novo assembly of Oxford Nanopore reads with Hifiasm,
including read clean-up and assembly evaluation.

## Steps

| Rule | What it does | Tool |
|---|---|---|
| `porechop` | Trims ONT adapters | Porechop |
| `remove_split_reads` | Keeps reads with a unique ID; sets aside reads whose ID occurs more than once (split at middle adapters) | awk, seqtk |
| `nanostat` | Read length and quality summary of the filtered reads | NanoStat |
| `hifiasm` | Polyploid-aware assembly in ONT mode (`--n-hap`, `--hg-size`) | Hifiasm |
| `gfa_to_fasta` | Writes the primary unitig graph segments as FASTA | awk |
| `compleasm_genome`, `parse_compleasm_genome` | Gene-space completeness, summarised as one BUSCO line | compleasm |
| `assembly_stats` | Contiguity statistics (N50, L50, sizes) | BBMap `stats.sh` |

The two evaluation steps come from `workflows/common/assembly_eval.smk`,
shared with workflows 03 and 04 so all assemblers are evaluated identically.

## Input

Raw ONT reads, one FASTQ per sample, found through the `reads` pattern,
e.g. `data/ont/{sample}.fastq.gz`.

## Configuration — `config/01_assembly_hifiasm_ont.yaml`

| Key | Meaning |
|---|---|
| `samples` | Sample names |
| `reads` | Path pattern to each sample's reads |
| `hifiasm.genome_size` | Estimated haploid genome size (`--hg-size`), e.g. `840m` |
| `hifiasm.n_hap` | Ploidy (`--n-hap`); 4 for an autotetraploid |
| `hifiasm.extra` | Other Hifiasm options; default `--ont -s 0.4` |
| `compleasm.lineage`, `compleasm.library` | BUSCO lineage and download folder |
| `threads.*`, `mem_mb.*` | Threads and memory per step |

## Output — `results/01_assembly_hifiasm_ont/<sample>/`

| Path | Content |
|---|---|
| `qc/<sample>.porechop.uniq.fastq` | Cleaned reads used for assembly |
| `qc/<sample>.nanostat.txt` | Read statistics |
| `assembly/<sample>.bp.p_utg.gfa` / `.fa` | Primary unitig assembly (graph and FASTA) |
| `eval/<sample>.hifiasm.compleasm.tsv` | One-line completeness summary |
| `eval/<sample>.hifiasm.stats.txt` | Contiguity statistics |
| `logs/`, `benchmarks/` | Log and run time / memory of each step |

## Run

```bash
snakemake -s workflows/01_assembly_hifiasm_ont/Snakefile \
          --configfile config/01_assembly_hifiasm_ont.yaml --use-conda --cores 40
```
