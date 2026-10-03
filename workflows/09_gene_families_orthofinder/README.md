# 09 · Gene families with OrthoFinder

Infers gene families (orthogroups), gene and species trees, hierarchical
orthogroups, orthologues and gene duplication events across genomes,
haplotypes or species with [OrthoFinder](https://github.com/davidemms/OrthoFinder).

## Steps

| Rule | What it does | Tool |
|---|---|---|
| `stage_proteome` | Copies each protein FASTA into one clean input folder as `<species>.fa`; OrthoFinder uses file names as species names | cp |
| `orthofinder` | All-vs-all protein search, orthogroup clustering, gene trees, rooted species tree, hierarchical orthogroups (HOGs), orthologues, duplication events and summary statistics | OrthoFinder |

The command run is

```bash
orthofinder -t <threads> -n <run_name> -y -z -f <proteomes folder>
```

| Option | Meaning |
|---|---|
| `-y` | Split paralogous clades below the root of a HOG into separate HOGs |
| `-z` | Do not trim the multiple sequence alignments |

## Input

One protein FASTA per genome, haplotype or species. The final protein sets of
workflow 06 (`<sample>_<hap>.all.nonTE.proteins.fa`) can be listed directly,
together with proteomes of reference or outgroup species.

For gene-family counts, consider using one protein per gene (e.g. the longest
isoform); otherwise alternative isoforms of the same gene are counted as
separate family members.

## Configuration — `config/09_gene_families_orthofinder.yaml`

| Key | Meaning |
|---|---|
| `proteomes` | Species name → protein FASTA |
| `orthofinder.run_name` | `-n`: suffix of the results folder |
| `orthofinder.extra` | Further options; default `-y -z` |
| `orthofinder.threads`, `orthofinder.mem_mb` | Resources |

## Output — `results/09_gene_families_orthofinder/proteomes/OrthoFinder/Results_<run_name>/`

| Path | Content |
|---|---|
| `Orthogroups/Orthogroups.tsv` | Genes in each orthogroup, per species |
| `Orthogroups/Orthogroups.GeneCount.tsv` | Gene counts per orthogroup and species (gene-family sizes) |
| `Orthogroups/Orthogroups_UnassignedGenes.tsv` | Species-specific genes not in any orthogroup |
| `Phylogenetic_Hierarchical_Orthogroups/N0.tsv` | Hierarchical orthogroups at the root of the species tree |
| `Species_Tree/SpeciesTree_rooted.txt` | Rooted species tree |
| `Orthologues/` | Pairwise orthologue tables |
| `Gene_Duplication_Events/` | Duplication events mapped onto the species tree |
| `Comparative_Genomics_Statistics/` | Summary statistics, e.g. `Statistics_Overall.tsv` |

`Orthogroups.GeneCount.tsv` is the usual starting point for gene-family
expansion and contraction analyses (for example with CAFE).

## Run

```bash
snakemake -s workflows/09_gene_families_orthofinder/Snakefile \
          --configfile config/09_gene_families_orthofinder.yaml --use-conda --cores 40
```
