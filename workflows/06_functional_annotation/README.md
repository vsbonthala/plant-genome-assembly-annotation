# 06 · Functional annotation, TE filtering and QC

Turns raw Helixer gene models into functionally annotated, TE-filtered gene
sets per haplotype, delivered in three tracks, with quality checks.

## Steps

| # | Rule | What it does | Tool / script |
|---|---|---|---|
| 0 | `feature_inventory` | Feature types, CDS phases, strands and genes per track before processing | `gff_feature_inventory.py` |
| 1 | `normalise_gff`, `sort_gff` | Tidies the gene hierarchy; orders genes Chr01…ChrNN → UA → unitigs | AGAT or GenomeTools; `sort_gff.py` |
| 2 | `extract_sequences` | Proteins, CDS and transcripts | gffread |
| 3 | `homology_search`, `best_hit` | Search against UniProt/Swiss-Prot; best hit per protein by bitscore | MMseqs2, DIAMOND or BLASTP; `best_hit.py` |
| 4 | `add_functional_descriptions` | `Note=` descriptions in the GFF3 and FASTA headers | `add_functional_descriptions.py` |
| 5 | `interproscan` | Domains and GO terms | InterProScan |
| 6 | `merge_interproscan` | `Dbxref` (Pfam, InterPro) and `Ontology_term` (GO) in the GFF3 | `merge_interproscan_gff.py` |
| 7 | `filter_short_models` | Drops genes with no protein ≥ `min_aa` | `filter_short_models.py` |
| 8 | `filter_te_pfam` | Drops genes with TE-associated Pfam domains | `filter_te_pfam.py` |
| 9 | `filter_te_keywords` | Drops genes whose annotation names TE functions | `filter_te_keywords.py` |
| 10 | `final_nonte_sequences` | Final non-TE proteins and transcripts | `final_nonte_sequences.py` |
| 11 | `split_tracks` | `chrom` (placed), `utg` (UA + unitigs), `all` | `split_tracks.py` |
| 12 | `filtering_summary`, `model_qc`, `compleasm_proteins`, `parse_compleasm_proteins` | Genes per filtering stage; model metrics; completeness | scripts + compleasm |
| 13 | `collect_annotation_qc` | Sample-level tables and `[ok]` / `[FLAG]` checks | `collect_annotation_qc.py` |

All scripts are in `scripts/`; shared helpers in `scripts/lib/genome_utils.py`.

### Design choices

- **Annotate once, split last.** Homology search and InterProScan run once
  on the full set; the track split happens after filtering, halving compute.
- **Exact ID matching.** Filters use exact ID sets and the GFF3 `Parent`
  attribute, so removing `gene1` never removes `gene10`.
- **Explicit best hit.** Several hits are kept per protein and the best
  bitscore wins, since `max-target-seqs 1` does not guarantee the top hit.

### QC flags (`<sample>_annotation_qc_flags.txt`)

| Check | Flagged when |
|---|---|
| Gene-count spread across haplotypes (`chrom`) | above `hap_spread_max` (default 15%) |
| Internal stop codons | above `stop_frac_max` (default 5%) |
| Pfam domain coverage (`chrom`) | below 40% of proteins |
| Genes on unplaced sequence | above 25% |

## Input

From workflow 05: `<sample>_<hap>.fa` and `<sample>_<hap>.Helixer.gff3`
(paths set by `genome_fasta` and `gene_models`).

## Configuration — `config/06_functional_annotation.yaml`

| Key | Meaning |
|---|---|
| `normaliser` | `agat`, `gt` or `none` |
| `homology.aligner` | `mmseqs`, `diamond` or `blastp` |
| `homology.database`, `homology.uniprot_fasta` | Search database and the FASTA it was built from |
| `interproscan.cmd`, `applications`, `write_analyses` | InterProScan call and which analyses go into the GFF3 (must include Pfam) |
| `min_aa` | Minimum protein length |
| `te_pfam_list`, `te_keywords` | TE filters |
| `compleasm.*`, `hap_spread_max`, `stop_frac_max` | QC settings |

Build the search database once from the same FASTA as `uniprot_fasta`
(commands in the config file).

## Output — `results/06_functional_annotation/<sample>/`

| Path | Content |
|---|---|
| `<hap>/<sample>_<hap>.{all,chrom,utg}.nonTE.gff3` | Final gene models per track |
| `<hap>/<sample>_<hap>.{track}.nonTE.proteins.fa` / `.transcripts.fa` | Final sequences |
| `<hap>/<sample>_<hap>.cds.fa` | CDS of all models (usable as panEDTA input, workflow 08) |
| `<sample>_annotation_qc.tsv` | Metrics per haplotype and track |
| `<sample>_filtering_summary.tsv` | Genes remaining after each filter |
| `<sample>_compleasm.tsv` | Completeness per haplotype and track |
| `<sample>_annotation_qc_flags.txt` | Checks |

## Requirements

InterProScan is installed separately. Replace
`resources/te_pfam_ids.example.txt` with your own curated list.
