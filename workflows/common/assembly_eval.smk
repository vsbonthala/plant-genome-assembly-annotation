# =============================================================================
# Shared assembly evaluation - included by workflows 01, 03 and 04
# =============================================================================
# Before `include:`, the including Snakefile defines:
#   ASSEMBLY  path pattern of the final assembly FASTA, containing {sample}
#   EVAL      path prefix pattern for evaluation outputs, containing {sample}
# and its config provides:
#   compleasm.lineage, compleasm.library, threads.compleasm
#
# Outputs per sample:
#   {EVAL}.compleasm/summary.txt   compleasm report
#   {EVAL}.compleasm.tsv           one-line BUSCO summary
#   {EVAL}.stats.txt               contiguity statistics (N50, L50, sizes)
# =============================================================================


rule compleasm_genome:
    """BUSCO-style gene-space completeness of the assembly with compleasm."""
    input:
        fa=ASSEMBLY,
    output:
        summary=f"{EVAL}.compleasm/summary.txt",
    params:
        outdir=lambda w, output: os.path.dirname(output.summary),
        lineage=config["compleasm"]["lineage"],
        library=config["compleasm"]["library"],
    log:
        f"{EVAL}.compleasm.log",
    threads: config["threads"]["compleasm"]
    conda:
        "../../envs/compleasm.yaml"
    shell:
        "compleasm run -a {input.fa} -o {params.outdir} -t {threads} "
        "-l {params.lineage} -L {params.library} &> {log}"


rule parse_compleasm_genome:
    """Fold the compleasm report into one standard BUSCO line."""
    input:
        summary=rules.compleasm_genome.output.summary,
    output:
        tsv=f"{EVAL}.compleasm.tsv",
    params:
        kind="Genome",
    script:
        "../../scripts/parse_compleasm.py"


rule assembly_stats:
    """Contiguity statistics (N50, L50, sizes) with BBMap stats.sh."""
    input:
        fa=ASSEMBLY,
    output:
        stats=f"{EVAL}.stats.txt",
    log:
        f"{EVAL}.stats.log",
    conda:
        "../../envs/assembly_stats.yaml"
    shell:
        "stats.sh in={input.fa} out={output.stats} addname=t extended=t "
        "format=3 -Xmx10g &> {log}"
