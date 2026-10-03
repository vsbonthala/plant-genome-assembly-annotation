# Container for running the workflows with Snakemake.
# Tools are installed per rule from envs/*.yaml when run with --use-conda.
FROM condaforge/miniforge3:latest

RUN mamba install -y -c conda-forge -c bioconda \
        snakemake snakemake-executor-plugin-slurm \
    && mamba clean -afy

WORKDIR /work
COPY . /opt/plant-genome-workflows

ENTRYPOINT ["snakemake"]
CMD ["--help"]
