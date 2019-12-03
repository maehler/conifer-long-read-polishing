# Conifer genome polishing workflow

The goal of this workflow is to enable solid polishing of conifer genome assemblies.

## Requirements

- snakemake
- conda ([miniconda](https://docs.conda.io/en/latest/miniconda.html))
- [cookiecutter](https://cookiecutter.readthedocs.io/en/latest/) (optional)

The recommended approach would be to create a snakemake conda environment that is used for running the workflow.

```sh
conda create -n snakemake snakemake
conda activate snakemake
```

