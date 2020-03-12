# Conifer genome long-read polishing workflow

The goal of this workflow is to enable solid and relatively fast long-read polishing of conifer genomes.

## Requirements

- snakemake
- conda

The recommended approach would be to create a snakemake conda environment that is used for running the workflow.

```sh
conda create -n snakemake snakemake
conda activate snakemake
```

## Installation

Either download and extract the latest release, or fork the repository and clone it in order to keep track of changes.
The latter approach is recommended.
If you fork the repository, then create a new branch where project specific changes are made.
Any changes that could be of more general value to the workflow can potentially be merged into this repository with a pull request.

If you choose to fork the repository, the following steps are recommended:

1. Fork the repository to your own account.
2. Create a separate branch where project-specific changes can be made.
3. Commit and push changes to the project-specific branch.
4. Any code that is more general and can be of value to the workflow itself can possibly be merged into the [upstream repository](https://github.com/maehler/conifer-long-read-polishing) with a [pull request](https://help.github.com/en/articles/creating-a-pull-request). This can be accomplished in a number of ways, but preferrably create a separate branch containing only these changes, and then do the pull request based on this.

## Configuration

In [`config.yaml`](config.yaml) there are a couple of entries that need special attention: `contig_fasta`,`read_metadata` and `iterations`.
`contig_fasta` should point to a FASTA file containing the contigs that should be polished with long reads.
`read_metadata` should point to a TSV file containing information on where to find the short read data. It should contain the following columns:

- `filename`: Path to the file containing the reads.
- `filetype`: Supported types are `fasta`, `fastq` and `bam`.
- `gzipped`: Whether the file is gzipped or not (boolean).
- `readtype`: The type of sequencing technology generated the reads, *e.g.* PacBio Sequel or Oxford Nanopore.
- `comment`: Any comments about the file (optional).

The parameter `iterations` controls how many iterations of long-read polishing should be performed.

### Cluster configuration

The workflow works best by running it in a cluster environment, and in order to minimise the amount of typing needed on the command line, some configuration is needed.
One of the more convenient ways of handling this is with a [Snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).
There are [existing templates](https://github.com/snakemake-profiles/doc) for the most common job schedulers available.
For convenience, a rule in the workflow is dedicated for setting up a profile, and it is handled by the following section in `config.yaml`:

```yaml
# Cluster configuration
cluster:
    # Cookiecutter parameters
    cookiecutter:
        url: https://github.com/Snakemake-Profiles/slurm
        profile_name: longread_polish
        cluster_config: cluster/slurm.yaml
        advanced_argument_conversion: false

    # Default Snakemake parameters
    snakemake:
        use-conda: true
        use-envmodules: true
        restart-times: 0
        jobs: 3000
        latency-wait: 120
```

This particular example sets up a profile for running jobs using Slurm.
The section `cookiecutter` defines the parameters that are specific for the cookiecutter template that is being used---in this case https://github.com/Snakemake-Profiles/slurm.
In the section `snakemake` some defaults for Snakemake are set up.
These parameters should be named as the long options for Snakemake.
In order to set up a profile, run

```sh
snakemake --use-conda cluster_config
```

The profile will be created under `$HOME/.config/snakemake`.

## Executing the workflow

It can be good to first run the workflow as a dry-run, *i.e.* just checking what rules would be executed but not actually execute them.
This is good in order to identify issues with the configuration and related files.

```sh
snakemake --profile <profile name> --dry-run
```

In order to run the full workflow from start to finish, run Snakemake with the profile generated

```sh
snakemake --profile <profile name>
```

## Output

The results from the workflow are stored in the directory `results`. The following subdirectories are created:

- `results/alignments`: All long read alignments generated.
- `results/racon_{iteration}`: Polished contigs from racon for each iteration.

## Acknowledgements

This workflow has been heavily inspired by the workflows in https://github.com/snakemake-workflows and the work of [Johannes Köster](https://github.com/johanneskoester) and others.
