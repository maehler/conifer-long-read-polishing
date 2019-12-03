from snakemake.utils import validate

configfile: 'config.yaml'
validate(config, 'schemas/config.schema.yaml')

localrules: all

rule all:
    input: config['contig_fasta']
