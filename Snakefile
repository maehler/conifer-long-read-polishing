import os
import pandas as pd
from pathlib import Path
import re
from snakemake.utils import validate, min_version

min_version('5.10.0')

configfile: 'config.yaml'
validate(config, 'schemas/config.schema.yaml')

read_metadata = pd.read_table(config['read_metadata']) \
    .set_index('filename', drop=False)
validate(read_metadata, 'schemas/read_metadata.schema.yaml')

localrules: all, bam_fofn, symlink_assembly, cluster_config,
            split_fasta, index_fasta,
            fasta_slice_fofn, racon_aggregate, polished_fofn

wildcard_constraints:
    iteration=r'\d+',
    slice=r'\d+'

rule all:
    input: 'data/contigs_racon_{iteration}.fasta'.format(iteration=config['iterations'])

include: 'rules/data_management.smk'
include: 'rules/polishing.smk'
include: 'rules/alignment.smk'
include: 'rules/cluster_config.smk'