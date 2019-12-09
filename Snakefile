import os
import pandas as pd
from pathlib import Path
import re
from snakemake.utils import validate

configfile: 'config.yaml'
validate(config, 'schemas/config.schema.yaml')

read_metadata = pd.read_table(config['read_metadata']) \
    .set_index('filename', drop=False)
validate(read_metadata, 'schemas/read_metadata.schema.yaml')

localrules: all, bam_fofn, symlink_assembly, cluster_config

rule all:
    input: 'results/arrow/aligned_subread_bams.fofn'

rule bam_fofn:
    '''
    Create a file of filenames (fofn) of all subreads
    aligned to the assembly.
    '''
    input: expand('results/alignments/{bam_name}_alignments.bam', \
                  bam_name=[Path(x).stem for x in read_metadata.filename[read_metadata.filetype == 'bam']])
    output: 'results/arrow/aligned_subread_bams.fofn'
    shell: 'printf "%s\\n" $(realpath -es {input}) > {output}'

def get_full_bam_path(wildcards):
    '''
    Get the full path to a bam file based on its basename.
    '''
    matches = [re.search(r'^.+{fname}.bam$'.format(fname=wildcards.bam_name), x) \
               for x in read_metadata.filename]
    matches = [x.group(0) for x in matches if x is not None]

    if len(matches) > 1:
        raise WorkflowError('multiple bam files with the name {fname}.bam' \
            .format(fname=wildcards.bam_name))

    return matches[0]

rule pbmm2_align:
    '''
    Align reads with pbmm2, i.e. the PacBio wrapper for minimap2.
    '''
    input:
        reference='reference/{reference}_subread.mmi' \
            .format(reference=Path(config['contig_fasta']).stem),
        query_bam=get_full_bam_path
    output:
        'results/alignments/{bam_name}_alignments.bam'
    threads: 10
    params:
        preset='SUBREAD',
    conda: 'envs/arrow.yaml'
    shell:
        '''
        pbmm2 align \
            --preset {params.preset} \
            --sort \
            -m 4G \
            -j {threads} \
            {input.reference} \
            {input.query_bam} \
            {output}
        '''

rule pbmm2_index:
    '''
    Create a minimap index for the assembly.
    '''
    input: 'reference/{reference}.fasta'
    output: 'reference/{reference}_subread.mmi'
    threads: 20
    params:
        preset='SUBREAD'
    conda: 'envs/arrow.yaml'
    shell:
        '''
        pbmm2 index \
            --num-threads {threads} \
            --preset {params.preset} \
            {input} \
            {output}
        '''

rule symlink_assembly:
    '''
    Symlink the assembly FASTA file to the project directory.
    '''
    input: Path(config['contig_fasta']).absolute()
    output: 'reference/{reference}.fasta' \
        .format(reference=Path(config['contig_fasta']).stem)
    shell:
        '''
        cd reference
        ln -s {input} $(basename {output})
        '''

rule cluster_config:
    '''
    Generate a cluster profile.
    '''
    output: directory( \
        '{home}/.config/snakemake/{profile_name}' \
            .format(home=Path.home(), \
                    profile_name=config['cluster']['profile_name']))
    params:
        url=config['cluster']['cookiecutter_url'],
        profile_name=config['cluster']['profile_name']
    conda: 'envs/cluster_config.yaml'
    shell: 'bash scripts/cluster_config.sh {params.url} {params.profile_name}'
