rule bam_fofn:
    """
    Create a file of filenames (fofn) of all subreads
    aligned to the assembly.
    """
    input: expand('results/alignments/subread_alignments_{{iteration}}/{bam_name}_alignments.bam', \
                  bam_name=[Path(x).stem for x in read_metadata.filename[read_metadata.filetype == 'bam']])
    output: 'results/alignments/aligned_subread_bams_{iteration}.fofn'
    shell: 'printf "%s\\n" $(realpath -es {input}) > {output}'

def get_full_bam_path(wildcards):
    """
    Get the full path to a bam file based on its basename.
    """
    matches = [re.search(r'^.+{fname}.bam$'.format(fname=wildcards.bam_name), x) \
               for x in read_metadata.filename]
    matches = [x.group(0) for x in matches if x is not None]

    if len(matches) > 1:
        raise WorkflowError('multiple bam files with the name {fname}.bam' \
            .format(fname=wildcards.bam_name))

    return matches[0]

rule pbmm2_align:
    """
    Align reads with pbmm2, i.e. the PacBio wrapper for minimap2.
    """
    input:
        reference='reference/contigs_racon_{iteration}_subread.mmi',
        query_bam=get_full_bam_path
    output: 'results/alignments/subread_alignments_{iteration}/{bam_name}_alignments.bam'
    threads: 10
    params:
        preset='SUBREAD',
    conda: '../envs/arrow.yaml'
    shell:
        """
        pbmm2 align \\
            --preset {params.preset} \\
            --sort \\
            -m 4G \\
            -j {threads} \\
            {input.reference} \\
            {input.query_bam} \\
            {output}
        """

rule pbmm2_index:
    """
    Create a minimap index for the assembly.
    """
    input: 'data/contigs_racon_{iteration}.fasta'
    output: 'reference/contigs_racon_{iteration}_subread.mmi'
    threads: 20
    params:
        preset='SUBREAD'
    conda: '../envs/arrow.yaml'
    shell:
        """
        pbmm2 index \\
            --num-threads {threads} \\
            --preset {params.preset} \\
            {input} \\
            {output}
        """