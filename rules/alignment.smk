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

rule minimap2_align:
    """
    Align subreads to the assembly with minimap2.

    The default parameters here are based on the SUBREAD preset for pbmm2.
    """
    input:
        reference='reference/contigs_racon_{iteration}.mmi',
        query_bam=get_full_bam_path
    output:
        bam='results/alignments/subread_alignments_{iteration}/{bam_name}_alignments.bam',
        tmpbam=temp('results/alignments/subread_alignments_{iteration}/{bam_name}_alignments.tmp.bam'),
        fasta=temp('results/alignments/subread_alignments_{iteration}/{bam_name}.fasta')
    threads: 10
    conda: '../envs/minimap2.yaml'
    envmodules: 'bioinfo-tools', 'minimap2/2.16', 'samtools/1.10'
    shell:
        """
        samtools fasta {input.query_bam} > {output.fasta}
        minimap2 \\
            -k 19 \\
            -w 10 \\
            -O 5,56 \\
            -E 4,1 \\
            -A 2 \\
            -B 5 \\
            -z 400,50 \\
            -r 2000 \\
            --lj-min-ratio 0.5 \\
            -g 5000 \\
            -a \\
            --secondary=no \\
            -t {threads} \\
            {input.reference} {output.fasta} | \\
            samtools view -hb -o {output.tmpbam}
        samtools sort -m 4G -@ $(({threads} - 1)) -o {output.bam} {output.tmpbam}
        """

rule minimap2_index:
    """
    Create a minimap2 index for the assembly.
    """
    input: 'data/contigs_racon_{iteration}.fasta'
    output: 'reference/contigs_racon_{iteration}.mmi'
    threads: 20
    conda: '../envs/minimap2.yaml'
    envmodules: 'bioinfo-tools', 'minimap2/2.16'
    shell:
        """
        minimap2 \\
            -I 500G \\
            -H \\
            -k 19 \\
            -w 10 \\
            -t {threads} \\
            -d {output} \\
            {input}
        """
