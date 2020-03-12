def get_racon_aggregate_input(wildcards):
    checkpoint_output = checkpoints.split_fasta \
        .get(iteration=int(wildcards.iteration)-1) \
        .output[0]
    fname_pattern = '{output_dir}/contigs_racon_{iteration}_{{slice}}.fasta' \
        .format(output_dir=checkpoint_output, iteration=int(wildcards.iteration)-1)
    gwc = glob_wildcards(fname_pattern)
    fnames = expand('results/racon_{iteration}/polished_slices/polished_slice_{slice}.fasta',
        iteration=wildcards.iteration, slice=gwc.slice)
    return sorted(fnames, key=lambda x: int(re.search(r'(\d+)\.fasta$', x).group(1)))

rule racon_aggregate:
    """
    Merge the polished FASTA slices from racon.
    """
    input: get_racon_aggregate_input
    output:
        fasta=protected('results/racon_{iteration}/contigs_racon_{iteration}.fasta'),
        linked_fasta='data/contigs_racon_{iteration}.fasta'
    priority: 10
    run:
        with open(output['fasta'], 'w') as f:
            print('\n'.join(str(Path(x).resolve()) for x in input), file=f)

rule racon:
    """
    Run one round of polishing using racon.
    """
    input:
        fastafiles=lambda wildcards: 'data/contigs_racon_{iteration}_slices.fofn' \
            .format(iteration=int(wildcards.iteration)-1),
        sam=lambda wildcards: 'results/alignments/alignment_slices_{iteration}/subread_alignments_slice_{{slice}}.sam.gz' \
            .format(iteration=int(wildcards.iteration)-1)
    output:
        fasta='results/racon_{iteration}/polished_slices/polished_slice_{slice}.fasta'
    threads: 10
    conda: '../envs/racon.yaml'
    shell:
        """
        slice_fasta=$(awk 'NR == {wildcards.slice} + 1' {input.fastafiles})

        samtools fastq {input.sam} | \\
            pigz -c -p {threads} > {input.sam}.fq.gz

        racon --include-unpolished --threads {threads} \\
            {input.sam}.fq.gz \\
            {input.sam} \\
            ${{slice_fasta}} > {output.fasta}
        """
