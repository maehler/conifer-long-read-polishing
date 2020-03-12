rule symlink_assembly:
    """
    Symlink the assembly FASTA file to the project directory.
    """
    input: Path(config['contig_fasta']).absolute()
    output: 'data/contigs_racon_0.fasta'
    priority: 50
    shell:
        """
        cd $(dirname {output})
        ln -s {input} $(basename {output})
        """

rule index_fasta:
    input: '{filename}.{extension}'
    output: '{filename}.{extension}.fai'
    wildcard_constraints:
        extension=r'fa|fasta'
    conda: '../envs/samtools.yaml'
    envmodules: 'bioinfo-tools', 'samtools/1.10'
    shell: 'samtools faidx {input}'

def get_fasta_slices(wildcards):
    checkpoint_output = checkpoints.split_fasta.get(iteration=wildcards.iteration).output[0]
    fname_pattern = '{output_dir}/contigs_racon_{iteration}_{{slice}}.fasta' \
        .format(output_dir=checkpoint_output, iteration=wildcards.iteration)
    gwc = glob_wildcards(fname_pattern)
    raw_files = expand(fname_pattern, slice=gwc.slice)
    files = sorted(raw_files, key=lambda x: int(re.search(r'(\d+)\.fasta$', x).group(1)))
    return files

rule fasta_slice_fofn:
    input: get_fasta_slices
    output: 'data/contigs_racon_{iteration}_slices.fofn'
    run:
        with open(output[0], 'w') as f:
            for fname in input:
                print(Path(fname).absolute(), file=f)

checkpoint split_fasta:
    input:
        fasta='data/contigs_racon_{iteration}.fasta'
    output:
        directory('data/contigs_racon_{iteration}_slices')
    params:
        chunk_size=config['slice-size']
    conda: '../envs/racon.yaml'
    shell:
        """
        mkdir {output}
        rampler split -o {output} {input.fasta} {params.chunk_size}
        """

rule merge_minimap_bams:
    input: 'results/alignments/aligned_subread_bams_{iteration}.fofn'
    output:
        bam='results/alignments/all_subread_alignments_{iteration}.bam'
    conda: '../envs/samtools.yaml'
    envmodules: 'bioinfo-tools', 'samtools/1.10'
    threads: 20
    shell:
        """
        samtools merge -@ $(({threads} - 1)) --write-index -b {input} {output.bam}
        """

rule bam_slice:
    """
    Get a slice of the alignments to use for polishing with racon.
    """
    input:
        fastafiles='data/contigs_racon_{iteration}_slices.fofn',
        bam='results/alignments/aligned_subread_bams_{iteration}.fofn'
    output:
        sam='results/alignments/alignment_slices_{iteration}/subread_alignments_slice_{slice}.sam.gz',
    wildcard_constraints:
        run=r'\d+'
    threads: 1
    conda: '../envs/samtools.yaml'
    envmodules: 'bioinfo-tools', 'samtools/1.10'
    shell:
        """
        slice_fasta=$(awk 'NR == {wildcards.slice} + 1' {input.fastafiles})
        samtools faidx ${{slice_fasta}}

        tmpdir="$(dirname {output.sam})/slice_{wildcards.slice}_tmp"
        mkdir -p ${{tmpdir}}

        region_bed="${{tmpdir}}/slice_{wildcards.slice}.bed"
        awk 'BEGIN {{OFS="\\t"}} {{print $1, 0, $2}}' ${{slice_fasta}}.fai > ${{region_bed}}

        while read -r bamfile; do
            echo ${{bamfile}}
            samtools view -@$(({threads} - 1)) -ML ${{region_bed}} -o ${{tmpdir}}/$(basename ${{bamfile}}).slice.bam ${{bamfile}}
        done < {input.bam}

        bam_fofn=${{tmpdir}}/alignment_bams.fofn
        find ${{tmpdir}} -type f -name "*.slice.bam" > ${{bam_fofn}}

        samtools merge -f -@$(({threads} - 1)) -b ${{bam_fofn}} ${{tmpdir}}/merged.bam
        samtools sort -@$(({threads} - 1)) -m 4G -n -o {output.sam} ${{tmpdir}}/merged.bam

        rm -r ${{tmpdir}}
        """
