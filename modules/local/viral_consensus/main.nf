process VIRAL_CONSENSUS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'https://depot.galaxyproject.org/singularity/viral_consensus:0.0.6--h43eeafb_1' :
       'biocontainers/viral_consensus:0.0.6--h43eeafb_1' }"

    input:
    tuple val(meta), path(bam), path(ref_fasta)
    val(min_qual)
    val(min_depth)
    val(min_freq)

    output:
    tuple val(meta), path("*.consensus.fa")         , emit: fasta
    tuple val(meta), path("*.position_counts.txt")  , emit: position_counts
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_base = ref_fasta.baseName
    """
    viral_consensus \
        -i $bam \
        -r $ref_fasta \
        -o ${prefix}.${ref_base}.consensus.fa \
        -op ${prefix}.${ref_base}.position_counts.txt \
        -q $min_qual \
        -d $min_depth \
        -f $min_freq \
        $args 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(echo \$(viral_consensus --version 2>&1) | sed 's/^viral_consensus /viral_consensus: /')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_base = ref_fasta.baseName
    """
    echo "" > ${prefix}.${ref_base}.consensus.fa 
    echo "" > ${prefix}.${ref_base}.position_counts.txtx 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(echo \$(viral_consensus --version 2>&1) | sed 's/^viral_consensus /viral_consensus: /')
    END_VERSIONS
    """
}

