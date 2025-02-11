process BCFTOOLS_INDEX {
    tag "$meta.id"
    // label "no_publish"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.csi"), optional:true, emit: csi
    tuple val(meta), path("*.gz"), optional:true, emit: gz
    tuple val(meta), path("*.gz", includeInputs: true), path("*.csi"), emit: gz_and_csi

    script:
    """
    bcftools index ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}