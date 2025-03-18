process LOFREQ_CALL {
    tag "$meta.id"
	// label "no_publish"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py38h588ecb2_4' :
        'biocontainers/lofreq:2.1.5--py38h588ecb2_4' }"

    input:
    tuple val(meta), path(input), path(fasta)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"              , emit: versions

    script:
    """
    lofreq call -f ${fasta} -o ${meta.id}.vcf ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """
}
