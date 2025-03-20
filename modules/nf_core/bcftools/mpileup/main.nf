process BCFTOOLS_MPILEUP {
    tag "$meta.id"
	// label "no_publish"
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(input), path(best10refseq)

    output:
    tuple val(meta), path("*vcf")     , emit: vcf
    
    script:
    """
    bcftools mpileup -f ${best10refseq} -Ou -o ${meta.id}.vcf ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
