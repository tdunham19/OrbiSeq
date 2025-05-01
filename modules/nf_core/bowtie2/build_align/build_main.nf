process BOWTIE2_BUILD {
    tag "$meta.id"
    label "no_publish"    

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.5.2--py39h6fed5c7_0' :
        'biocontainers/bowtie2:2.5.2--py39h6fed5c7_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path('bowtie2')                 , emit: index
    path  "versions.yml"                             , emit: versions

    script:
    """
    mkdir bowtie2
    bowtie2-build $input bowtie2/${input.baseName}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
