process COUNT_FASTQ {
    tag "$meta.id"
    label 'process_low'
    label 'no_publish'

    // use pigz to decompress fastq & pipe to wc
    conda (params.enable_conda ? "conda-forge::pigz=2.3.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'quay.io/biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(files), val(count_type)

    output:
    path("*.count.txt") ,  emit: count_file
    path "versions.yml" ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    def prefix = task.ext.prefix ?: "${meta.id}"

    '''
    pigz -dcf !{files[0]} | wc -l | awk '{print "!{meta.id}" "\t" "!{count_type}" "\t" $0}' > "!{meta.id}.!{count_type}.count.txt"

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
            pigz: $( pigz --version | sed -e "s/pigz //g" )
    END_VERSIONS
    '''
}
