process BOWTIE2_ALIGN {
    tag "$meta.id"
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0' :
        'biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0' }"

    input:
    tuple val(meta) , path(reads), path(index)
    
    output:
    tuple val(meta), path("*.sam")      , emit: sam
    path  "versions.yml"                , emit: versions

    script:
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/.rev.1.bt2//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/.rev.1.bt2l//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1
    
    bowtie2 -x \$INDEX -q -1 ${reads} -2 ${reads} --local --qc-filter --score-min C,120,1 --maxins 700 --time --no-unal --al-conc ${meta.id}.conc_hits.fastq --threads 24 -S ${meta.id}_new_draft_seq.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}