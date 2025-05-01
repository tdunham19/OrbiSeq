process RENAME_ONE_ALN {
    tag "$meta.id"
    label "no_publish"
    
    conda "${moduleDir}/environment.yml"
	container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"
        
	input: 
    tuple val(meta), path(bam_ref), path(new_bam)
    val (suffix) // filename suffix

	output: 
	tuple val(meta), path("*.bam") , emit: bam

	script: 
	"""
	rename_one_aln $meta.id $bam_ref $new_bam ${suffix}
	"""
	}