process IDENTIFY_BEST_SEGMENTS_FROM_SAM {
    tag "$meta.id"
    // label "no_publish"
    
    conda "${moduleDir}/environment.yml"
        
	input: 
	tuple val(meta), path(input)
    tuple val(meta2), path(reference)
	
	output: 
	tuple val(meta), path("*.fa") , emit: fa

	script: 
	"""
	identify_best_segments_from_sam ${reference} ${input} > ${meta.id}_best10_refseq.fa
	"""
	}