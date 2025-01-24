process CREATE_MASK_FILE {
    tag "$meta.id"
    // label "no_publish"
    
    conda "${moduleDir}/environment.yml"
        
	input: 
    tuple val(meta), path(input)

	output: 
	tuple val(meta), path("*.mask") , emit: mask

	script: 
	"""
    create_mask_file ${input} > ${meta.id}.mask
	"""
	}