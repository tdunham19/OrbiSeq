process REMOVE_TRAILING_FASTA_NS {
    tag "$meta.id"
    label "no_publish"
    
    conda "${moduleDir}/environment.yml"
        
	input: 
    tuple val(meta), path(input)

	output: 
	tuple val(meta), path("*.fa") , emit: fa

	script: 
	"""
	remove_trailing_fasta_Ns ${input} > ${meta.id}_RM_Ns.fa
	"""
	}