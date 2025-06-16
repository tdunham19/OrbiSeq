process REMOVE_TRAILING_FASTA_NS {
   tag "$meta.id"
   label "no_publish"
    
   // conda "${moduleDir}/environment.yml"
        
	input: 
   tuple val(meta), path(input)

	output: 
	tuple val(meta), path(input) , emit: fa

	script: 
	"""
	mv $input ${input}.with_trailing_N
	remove_trailing_fasta_Ns ${input}.with_trailing_N > ${input}
	"""
	}
