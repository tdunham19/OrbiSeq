process RENAME_ONE_FASTA {
    tag "$meta.id"
    // label "no_publish"
    
    conda "${moduleDir}/environment.yml"
        
	input: 
    tuple val(meta), path(new_fasta), path(refseq_fasta)

	output: 
	tuple val(meta), path("*.fasta") , emit: fasta

	script: 
	"""
	rename_one_fasta $meta.id $refseq_fasta $new_fasta > ${new_fasta}.renamed.fasta
	"""
	}