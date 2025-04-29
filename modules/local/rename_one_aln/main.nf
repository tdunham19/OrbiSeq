process RENAME_ONE_ALN {
    tag "$meta.id"
    label "no_publish"
    
    conda "${moduleDir}/environment.yml"
        
	input: 
    tuple val(meta), path(refseq_fasta), path(new_aln)
    val (suffix) // filename suffix

	output: 
	tuple val(meta), path("*.sam") , emit: sam
	tuple val(meta), path("*.bam") , emit: bam

	script: 
	"""
	rename_one_aln $meta.id $refseq_fasta $new_aln > ${new_aln}.renamed.${suffix}
	"""
	}