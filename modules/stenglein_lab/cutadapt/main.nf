/*
 This module is derived from the nf-core fastqc module:
 https://github.com/nf-core/modules/tree/master/modules/nf-core/fastqc
 
 I am modifying it because of two issues I ran into when trying to use it 
 in pipelines I was developing.  

 **Issue 1: This process emits trimmed files using path('*.trim.fastq.gz')
 
 The issue is that I'd like this pipeline to work seamlessly for either 
 single-end or paired-end data (or a mixture of the 2)..  

 path() behaves differently if there is one file (single-end case) or >1 file 
 (paired-end case).  When there is one file, it outputs a plain value.  When
 there is > 1 file, it emits a list.  This means downstream handling of this
 output can be challenging.

 This inconsistency is discussed in this issue: 
 https://github.com/nextflow-io/nextflow/issues/2425
 
 The solution was to create a files() method, which always returns a list.
 The corresponding commit is:
https://github.com/nextflow-io/nextflow/commit/00bb8896a18d23bf82df66b267b4d03f0a42104d
 
 So the main change I made was to change path() -> files() in the reads output
 This causes the process to always output a list, even if the list only has
 one element..

 A more backwards-compatible solution would be to declare a second output
 channel (reads_list or something) that used files() instead of path. 

 ** Issue 2: This process doesn't dynamically handle either single-end or paired-
 end data.  

 It uses ext.args to provide options to cutadapt.  The problem is 
 that cutadapt doesn't by default apply single-end trimming parameters to 
 paired-end data.  So if you only provided single-end trimming options (e.g. -a)
 then the paired reads wouldn't be trimmed.  And if you provided paired-end 
 options (e.g. -A) when there was only single-end data, cutadapt errors and quits.  

 A possible solution is to use args2 to provide paired-end trimming options
 and then only use args2 when meta.single_end is false. 

 This work-around is perhaps not ideal, since arg2 is documented as being for 
 the second tool in multi-tool processes, but it seems like this is a second 
 case of a single-tool process, and so more or less lives up to the spirit of
 what args2 is for (extra options when needed).

*/

process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::cutadapt=3.4' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:3.4--py39h38f01e4_1' :
        'quay.io/biocontainers/cutadapt:3.4--py39h38f01e4_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fastq.gz') , emit: reads
    // reads_list will always emit a list of trimmed files
    // even if only one file
    // tuple val(meta), files('*.trim.fastq.gz'), emit: reads_list 
    tuple val(meta), path('*.log')           , emit: log
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args ?: ''
    // use arg2 as PE trimming options
    def paired_end_args  = meta.single_end ? '': task.ext.args2 
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def trimmed          = meta.single_end ? "-o ${prefix}.trim.fastq.gz" : "-o ${prefix}_1.trim.fastq.gz -p ${prefix}_2.trim.fastq.gz"
    """
    cutadapt \\
        --cores $task.cpus \\
        $args \\
        $paired_end_args \\
        $trimmed \\
        $reads \\
        > ${prefix}.cutadapt.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
