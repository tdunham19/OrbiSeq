process PROCESS_FASTQ_COUNTS {
    tag "$fastq_counts_tsv"
    label 'process_low'

    conda (params.enable_conda ? 'conda-forge::r-tidyverse=1.3.1' : null) 
    // why aren't singularity biocontainers updated to a newer tidyverse version?
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1' : 
        'quay.io/biocontainers/r-tidyverse:1.2.1' }"

    input:
    path(fastq_counts_tsv)                                                

    output:
    path "all_read_counts.txt" ,         emit: all_read_counts
    path "versions.yml" ,                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # make a copy of input file so can be output to results dir
    cp $fastq_counts_tsv all_fastq_counts_tidy.txt

    # make a slightly untidy table too
    process_fastq_counts.R $fastq_counts_tsv 
    
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        R: \$(echo \$(R --version 2>&1))
    END_VERSIONS
    """
}
