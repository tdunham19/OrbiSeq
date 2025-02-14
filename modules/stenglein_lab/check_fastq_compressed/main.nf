process CHECK_FASTQ_COMPRESSED {
    tag "$meta.id"
    label 'process_low'

    // this is all groovy so no containers needed

    input:
    tuple val(meta), path(files)                                                

    output:
    // path "versions.yml" ,             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def r1 = files[0]
    def r2 = files[1]

    // simply check that file names end with .gz
    r1_not_gz   = !(r1.toString().endsWith('.gz'))
    r2_not_gz   = !(r2 && r2.toString().endsWith('.gz'))

    single_end_command = r1_not_gz ? """
        >&2 echo "Error: Input file $r1 is not gzip compressed."
        >&2 echo "       Please use gzip to compress all input fastq."
        exit 2
    """ : ""

    paired_end_command = r2 && r2_not_gz ? """
        >&2 echo "Error: Input file $r2 is not gzip compressed."
        >&2 echo "       Please use gzip to compress all input fastq."
        exit 2
    """ : ""

    // since this process just groovy, version info 
    // already included in nextflow version info
    def version_command = """
    """

    // this will output the appropriate combination of bash commands
    if (meta.single_end) {                                                      
      single_end_command + version_command
    } else {
      single_end_command + paired_end_command + version_command
    }
}
