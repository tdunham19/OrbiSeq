/*
  Concatenate files 
 */
process CONCATENATE_FILES {
    tag "$file_to_save"
    label 'process_low'
    label "no_publish"

    input:
	 tuple val(meta), path(files)
	 val (filename_suffix)

    output:
    tuple val (meta), path("*${filename_suffix}"), emit: file

    when:
    task.ext.when == null || task.ext.when

    script:
	 def output_file = "${meta.id}${filename_suffix}"
    """
	 cat $files > $output_file
    """

}
