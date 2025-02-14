/*
  This is just a dummy process to enable saving of output files         
  that aren't otherwise output.  For instance, files created by
  collectFile() in a workflow context.
 */
process SAVE_OUTPUT_FILE {
    tag "$file_to_save"
    label 'process_low'

    input:
    path(file_to_save, stageAs: "input/*")

    output:
    path("*", includeInputs: true) 

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
	 def link_name  = file_to_save.fileName.name
    """
    # force a hard link with ln -L
    ln -L $file_to_save $link_name
    """

}
