/*
  This is just a dummy process to enable saving of output files         
  that aren't otherwise output.  For instance, files created by
  collectFile() in a workflow context.
 */
process SAVE_OUTPUT_FILES {
    tag "$file_to_save"
    label 'process_low'

    input:
    path(files_to_save, stageAs: "input/*")

    output:
    path("*", includeInputs: true) 

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    // def link_name  = file_to_save.fileName.name

	 def complete_command    = ""

    // handle multiple files
	 // construct a complete command consisting of a ln command for each file
	 for( file in file_to_save){
      def link_name  = file.fileName.name
		// use -L to force a hard link
		complete_command = ${complete_command} + "ln -L $file $link_name; "

    }

	 // this will be the complete command, as last expression in this process block
	 command
}
