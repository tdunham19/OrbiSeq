process COLLAPSE_DUPLICATE_READS {
  label 'lowmem_threaded'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/cd-hit-auxtools:4.8.1--h7d875b9_1"
  } else {
    container "quay.io/biocontainers/cd-hit-auxtools:4.8.1--h7d875b9_1"
  }

  input:
  tuple val(meta), path(input_fastq) 

  output:
  tuple val(meta), path("*.uniq.fastq.gz") , emit: reads
  path "versions.yml"                      , emit: versions

  script:

  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"


  // this handles paired-end data, in which case must specify a paired output file
  def r1_gz           = input_fastq[0]
  def r2_gz           = input_fastq[1] ? input_fastq[1] : ""
  def r1              = input_fastq[0].getBaseName()
  def r2              = input_fastq[1] ? input_fastq[1].getBaseName() : ""
  def single_input    = "-i $r1"
  def single_output   = input_fastq[1] ? "-o ${prefix}_1.trim.uniq.fastq" : "-o ${prefix}.trim.uniq.fastq"
  def paired_input    = input_fastq[1] ? "-i2 $r2" : ""
  def paired_output   = input_fastq[1] ? "-o2 ${prefix}_2.trim.uniq.fastq" : ""

  """
  # decompress input because cd-hit-dup can't handle compressed input :|
  # gzip in cd-hit-auxtools singularity image
  gunzip $r1_gz $r2_gz 

  cd-hit-dup \
   $args \
   $single_input \
   $paired_input \
   $single_output \
   $paired_output \

  # recompress everything
  gzip *.fastq 

  # get rid of cd-hit .clstr output files: big and not needed
  rm *.clstr

  cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       cd-hit-dup: 4.8.1
  END_VERSIONS
  """
}
