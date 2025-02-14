process CD_HIT_EST {
  tag "$meta.id"
  label 'lowmem_threaded'

  conda (params.enable_conda ? 'bioconda::cd-hit=4.8.1' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/cd-hit:4.8.1--h5b5514e_7" 
  } else {
    container "quay.io/biocontainers/cd-hit:4.8.1--h5b5514e_7" 
  }

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path("*.uniq.fastq.gz") , emit: reads
  path "versions.yml"                      , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''
  def prefix           = task.ext.prefix ?: "${meta.id}"

  def r1 = reads[0].toString()
  def r2 = reads[1].toString()

  def r1_out = r1.replaceAll('.fastq.gz$', ".uniq.fastq")
  def r2_out = r2 ? r2.replaceAll('.fastq.gz$', ".uniq.fastq") : ""

  def paired_input    = reads[1] ? "-j $r2" : ""
  def paired_output   = reads[1] ? "-op $r2_out" : ""


  """
  cd-hit-est \
   -T ${task.cpus} \
   $args \
   -i $r1 \
   $paired_input \
   -o $r1_out \
   $paired_output 

  # cd-hit-est doesn't output compressed fastq 
  gzip *.uniq.fastq

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      cd-hit-est: 4.8.1
  END_VERSIONS
  """

}
