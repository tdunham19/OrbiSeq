process MINIMAP2_ALIGN_TO_EXISTING {
	tag "$meta.id"
    label "no_publish"

conda "${moduleDir}/environment.yml"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

input:
tuple val(meta), path(reads)
tuple val(meta2), path(reference)
val (suffix) // filename suffix

output:
  tuple val(meta), path("*.sam")                       , emit: sam
  tuple val(meta), path("*.bam")                       , emit: bam
  tuple val(meta), path(reference)                     , emit: reference
  tuple val(meta), path("*.bam"), path(reference)      , emit: bam_reference
  path ("versions.yml")                                , emit: versions

  script:
  """
  minimap2 -ax map-ont $reference $reads | samtools view -h -F 4 > ${meta.id}.${suffix}.sam 
  samtools sort -O bam ${meta.id}.${suffix}.sam  > ${meta.id}.${suffix}.bam

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
  """
}

process MINIMAP2_ALIGN {
	tag "$meta.id"
	label "no_publish"

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

  input:
  tuple val(meta),  path(reads), path(refseq)
  val (suffix) // filename suffix

  output:
  tuple val(meta), path("*.sam")                       , emit: sam
  tuple val(meta), path("*.bam")                       , emit: bam
  tuple val(meta), path(refseq)                        , emit: refseq
  tuple val(meta), path("*.bam"), path(refseq)         , emit: bam_refseq
  path ("versions.yml")                                , emit: versions

  script:
  """
  minimap2 -ax map-ont $refseq $reads | samtools view -h -F 4 > ${meta.id}_${suffix}.sam 
  samtools sort -O bam ${meta.id}_${suffix}.sam  > ${meta.id}_${suffix}.bam

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
  """
}
