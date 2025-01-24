include { NANOPORE_CONSENSUS } from './workflows/nanopore.nf'
include { ILLUMINA_CONSENSUS } from './workflows/illumina.nf'

println "Selected platform: ${params.platform}"

workflow {
  if (params.platform == 'nanopore') {
      NANOPORE_CONSENSUS()
  } else if (params.platform == 'illumina') {
      ILLUMINA_CONSENSUS()
  } else {
      exit 1, "Error: Unknown platform specified. Use 'nanopore' or 'illumina'."
  }
}