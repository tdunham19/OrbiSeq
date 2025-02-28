include { NANOPORE_CONSENSUS } from './workflows/nanopore.nf'
include { ILLUMINA_CONSENSUS } from './workflows/illumina.nf'

println "Selected platform: ${params.platform}"
println "Selected reference: ${params.reference}"

println "fastq directory: ${params.fastq_dir}"
println "Output directory: ${params.outdir}"

workflow {  
  
  if (params.platform == 'nanopore') {
      NANOPORE_CONSENSUS()
  } else if (params.platform == 'illumina') {
      ILLUMINA_CONSENSUS()
  } else {
      exit 1, "Error: Unknown platform specified. Use --platform 'nanopore' or 'illumina'."
  }
  
}