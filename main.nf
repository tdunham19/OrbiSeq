include { NANOPORE_CONSENSUS } from './workflows/nanopore.nf'
include { ILLUMINA_CONSENSUS } from './workflows/illumina.nf'

println "Selected platform: ${params.platform}"
println "Selected reference: ${params.reference}"

workflow {  
  
  if (params.reference == 'BTV') {
     params.BTV_reference
  } else if (params.reference == 'EHDV') {
     params.EHDV_reference
  } else if (params.reference == 'custom') {
     params.custom_reference
  } else {
    exit 1, "Error: Unknown reference provided. Use --reference 'BTV', 'EHDV', or 'custom'."
  }
  
  if (params.platform == 'nanopore') {
      NANOPORE_CONSENSUS()
  } else if (params.platform == 'illumina') {
      ILLUMINA_CONSENSUS()
  } else {
      exit 1, "Error: Unknown platform specified. Use --platform 'nanopore' or 'illumina'."
  }
  
}