include { NANOPORE_CONSENSUS } from './workflows/nanopore.nf'
include { ILLUMINA_CONSENSUS } from './workflows/illumina.nf'

println "Selected platform: ${params.platform}"
println "Selected reference: ${params.reference}"

workflow {  
  if (params.reference == 'BTV') {
    path_to_reference = params.btv_reference
  } else if (params.reference == 'EHD') {
    path_to_reference = params.ehdv_reference
  } else if (params.reference == 'custom') {
    path_to_reference = params.custom_reference
  } else {
    error "No reference provided. Please specify --reference (BTV, EHDV, or custom)."
  }

  if (params.platform == 'nanopore') {
      NANOPORE_CONSENSUS()
  } else if (params.platform == 'illumina') {
      ILLUMINA_CONSENSUS()
  } else {
      exit 1, "Error: Unknown platform specified. Use 'nanopore' or 'illumina'."
  }
}