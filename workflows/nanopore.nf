#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include modules: local, nf-core, and Stenglein lab
include { PYCOQC										 	 	  } from '../modules/nf_core/pycoqc/main.nf'
include { NANOPLOT										 	 	  } from '../modules/nf_core/nanoplot/main.nf'
include { MINIMAP2_ALIGN_TO_EXISTING  	 				     	  } from '../modules/nf_core/minimap2/align/main.nf'
include { IDENTIFY_BEST_SEGMENTS_FROM_SAM     				 	  } from '../modules/local/identify_best_segments_from_sam/main.nf'
include { CALL_INDIVIDUAL_CONSENSUS_NANOPORE              	 	  } from '../subworkflows/call_individual_consensus_nanopore.nf'
include { RENAME_ONE_ALN									 	  } from '../modules/local/rename_one_aln/main.nf'
include { RENAME_ONE_FASTA as RENAME_ONE_FASTA_VC			 	  } from '../modules/local/rename_one_fasta/main.nf'
include { CONCATENATE_FILES as CONCATENATE_VC_FILES         	  } from '../modules/stenglein_lab/concatenate_files/main.nf'
include { REMOVE_TRAILING_FASTA_NS as REMOVE_TRAILING_FASTA_NS_VC } from '../modules/local/remove_trailing_fasta_ns/main.nf'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_FINAL_VC		  		  } from '../modules/nf_core/minimap2/align/main.nf'


workflow NANOPORE_CONSENSUS {

  ch_versions = Channel.empty()                                               

  // fastq input files

  Channel.fromFilePairs("${params.fastq_dir}/${params.nanopore_input_pattern}", size: -1, checkIfExists: true, maxDepth: 1)
  .map{ name, reads ->
         def matcher = name =~ /^([^.]+)/
         def meta = [:]
         if (matcher.find()) {
             meta.id = matcher.group(0)
         } else {
             meta.id = "UNKNOWN" 
         }
         [ meta, reads ]
     }
  .set { ch_reads }
  
  // refseq input files

    Channel.fromPath("${params.reference}")
    .collect()
    .map { reference ->
            def meta2 = [:]
            meta2.id = "reference"
            [meta2, reference]
        }
    .set { ch_reference }
    
    // summary input file
    
    Channel.fromPath("${params.summary_file}")
    .collect()
    .map { summary ->
            def meta3 = [:]
            meta3.id = "summary"
            [meta3, summary]
        }
    .set { ch_summary }
    
  // run PycoQC on sequencing summary file
  PYCOQC ( ch_summary )
  
  // run Nanoplot on input fastq files
  NANOPLOT ( ch_reads )
    
  // align input reads using minimap2
  MINIMAP2_ALIGN_TO_EXISTING ( ch_reads, ch_reference, "existing_refseq" )
  
  // extract new fasta file containing best aligned-to seqs for this dataset
  IDENTIFY_BEST_SEGMENTS_FROM_SAM ( MINIMAP2_ALIGN_TO_EXISTING.out.sam, ch_reference )
  
  // split up best10 segments into individual sequences because virus-focused 
  // consensus callers (namely iVar and viral_consensus) only work on one
  // sequence at a time
  // see: https://www.nextflow.io/docs/latest/reference/operator.html#splitfasta  
  IDENTIFY_BEST_SEGMENTS_FROM_SAM.out.fa
    .splitFasta(by: 1, file: true, elem: 1)
    .set { ch_best10_individual_fasta }
  
  // this uses the nextflow combine operator to create a new channel
  // that contains the reads for each dataset and all individual fasta files 
  individual_fasta_ch = ch_reads.combine(ch_best10_individual_fasta, by: 0)
  
  // parameters related consensus calling: min depth, basecall quality, frequency for consensus calling
  min_depth_ch = Channel.value(params.nanopore_min_depth)
  min_qual_ch  = Channel.value(params.nanopore_min_qual)
  min_freq_ch  = Channel.value(params.nanopore_min_freq)
  CALL_INDIVIDUAL_CONSENSUS_NANOPORE(individual_fasta_ch, min_depth_ch, min_qual_ch, min_freq_ch)
  
  // rename best10 alignments with unique id and segment number
  RENAME_ONE_ALN ( CALL_INDIVIDUAL_CONSENSUS_NANOPORE.out.minimap2_aln_bam_ref, "_best10_alignment")
  
  // rename file headers with unique id and segment number
  RENAME_ONE_FASTA_VC ( CALL_INDIVIDUAL_CONSENSUS_NANOPORE.out.viral_consensus_refseq_and_new, "_vc")

  // collect individual consensus sequences and combine into single files  
  collected_vc_fasta_ch   = RENAME_ONE_FASTA_VC.out.fasta.groupTuple()
  CONCATENATE_VC_FILES  (collected_vc_fasta_ch,   "_viral_consensus.fasta")
  
  // pipe output through remove_trailing_fasta_Ns to strip N characters from beginning and ends of seqs
  REMOVE_TRAILING_FASTA_NS_VC ( CONCATENATE_VC_FILES.out.file )
  
  // filter out empty fasta files
  FINAL_CONSENSUS_SEQUENCE_FILTERED_VC_NANOPORE = REMOVE_TRAILING_FASTA_NS_VC.out.filter { meta, fasta ->
    fasta.text.readLines().find { it && !it.startsWith(">") } != null 
    }

  // re-align data against the new draft sequence (ie. final consensus sequence) using minimap2.
  MINIMAP2_ALIGN_FINAL_VC ( ch_reads.join(FINAL_CONSENSUS_SEQUENCE_FILTERED_VC_NANOPORE), "viral_consensus")

  }
  
  // specify the entry point for the workflow
workflow {
  NANOPORE_CONSENSUS()
}
