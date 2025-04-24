#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include modules: local, nf-core, and Stenglein lab
include { PYCOQC										 	 } from '../modules/nf_core/pycoqc/main.nf'
include { MINIMAP2_ALIGN_TO_EXISTING  	 				     } from '../modules/nf_core/minimap2/align/main.nf'
include { IDENTIFY_BEST_SEGMENTS_FROM_SAM     				 } from '../modules/local/identify_best_segments_from_sam/main.nf'
include { MINIMAP2_ALIGN_TO_BEST10     						 } from '../modules/nf_core/minimap2/align/main.nf'

include { CALL_INDIVIDUAL_CONSENSUS_NANOPORE              	 } from '../subworkflows/call_individual_consensus_nanopore.nf'

include { RENAME_ONE_FASTA 									 } from '../modules/local/rename_one_fasta/main.nf'
include { CONCATENATE_FILES as CONCATENATE_VC_FILES          } from '../modules/stenglein_lab/concatenate_files/main.nf'
include { CONCATENATE_FILES as CONCATENATE_IVAR_FILES        } from '../modules/stenglein_lab/concatenate_files/main.nf'
include { REMOVE_TRAILING_FASTA_NS					 		 } from '../modules/local/remove_trailing_fasta_ns/main.nf'
include { SED as FINAL_CONSENSUS_SEQUENCE					 } from '../modules/local/sed/main.nf'
include { MINIMAP2_ALIGN_TO_FINAL		  					 } from '../modules/nf_core/minimap2/align/main.nf'

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
    
  // run PycoQC
  PYCOQC ( ch_summary )
    
  // align input reads using minimap2
  MINIMAP2_ALIGN_TO_EXISTING ( ch_reads, ch_reference )
  
  // extract new fasta file containing best aligned-to seqs for this dataset
  IDENTIFY_BEST_SEGMENTS_FROM_SAM ( MINIMAP2_ALIGN_TO_EXISTING.out.sam, ch_reference )
  
  // re-align data against best 10 BTV ref seqs using minimap2.
  MINIMAP2_ALIGN_TO_BEST10 ( ch_reads.join(IDENTIFY_BEST_SEGMENTS_FROM_SAM.out.fa))
  
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
  
  // rename file headers with unique id and segment number
  // RENAME_ONE_FASTA ( CALL_INDIVIDUAL_CONSENSUS_NANOPORE.out.viral_consensus_fasta.join(ch_best10_individual_fasta) )
  RENAME_ONE_FASTA ( CALL_INDIVIDUAL_CONSENSUS_NANOPORE.out.ivar_fasta.join(ch_best10_individual_fasta) )

  // collect individual consensus sequences and combine into single files  
  collected_vc_fasta_ch   = RENAME_ONE_FASTA.out.fasta.groupTuple()
  CONCATENATE_VC_FILES  (collected_vc_fasta_ch,   ".viral_consensus.fasta")
  
  collected_ivar_fasta_ch = RENAME_ONE_FASTA.out.fasta.groupTuple()
  CONCATENATE_IVAR_FILES(collected_ivar_fasta_ch, ".ivar_consensus.fasta")
  
  // pipe output through remove_trailing_fasta_Ns to strip N characters from beginning and ends of seqs
  // REMOVE_TRAILING_FASTA_NS ( CONCATENATE_VC_FILES.out.file )
  REMOVE_TRAILING_FASTA_NS ( CONCATENATE_IVAR_FILES.out.file )
  
  // pipe output through a sed to append new_X_draft_sequence to name of fasta record
  FINAL_CONSENSUS_SEQUENCE ( REMOVE_TRAILING_FASTA_NS.out.fa )
  
  // re-align data against the new draft sequence (ie. final consensus sequence) using minimap2.
  MINIMAP2_ALIGN_TO_FINAL ( ch_reads.join(FINAL_CONSENSUS_SEQUENCE.out.fa))
  
  }
  
  // specify the entry point for the workflow
workflow {
  NANOPORE_CONSENSUS()
}
