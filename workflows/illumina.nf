#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include Stenglein lab read_preprocessing pipeline. 
include { PREPROCESS_READS 									 } from '../subworkflows/stenglein_lab/preprocess_reads.nf'

// include other modules: local, nf-core, and Stenglein lab
include { BOWTIE2_BUILD as BOWTIE2_BUILD_INDEX_EXISTING      } from '../modules/nf_core/bowtie2/build/main.nf'
include { BOWTIE2_ALIGN_TO_EXISTING  		 				 } from '../modules/nf_core/bowtie2/align/main.nf'
include { IDENTIFY_BEST_SEGMENTS_FROM_SAM     				 } from '../modules/local/identify_best_segments_from_sam/main.nf'

include { CALL_INDIVIDUAL_CONSENSUS_ILLUMINA              	 } from '../subworkflows/call_individual_consensus_illumina.nf'

include { RENAME_ONE_FASTA 									 } from '../modules/local/rename_one_fasta/main.nf'
include { CONCATENATE_FILES as CONCATENATE_VC_FILES          } from '../modules/stenglein_lab/concatenate_files/main.nf'
include { CONCATENATE_FILES as CONCATENATE_IVAR_FILES        } from '../modules/stenglein_lab/concatenate_files/main.nf'
include { REMOVE_TRAILING_FASTA_NS					 		 } from '../modules/local/remove_trailing_fasta_ns/main.nf'
include { SED as FINAL_CONSENSUS_SEQUENCE					 } from '../modules/local/sed/main.nf'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_INDEX_FINAL	     } from '../modules/nf_core/bowtie2/build/main.nf'
include { BOWTIE2_ALIGN_TO_FINAL           				     } from '../modules/nf_core/bowtie2/align/main.nf'

workflow ILLUMINA_CONSENSUS {

    ch_versions = Channel.empty()                                              

  // refseq input files

    Channel.fromPath("${params.reference}")
    .collect()
    .map { reference ->
            def meta2 = [:]
            meta2.id = "reference"
            [meta2, reference]
        }
    .set { ch_reference } 
  
   // use Stenglein lab read_preprocessing pipeline 
  PREPROCESS_READS(params.fastq_dir, params.input_pattern, params.collapse_duplicate_reads)
  
  // build bowtie2 index
  BOWTIE2_BUILD_INDEX_EXISTING ( ch_reference )
    
  // run bowtie2-align on input reads with large reference  
  ch_processed_reads = PREPROCESS_READS.out.reads
  BOWTIE2_ALIGN_TO_EXISTING (ch_processed_reads, BOWTIE2_BUILD_INDEX_EXISTING.out.index )
  
  // extract new fasta file containing best aligned-to seqs for this dataset
  IDENTIFY_BEST_SEGMENTS_FROM_SAM ( BOWTIE2_ALIGN_TO_EXISTING.out.sam, ch_reference )
  
  // split up best10 segments into individual sequences because virus-focused 
  // consensus callers (namely iVar and viral_consensus) only work on one
  // sequence at a time
  // see: https://www.nextflow.io/docs/latest/reference/operator.html#splitfasta
  IDENTIFY_BEST_SEGMENTS_FROM_SAM.out.fa
    .splitFasta(by: 1, file: true, elem: 1)
    .set { ch_best10_individual_fasta }
    
  // this uses the nextflow combine operator to create a new channel 
  // that contains the reads for each dataset and all individual fasta files 
  individual_fasta_ch = ch_processed_reads.combine(ch_best10_individual_fasta, by: 0)
  
  // parameters related consensus calling: min depth, basecall quality, frequency for consensus calling
  min_depth_ch = Channel.value(params.illumina_min_depth)
  min_qual_ch  = Channel.value(params.illumina_min_qual)
  min_freq_ch  = Channel.value(params.illumina_min_freq)
  CALL_INDIVIDUAL_CONSENSUS_ILLUMINA(individual_fasta_ch, min_depth_ch, min_qual_ch, min_freq_ch)
  
  // rename file headers with unique id and segment number
  // RENAME_ONE_FASTA ( CALL_INDIVIDUAL_CONSENSUS_ILLUMINA.out.viral_consensus_refseq.join(CALL_INDIVIDUAL_CONSENSUS_ILLUMINA.out.viral_consensus_fasta) )
  RENAME_ONE_FASTA ( CALL_INDIVIDUAL_CONSENSUS_ILLUMINA.out.viral_consensus_refseq.join(CALL_INDIVIDUAL_CONSENSUS_ILLUMINA.out.ivar_fasta) )

  // collect individual consensus sequences and combine into single files
  // collected_vc_fasta_ch   = RENAME_ONE_FASTA.out.fasta.groupTuple()
  // CONCATENATE_VC_FILES  (collected_vc_fasta_ch,   ".viral_consensus.fasta")

  collected_ivar_fasta_ch = RENAME_ONE_FASTA.out.fasta.groupTuple()
  CONCATENATE_IVAR_FILES(collected_ivar_fasta_ch, ".ivar_consensus.fasta")
  
  // pipe output through remove_trailing_fasta_Ns to strip N characters from beginning and ends of seqs
  // REMOVE_TRAILING_FASTA_NS ( CONCATENATE_VC_FILES.out.file )
  REMOVE_TRAILING_FASTA_NS ( CONCATENATE_IVAR_FILES.out.file )
  
  // pipe output through a sed to append new_X_draft_sequence to name of fasta record
  FINAL_CONSENSUS_SEQUENCE ( REMOVE_TRAILING_FASTA_NS.out.fa )
  
   // filter out empty fasta files
  FINAL_CONSENSUS_SEQUENCE_FILTERED = FINAL_CONSENSUS_SEQUENCE.out.filter { meta, fasta ->
    fasta.text.readLines().find { it && !it.startsWith(">") } != null 
    }
  
  // make bowtie2 index for final alignment
  BOWTIE2_BUILD_INDEX_FINAL ( FINAL_CONSENSUS_SEQUENCE_FILTERED )
   
  // re-align data against the new draft sequence (ie. final consensus sequence) using bowtie2. 
  BOWTIE2_ALIGN_TO_FINAL ( ch_processed_reads.join(BOWTIE2_BUILD_INDEX_FINAL.out.index) )
  
  }
  
  // specify the entry point for the workflow
workflow {
  ILLUMINA_CONSENSUS()
}
