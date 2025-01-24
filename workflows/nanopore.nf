include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_TO_EXISTING  	 } from '../modules/local/minimap2/align/main.nf'
include { SAMTOOLS_VIEW	 as SAMTOOLS_VIEW_ALIGNMENT    		 } from '../modules/local/samtools/view/main.nf'
include { SAMTOOLS_SORT  as SAMTOOLS_SORT_ALIGNMENT    		 } from '../modules/local/samtools/sort/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGNMENT  		 } from '../modules/local/samtools/index/main.nf'
include { IDENTIFY_BEST_SEGMENTS_FROM_SAM     				 } from '../modules/local/identify_best_segments_from_sam/main.nf'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_TO_NEW_DRAFT  	 } from '../modules/local/minimap2/align/main2.nf'
include { SAMTOOLS_VIEW	 as SAMTOOLS_VIEW_BEST10_ALIGNMENT   } from '../modules/local/samtools/view/main.nf'
include { SAMTOOLS_SORT  as SAMTOOLS_SORT_BEST10_ALIGNMENT   } from '../modules/local/samtools/sort/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BEST10_ALIGNMENT  } from '../modules/local/samtools/index/main.nf'
include { SAMTOOLS_FAIDX  							     	 } from '../modules/local/samtools/faidx/main.nf'
include { BCFTOOLS_MPILEUP	 					   			 } from '../modules/local/bcftools/mpileup/main.nf'
include { CREATE_MASK_FILE				       			 	 } from '../modules/local/create_mask_file/main.nf'
include { BCFTOOLS_VIEW	 					   			 	 } from '../modules/local/bcftools/view/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_BCF	 			 } from '../modules/local/bcftools/index/main.nf'
include { BCFTOOLS_CALL 	 				   			 	 } from '../modules/local/bcftools/call/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_CONS	 			 } from '../modules/local/bcftools/index/main.nf'
include { BCFTOOLS_CONSENSUS					 			 } from '../modules/local/bcftools/consensus/main.nf'
include { REMOVE_TRAILING_FASTA_NS					 		 } from '../modules/local/remove_trailing_fasta_ns/main.nf'
include { SED as FINAL_CONSENSUS_SEQUENCE					 } from '../modules/local/sed/main.nf'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_TO_FINAL		  	 } from '../modules/local/minimap2/align/main3.nf'

workflow NANOPORE_CONSENSUS {

  ch_versions = Channel.empty()                                               

  // fastq input files

  Channel.fromFilePairs("${params.fastq_dir}/${params.input_pattern}", size: -1, checkIfExists: true, maxDepth: 1)
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

    Channel.fromPath("${params.reference_fasta}")
    .collect()
    .map { reference ->
            def meta2 = [:]
            meta2.id = "reference"
            [meta2, reference]
        }
    .set { ch_reference }
    
  // run minimap2 on input reads
  MINIMAP2_ALIGN_TO_EXISTING ( ch_reads, ch_reference )
  
  // run samtools to process minimap2 alignment - view
  SAMTOOLS_VIEW_ALIGNMENT ( MINIMAP2_ALIGN_TO_EXISTING.out.sam )
  
  // run samtools to process minimap2 alignment - sort
  SAMTOOLS_SORT_ALIGNMENT ( SAMTOOLS_VIEW_ALIGNMENT.out.bam )
  
  // run samtools to process minimap2 alignment - index
  SAMTOOLS_INDEX_ALIGNMENT ( SAMTOOLS_SORT_ALIGNMENT.out.bam )
  
  // extract new fasta file containing best aligned-to seqs for this dataset
  IDENTIFY_BEST_SEGMENTS_FROM_SAM ( MINIMAP2_ALIGN_TO_EXISTING.out.sam, ch_reference )
  
  // re-minimap data against best 10 BTV ref seqs.
  MINIMAP2_ALIGN_TO_NEW_DRAFT ( ch_reads.join(IDENTIFY_BEST_SEGMENTS_FROM_SAM.out.fa))
  
  // run samtools to process minimap2 alignment again - view
  SAMTOOLS_VIEW_BEST10_ALIGNMENT ( MINIMAP2_ALIGN_TO_NEW_DRAFT.out.sam )
  
  // run samtools to process minimap2 alignment again - sort
  SAMTOOLS_SORT_BEST10_ALIGNMENT ( SAMTOOLS_VIEW_BEST10_ALIGNMENT.out.bam )
  
   // run samtools to process minimap2 alignment - index
  SAMTOOLS_INDEX_BEST10_ALIGNMENT ( SAMTOOLS_SORT_BEST10_ALIGNMENT.out.bam )
  
  
  // commands below create "new" consensus sequences
  
  
  // have to make a .fai file to make mpileup happy
  SAMTOOLS_FAIDX ( IDENTIFY_BEST_SEGMENTS_FROM_SAM.out.fa )
  
  // lofreq calls variants -> output is a vcf file
  BCFTOOLS_MPILEUP ( SAMTOOLS_SORT_BEST10_ALIGNMENT.out.bam.join(IDENTIFY_BEST_SEGMENTS_FROM_SAM.out.fa))
  
  // this script creates a mask file
  // this is necessary because otherwise bcftools consensus doesn't hanlde positions with no coverage well
  CREATE_MASK_FILE ( BCFTOOLS_MPILEUP.out.vcf )
  
  // convert vcf -> compressed vcf to make bcftools happy
  BCFTOOLS_VIEW ( BCFTOOLS_MPILEUP.out.vcf ) 
  
  // need to make an indexed vcf file to make other bcftools commands happy
  BCFTOOLS_INDEX_BCF ( BCFTOOLS_VIEW.out.gz )
  
  // bcftools call creates a new vcf file that has consensus bases 
  BCFTOOLS_CALL ( BCFTOOLS_INDEX_BCF.out.gz_and_csi )
  
  // need to make an indexed cons.vcf.gz file to make other bcftools commands happy
  BCFTOOLS_INDEX_CONS ( BCFTOOLS_CALL.out.vcf )
  
  // bcftools consensus will output a fasta file containing new draft consensus sequence based on called variants
  BCFTOOLS_CONSENSUS ( CREATE_MASK_FILE.out.mask.join(IDENTIFY_BEST_SEGMENTS_FROM_SAM.out.fa).join(BCFTOOLS_INDEX_CONS.out.gz_and_csi)) 
  
  // pipe output through remove_trailing_fasta_Ns to strip N characters from beginning and ends of seqs
  REMOVE_TRAILING_FASTA_NS ( BCFTOOLS_CONSENSUS.out.fa )
  
  // pipe output through a sed to append new_X_draft_sequence to name of fasta record
  FINAL_CONSENSUS_SEQUENCE ( REMOVE_TRAILING_FASTA_NS.out.fa )
  
  // re-minimap data against the new draft sequence (ie. final consensus sequence)
  MINIMAP2_ALIGN_TO_FINAL ( ch_reads.join(BCFTOOLS_CONSENSUS.out.fa))
  
  }
  
  // specify the entry point for the workflow
workflow {
  NANOPORE_CONSENSUS()
}
