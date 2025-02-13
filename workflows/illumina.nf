// clone and include read preprocessing subworkflow
// include { CLONE_SUBWORKFLOW_REPO 							 } from '../modules/local/clonerepo/main.nf'
// include { READ_PREPROCESSING     							 } from '../subworkflows/main.nf'

// include other nf-core and Stenglein lab/local modules
include { BOWTIE2_BUILD as BOWTIE2_BUILD_INDEX_EXISTING      } from '../modules/nf_core/bowtie2/build/main.nf'
include { BOWTIE2_ALIGN_TO_EXISTING  		 				 } from '../modules/nf_core/bowtie2/align/main.nf'
include { SAMTOOLS_VIEW	 as SAMTOOLS_VIEW_ALIGNMENT    		 } from '../modules/nf_core/samtools/view/main.nf'
include { SAMTOOLS_SORT  as SAMTOOLS_SORT_ALIGNMENT    		 } from '../modules/nf_core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGNMENT  		 } from '../modules/nf_core/samtools/index/main.nf'
include { IDENTIFY_BEST_SEGMENTS_FROM_SAM     				 } from '../modules/local/identify_best_segments_from_sam/main.nf'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_INDEX_BEST10        } from '../modules/nf_core/bowtie2/build/main2.nf'
include { BOWTIE2_ALIGN_TO_NEW_DRAFT   	 				     } from '../modules/nf_core/bowtie2/align/main.nf'
include { SAMTOOLS_VIEW	 as SAMTOOLS_VIEW_BEST10_ALIGNMENT   } from '../modules/nf_core/samtools/view/main.nf'
include { SAMTOOLS_SORT  as SAMTOOLS_SORT_BEST10_ALIGNMENT   } from '../modules/nf_core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BEST10_ALIGNMENT  } from '../modules/nf_core/samtools/index/main.nf'
include { SAMTOOLS_FAIDX  							     	 } from '../modules/nf_core/samtools/faidx/main.nf'
include { BCFTOOLS_MPILEUP	 					   			 } from '../modules/nf_core/bcftools/mpileup/main.nf'
include { CREATE_MASK_FILE				       			 	 } from '../modules/local/create_mask_file/main.nf'
include { BCFTOOLS_VIEW	 					   			 	 } from '../modules/nf_core/bcftools/view/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_BCF	 			 } from '../modules/nf_core/bcftools/index/main.nf'
include { BCFTOOLS_CALL 	 				   			 	 } from '../modules/nf_core/bcftools/call/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_CONS	 			 } from '../modules/nf_core/bcftools/index/main.nf'
include { BCFTOOLS_CONSENSUS					 			 } from '../modules/nf_core/bcftools/consensus/main.nf'
include { REMOVE_TRAILING_FASTA_NS					 		 } from '../modules/local/remove_trailing_fasta_ns/main.nf'
include { SED as FINAL_CONSENSUS_SEQUENCE					 } from '../modules/local/sed/main.nf'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_INDEX_FINAL	     } from '../modules/nf_core/bowtie2/build/main3.nf'
include { BOWTIE2_ALIGN_TO_FINAL           				     } from '../modules/nf_core/bowtie2/align/main.nf'
// include { MULTIQC 											 } from './modules/nf_core/multiqc/main.nf'

workflow ILLUMINA_CONSENSUS {

    ch_versions = Channel.empty()                                              

  // fastq input files  
  
Channel.fromFilePairs("${params.fastq_dir}/${params.input_pattern}", size: -1, checkIfExists: true, maxDepth: 1)
  .map{ name, reads ->
         def meta        = [:]
         meta.id         = name.replaceAll( /.gz$/ ,"")
         meta.id         = meta.id.replaceAll( /.fastq$/ ,"")
         meta.id         = meta.id.replaceAll( /.fq$/ ,"")
         meta.id         = meta.id.replaceAll( /.uniq$/ ,"")
         meta.id         = meta.id.replaceAll( /.trim$/ ,"")
         meta.id         = meta.id.replaceFirst( /_001$/ ,"")
         meta.id         = meta.id.replaceFirst( /_R[12]$/ ,"")
         meta.id         = meta.id.replaceFirst( /_S\d+$/ ,"")
         meta.single_end = reads[1] ? false : true
         [ meta, reads ] }

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
    
  // clone Stenglein lab read preprocessing pipeline
  // CLONE_SUBWORKFLOW_REPO ()  
      
  // use Stenglein lab read preprocessing pipeline to do fastqc and trim reads. 
  // READ_PREPROCESSING ( ch_reads, collapse_duplicate_reads: params.collapse_duplicate_reads, fastq_dir: params.fastq_dir )

  // Build bowtie2 index
  BOWTIE2_BUILD_INDEX_EXISTING ( ch_reference )
    
  // run bowtie2-align on input reads with large reference
  BOWTIE2_ALIGN_TO_EXISTING ( ch_reads, BOWTIE2_BUILD_INDEX_EXISTING.out.index )
 
// Trying to get the suffix added in for one main.nf
  // BOWTIE2_ALIGN_TO_EXISTING ( READ_PREPROCESSING.out.gz, BOWTIE2_BUILD_INDEX_EXISTING.out.index )
  // BOWTIE2_ALIGN_TO_EXISTING ( ch_reads, BOWTIE2_BUILD_INDEX_EXISTING.out.index, "initial" )
  // BOWTIE2_ALIGN_TO_EXISTING ( ch_reads.join(BOWTIE2_BUILD_INDEX_EXISTING.out.index).map { meta, reads, index, file_suffix -> [meta, reads, index, "initial"] } )
  // ch_align_to_existing = ch_reads.join(BOWTIE2_BUILD_INDEX_EXISTING.out.index).map {meta, reads, index -> [meta, reads, index, "initial"] }
  // BOWTIE2_ALIGN_TO_EXISTING( ch_align_to_existing )
  // def ch_initial = Channel.from(ch_reads, BOWTIE2_BUILD_INDEX_EXISTING.out.index)
  // Channel.from( [ch_reads], [BOWTIE2_BUILD_INDEX_EXISTING.out.index] ).collect().set { ch_initial_reads_index }
  // ch_initial = ( ch_reads.combine(BOWTIE2_BUILD_INDEX_EXISTING.out.index) )
  // BOWTIE2_ALIGN_TO_EXISTING( ch_initial.map {meta, initial, file_suffix -> [meta, initial, "initial"] } )
  
  // run samtools to process bowtie2 alignment - view converts .sam to .bam 
  SAMTOOLS_VIEW_ALIGNMENT ( BOWTIE2_ALIGN_TO_EXISTING.out.sam )
  
  // run samtools to process bowtie2 alignment - sort will sort the .bam alignment
  SAMTOOLS_SORT_ALIGNMENT ( SAMTOOLS_VIEW_ALIGNMENT.out.bam )
  
  // run samtools to process bowtie2 alignment - index allows for fast random access for alignment files
  SAMTOOLS_INDEX_ALIGNMENT ( SAMTOOLS_SORT_ALIGNMENT.out.bam )
  
  // extract new fasta file containing best aligned-to seqs for this dataset
  IDENTIFY_BEST_SEGMENTS_FROM_SAM ( BOWTIE2_ALIGN_TO_EXISTING.out.sam, ch_reference )
  
  // build bowtie index for best 10 BTV ref seqs.
  BOWTIE2_BUILD_INDEX_BEST10 ( IDENTIFY_BEST_SEGMENTS_FROM_SAM.out.fa )
  
  // re-align data against best 10 BTV ref seqs.
  BOWTIE2_ALIGN_TO_NEW_DRAFT ( ch_reads.join(BOWTIE2_BUILD_INDEX_BEST10.out.index) )
  // BOWTIE2_ALIGN_TO_NEW_DRAFT ( ch_reads.join(BOWTIE2_BUILD_INDEX_BEST10.out.index).map {meta, reads, index -> tuple(meta, reads, index, "best10") } )

  ch_best10=ch_reads.join(BOWTIE2_BUILD_INDEX_BEST10.out.index)

  // run samtools to process bowtie2 alignment again - view converts .sam to .bam 
  SAMTOOLS_VIEW_BEST10_ALIGNMENT ( BOWTIE2_ALIGN_TO_NEW_DRAFT.out.sam )
  
  // run samtools to process bowtie2 alignment again - sort will sort the .bam alignment
  SAMTOOLS_SORT_BEST10_ALIGNMENT ( SAMTOOLS_VIEW_BEST10_ALIGNMENT.out.bam )
  
   // run samtools to process bowtie2 alignment - index allows for fast random access for alignment files
  SAMTOOLS_INDEX_BEST10_ALIGNMENT ( SAMTOOLS_SORT_BEST10_ALIGNMENT.out.bam )
  
  
  // commands below create "new" consensus sequences
  
  
  // have to make a .fai file to make mpileup happy - faidx allows for fast random access for reference files in fasta format 
  SAMTOOLS_FAIDX ( IDENTIFY_BEST_SEGMENTS_FROM_SAM.out.fa )
  
  // bcftools mpileup calls variants -> output is a vcf file
  BCFTOOLS_MPILEUP ( SAMTOOLS_SORT_BEST10_ALIGNMENT.out.bam.join(IDENTIFY_BEST_SEGMENTS_FROM_SAM.out.fa))

  // this script creates a mask file which is necessary because otherwise bcftools consensus doesn't hanlde positions with no coverage well
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
  
  // make bowtie2 index for final alignment
  BOWTIE2_BUILD_INDEX_FINAL ( BCFTOOLS_CONSENSUS.out.fa )
   
  // re-align data against the new draft sequence (ie. final consensus sequence)
  BOWTIE2_ALIGN_TO_FINAL ( ch_reads.join(BOWTIE2_BUILD_INDEX_FINAL.out.index) )
  // BOWTIE2_ALIGN_TO_FINAL ( ch_reads.join(BOWTIE2_BUILD_INDEX_FINAL.out.index).map {meta, reads, index -> tuple(meta, reads, index, "new_draft_seq") } )

  // run multiqc on final output
  // MULTIQC
  
  }
  
  // specify the entry point for the workflow
workflow {
  ILLUMINA_CONSENSUS()
}
