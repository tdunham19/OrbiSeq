include { MARSHAL_FASTQ                    }         from '../../subworkflows/stenglein_lab/marshal_fastq'

include { FASTQC   as FASTQC_PRE           }         from '../../modules/nf_core/fastqc/main'
include { FASTQC   as FASTQC_POST_TRIM     }         from '../../modules/nf_core/fastqc/main'
include { FASTQC   as FASTQC_POST_COLLAPSE }         from '../../modules/nf_core/fastqc/main'

include { MULTIQC  as MULTIQC_PRE     }              from '../../modules/nf_core/multiqc/main' 
include { MULTIQC  as MULTIQC_POST    }              from '../../modules/nf_core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS }              from '../../modules/nf_core/custom/dumpsoftwareversions/main'

include { CHECK_FASTQ_COMPRESSED      }              from '../../modules/stenglein_lab/check_fastq_compressed/main'
include { CUTADAPT                    }              from '../../modules/stenglein_lab/cutadapt/main'
include { BBMAP_BBDUK                 }              from '../../modules/nf_core/bbmap/bbduk/main'
include { COLLAPSE_DUPLICATE_READS    }              from '../../modules/stenglein_lab/collapse_duplicate_reads/main'

include { COUNT_FASTQ as COUNT_FASTQ_INITIAL       } from '../../modules/stenglein_lab/count_fastq/main'
include { COUNT_FASTQ as COUNT_FASTQ_POST_TRIM     } from '../../modules/stenglein_lab/count_fastq/main'
include { COUNT_FASTQ as COUNT_FASTQ_POST_COLLAPSE } from '../../modules/stenglein_lab/count_fastq/main'

include { SAVE_OUTPUT_FILE as SAVE_FASTQ_OUTPUT       } from '../../modules/stenglein_lab/save_output_file/main'
include { SAVE_OUTPUT_FILE as SAVE_FASTQ_DEDUP_OUTPUT } from '../../modules/stenglein_lab/save_output_file/main'

include { PROCESS_FASTQ_COUNTS        }              from '../../modules/stenglein_lab/process_fastq_counts/main'

workflow PREPROCESS_READS {

 take:
  input_fastq_dir        // the path to a directory containing fastq file(s) or a comma-separated list of dirs
  fastq_pattern          // the regex that will be matched to identify fastq
  collapse_duplicates    // boolean: collapse duplicate reads?

 main:

  // define some empty channels for keeping track of stuff
  ch_versions         = Channel.empty()                                               
  ch_fastq_counts     = Channel.empty()                                               
  ch_processed_fastq  = Channel.empty()                                               

  // make sure all the fastq are in order
  MARSHAL_FASTQ(input_fastq_dir, fastq_pattern)
  ch_reads        = MARSHAL_FASTQ.out.reads
  ch_fastq_counts = ch_fastq_counts.mix(MARSHAL_FASTQ.out.fastq_counts)
  ch_versions     = ch_versions.mix ( MARSHAL_FASTQ.out.versions )      
  
  // run fastqc on input reads
  FASTQC_PRE ( ch_reads )
  ch_versions = ch_versions.mix ( FASTQC_PRE.out.versions )      

  // trim low quality bases and adapters from reads
  CUTADAPT ( ch_reads )
  ch_versions = ch_versions.mix ( CUTADAPT.out.versions ) 

  // BBDUK to remove adapter-containing reads that made it past cutadapt adapter trimming
  // this can happen, for instance, if the adapter sequence is truncated on its 5' end 
  // (perhaps resulting from ligation of a partially-degraded adapter molecule?)
  BBMAP_BBDUK(CUTADAPT.out.reads, params.bbduk_adapters)
  ch_versions = ch_versions.mix ( BBMAP_BBDUK.out.versions ) 

  // count numbers of reads after trimming
  COUNT_FASTQ_POST_TRIM ( BBMAP_BBDUK.out.reads.map{ meta, reads -> [ meta, reads, "post_trimming"] } )
  ch_fastq_counts = ch_fastq_counts.mix(COUNT_FASTQ_POST_TRIM.out.count_file)

  // run fastqc on post trimmed reads
  FASTQC_POST_TRIM ( BBMAP_BBDUK.out.reads )

  ch_processed_reads = BBMAP_BBDUK.out.reads
  SAVE_FASTQ_OUTPUT(ch_processed_reads.map{meta, reads -> reads}.flatten())

  // optionally collapse duplicates reads 
  if (collapse_duplicates) {
    COLLAPSE_DUPLICATE_READS ( BBMAP_BBDUK.out.reads.filter{ it[1]*.getAt(0).size() > 0} ) 
    ch_versions = ch_versions.mix ( COLLAPSE_DUPLICATE_READS.out.versions ) 
  
    // run fastqc on post collapsed reads
    FASTQC_POST_COLLAPSE ( COLLAPSE_DUPLICATE_READS.out.reads )
  
    // count reads after collapsing
    COUNT_FASTQ_POST_COLLAPSE ( COLLAPSE_DUPLICATE_READS.out.reads.map{ meta, reads -> [ meta, reads, "post_collapse"] } )
    ch_fastq_counts = ch_fastq_counts.mix(COUNT_FASTQ_POST_COLLAPSE.out.count_file)

    ch_processed_reads = COLLAPSE_DUPLICATE_READS.out.reads
    SAVE_FASTQ_DEDUP_OUTPUT(ch_processed_reads.map{meta, reads -> reads}.flatten())
  }


  // MultiQC                                                          
  // workflow_summary    = WorkflowRnaseq.paramsSummaryMultiqc(workflow, summary_params)
  // ch_workflow_summary = Channel.value(workflow_summary)                   

  // collect files that will be input to multiqc
  ch_multiqc_pre_files = Channel.empty()
  ch_multiqc_pre_files = ch_multiqc_pre_files.mix(FASTQC_PRE.out.zip.collect{it[1]}.ifEmpty([]))
  ch_multiqc_pre_files = ch_multiqc_pre_files.mix(CUTADAPT.out.log.collect{it[1]}.ifEmpty([]))
  ch_multiqc_pre_files = ch_multiqc_pre_files.mix(BBMAP_BBDUK.out.log.collect{it[1]}.ifEmpty([]))
  if (collapse_duplicates) {
    ch_multiqc_pre_files = ch_multiqc_pre_files.mix(FASTQC_POST_COLLAPSE.out.zip.collect{it[1]}.ifEmpty([]))
  }
  // run multiqc on pre-trimmed files and trimming output
  MULTIQC_PRE (
      ch_multiqc_pre_files.collect(), [], [], []
  )
  ch_versions    = ch_versions.mix(MULTIQC_PRE.out.versions)

  // run multiqc on pre-trimmed files and trimming output
  ch_multiqc_post_files = Channel.empty()
  ch_multiqc_post_files = ch_multiqc_post_files.mix(FASTQC_POST_TRIM.out.zip.collect{it[1]}.ifEmpty([]))
  MULTIQC_POST (
      ch_multiqc_post_files.collect(), [], [], []
  )
  ch_versions    = ch_versions.mix(MULTIQC_POST.out.versions)

  // combine all multiqc input files
  ch_multiqc_files = Channel.empty()
  ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_pre_files)
  ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_post_files)

  // ch_multiqc_pre_files = ch_multiqc_pre_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

  /*
     Process all fastq counts and create PDF output
   */
  PROCESS_FASTQ_COUNTS(ch_fastq_counts.collectFile(name: "all_fastq_counts.txt"))


 emit: 
  versions        = ch_versions
  multiqc_files   = ch_multiqc_files 
  fastq_counts    = ch_fastq_counts 
  reads           = ch_processed_reads

}
