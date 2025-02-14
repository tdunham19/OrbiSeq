include { CHECK_FASTQ_COMPRESSED   } from '../../modules/stenglein_lab/check_fastq_compressed/main'
include { COUNT_FASTQ              } from '../../modules/stenglein_lab/count_fastq/main'


/*
 Create a channel with input fastq, possibly from multiple directories (specified as a comma-separated list)
 also returns IDs, determeind from filenames stripped of extensions and Illumina-added stuff like _001_ _S?_ etc

 Asumptions about input fastq:

 - gzip compressed
 - filenames end in .fastq.gz or .fq.gz
 - single-end or paired en OK
 - filenames matchable using param fastq_pattern
 - in one or more directories specified by param input_fastq_dir
 */
workflow MARSHAL_FASTQ {

 take:
 input_fastq_dir        // the path to a directory containing fastq file(s) or a comma-separated list of dirs
 fastq_pattern          // the regex that will be matched to identify fastq

 main:

  ch_versions         = Channel.empty()
  ch_fastq_counts     = Channel.empty()

  // User can specify multiple directories containing input fastq
  // In this case, the directories should be provided as a 
  // comma-separated list (no spaces)
  def fastq_dirs = input_fastq_dir.tokenize(',')

  // construct list of directories in which to find fastq
  fastq_dir_list = []
  for (dir in fastq_dirs){
     def file_pattern = "${dir}/${fastq_pattern}"
     fastq_dir_list.add(file_pattern)
  }

  /*
   These fastq files are the main input to this workflow
  */
  Channel
  .fromFilePairs(fastq_dir_list, size: -1, checkIfExists: true, maxDepth: 1)
  .map{ name, reads ->

         // define a new empty map named meta for each sample
         // and populate it with id and single_end values
         // for compatibility with nf-core module expected parameters
         // reads are just the list of fastq
         def meta        = [:]

         // this map gets rid of any of the following at the end of sample IDs:
         // .gz
         // .fastq
         // .fq
         // _001
         // _R[12]
         // _S\d+ 
         // E.g. strip _S1 from the end of a sample ID..  
         // This is typically sample #s from Illumina basecalling.
         // could cause an issue if sequenced the same sample with 
         // multiple barcodes so was repeated on a sample sheet. 
         meta.id         = name.replaceAll( /.gz$/ ,"")
         meta.id         = meta.id.replaceAll( /.fastq$/ ,"")
         meta.id         = meta.id.replaceAll( /.fq$/ ,"")
         meta.id         = meta.id.replaceAll( /.uniq$/ ,"")
         meta.id         = meta.id.replaceAll( /.trim$/ ,"")
         meta.id         = meta.id.replaceFirst( /_001$/ ,"")
         meta.id         = meta.id.replaceFirst( /_R[12]$/ ,"")
         meta.id         = meta.id.replaceFirst( /_S\d+$/ ,"")

         // if 2 fastq files then paired end data, so single_end is false
         meta.single_end = reads[1] ? false : true

         // this last statement in the map closure is the return value
         [ meta, reads ] }

  .set { ch_reads }

  // double check input is compressed - will error if not
  CHECK_FASTQ_COMPRESSED ( ch_reads ) 

  // count # of reads in each fastq file (or file pair)
  COUNT_FASTQ ( ch_reads.map{ meta, reads -> [ meta, reads, "initial"] } )
  ch_fastq_counts = ch_fastq_counts.mix(COUNT_FASTQ.out.count_file)
  
 emit: 
  reads           = ch_reads
  versions        = ch_versions
  fastq_counts    = ch_fastq_counts

}
