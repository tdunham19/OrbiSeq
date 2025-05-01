include { BOWTIE2_BUILD_ALIGN } from '../modules/nf_core/bowtie2/build_align/main.nf'
include { VIRAL_CONSENSUS    } from '../modules/local/viral_consensus/main.nf'
include { IVAR_CONSENSUS     } from '../modules/nf_core/ivar/consensus/main.nf'

workflow CALL_INDIVIDUAL_CONSENSUS_ILLUMINA {

 take:
  reads_refseq    // tuple val (meta), path(reads), path(fasta)
  min_depth
  min_qual
  min_freq

 main:

  // define some empty channels for keeping track of stuff
  ch_versions         = Channel.empty()                                               

  // define some inputs to BT2 align
  suffix         = Channel.value("individual_refseq")
  save_unaligned = Channel.value(false)
  sort_bam       = Channel.value(true)

  // build bt2 index and align reads
  BOWTIE2_BUILD_ALIGN (reads_refseq, suffix, save_unaligned, sort_bam)
  ch_versions = ch_versions.mix ( BOWTIE2_BUILD_ALIGN.out.versions )

  // call consensus using viral_consensus
  VIRAL_CONSENSUS(BOWTIE2_BUILD_ALIGN.out.bam_fasta, min_qual, min_depth, min_freq)
  ch_versions = ch_versions.mix ( VIRAL_CONSENSUS.out.versions )

  // call consensus using ivar
  // def save_mpileup = true
  // IVAR_CONSENSUS(BOWTIE2_BUILD_ALIGN.out.bam_fasta, min_qual, min_depth, min_freq, save_mpileup)
  // ch_versions = ch_versions.mix ( IVAR_CONSENSUS.out.versions )      



 emit: 
  versions                        = ch_versions
  viral_consensus_refseq_and_new  = VIRAL_CONSENSUS.out.refseq_and_new
  viral_consensus_fasta           = VIRAL_CONSENSUS.out.fasta
  viral_consensus_position_counts = VIRAL_CONSENSUS.out.position_counts
  viral_consensus_refseq          = VIRAL_CONSENSUS.out.refseq
  bowtei2_build_align_bam_ref   = BOWTIE2_BUILD_ALIGN.out.bam_fasta
  // ivar_refseq_and_new             = IVAR_CONSENSUS.out.refseq_and_new  
  // ivar_fasta                      = IVAR_CONSENSUS.out.fasta
  // ivar_qual                       = IVAR_CONSENSUS.out.qual
  // ivar_mpileup                    = IVAR_CONSENSUS.out.mpileup
}