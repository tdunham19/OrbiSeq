include { MINIMAP2_ALIGN     } from '../modules/nf_core/minimap2/align/main.nf'
include { VIRAL_CONSENSUS    } from '../modules/local/viral_consensus/main.nf'
// include { IVAR_CONSENSUS     } from '../modules/nf_core/ivar/consensus/main.nf'

workflow CALL_INDIVIDUAL_CONSENSUS_NANOPORE {

 take:
  reads_refseq       // tuple val (meta), path(reads), path(refseq)
  min_depth
  min_qual
  min_freq
  
 main:

  // define some empty channels for keeping track of stuff
  ch_versions         = Channel.empty()                                               

  // map reads to refseq
  MINIMAP2_ALIGN(reads_refseq, "individual_refseq")
  ch_versions = ch_versions.mix ( MINIMAP2_ALIGN.out.versions )      

  // call consensus using viral_consensus
  VIRAL_CONSENSUS(MINIMAP2_ALIGN.out.bam_refseq, min_qual, min_depth, min_freq)
  ch_versions = ch_versions.mix ( VIRAL_CONSENSUS.out.versions )
  
  // call consensus using ivar
  // def save_mpileup = true
  // IVAR_CONSENSUS(MINIMAP2_ALIGN.out.bam_refseq, min_qual, min_depth, min_freq, save_mpileup)
  // ch_versions = ch_versions.mix ( IVAR_CONSENSUS.out.versions )      

 emit: 
  versions                        = ch_versions
  viral_consensus_fasta           = VIRAL_CONSENSUS.out.fasta
  viral_consensus_position_counts = VIRAL_CONSENSUS.out.position_counts
  viral_consensus_refseq          = VIRAL_CONSENSUS.out.refseq
  // ivar_fasta                      = IVAR_CONSENSUS.out.fasta
  // ivar_qual                       = IVAR_CONSENSUS.out.qual
  // ivar_mpileup                    = IVAR_CONSENSUS.out.mpileup
}
