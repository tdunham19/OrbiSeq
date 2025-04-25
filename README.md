# OrbiSeq
OrbiSeq is a Nextflow pipeline to analyze short and long sequences for *Orbiviruses*.

## Contents
- [Pipeline Overview](#Pipeline-Overview)
- [Input Files](#Input-Files)
- [Output Files](#Output-Files)
- [Workflow Steps](#Workflow-Steps)
- [Running the Pipeline](#Running-the-Pipeline)
	- [Optional Deduplication](#Optional-Deduplication)
- [Testing](#Testing)
- [Stopping and Resuming](#Stopping-and-Resuming)
- [Dependencies](#Dependencies)


## Pipeline Overview
OrbiSeq is a nextflow pipeline that will perform alignment and create consensus sequences from Illumina or Nanopore sequence data for any *Orbivirus* with 10 segments. 

### Platform:
This pipeline has the capability to run either Illumina or Nanopore sequencing data. When running this pipeline the user must specify which platform their data came from by using the --platform parameter. 


## Input Files

### Data:
- Users will direct the pipeline where fastq files are located by using the --fastq_dir parameter. Illumina and Nanopore data should not be run through the pipeline at the same time. 
	- Illumina: 
		- Paired-end (R1 and R2) or single-end reads.
	- Nanopore: 
		- Users must concatenate raw read files into a single fastq.gz file (one file per barcode). This can be done using the cat command. 
		- The sequence summary file should be placed in the ./summary directory and begin with {sequencing_summary_}*.txt. 

### Reference:

This pipeline can utilize any reference genome from *Orbiviruses* with 10 segments. To include a reference the user must specify where the file is located using the parameter --reference. 
- Premade Reference Sequences 
	- There are premade reference files for Bluetongue virus (BTV) and Epizootic hemorrhagic disease virus (EHDV). More information on the creation of these references can be found in [cite publication]. 
		- The reference files for BTV and EHDV can be found in ./reference/BTV and ./reference/EHDV respectively. 
- User Included Reference Sequences	
	- The user is also able to upload their own custom reference file. 
	- When uploading custom reference sequences the user must ensure that the file is formatted correctly. It should be formatted as: segment#_sample_id (ex. s1_Genbank_Acession). Additionally there should be NO BLANK LINES throughout the document. 

## Output files

- The pipeline outputs a consensus sequence, alignment, and variant calling files. The files are located in the following results directories:
	
	- Illumina 
		- Quality Assessment: 
			- FastQC: {sample.id}_fastqc.html
			- MultiQC: *_multiqc_report.html in ./results/multiqc 
		- Consensus sequences: 
			- ivar - {sample_id}.ivar_consensus.fasta in ./results/concatenate
			- ViralConsensus - {sample_id}._new_draft_seq.fa in ./results/final
		- Final alignment: minimap2 - {sample_id}_new_draft_seq.sam file in ./results/bowtie2
			
	- Nanopore 
		- Quality Assessment: PycoQC - summary.html in ./results/pycoqc
		- Consensus sequence: ViralConsensus - {sample_id}.new_draft_seqs.fa in ./results/final
		- Final alignment: minimap2 - {sample_id}_new_draft_seq.sam file in ./results/minimap2

## Workflow Steps

### Illumina workflow 
- Preprocessing and quality assessment of input reads : Stenglein Lab read_preprocessing pipeline
	- Optional collapse duplicate reads 
- Align input reads to large Orbi RefSeq : bowtie2 build & align 
- Choose best 10 segments from initial alignment
- Align input reads to best10 RefSeq : bowtie2 build & align 
- Call variants & create consensus sequences: iVar & ViralConsensus
- Align input reads to final consensus sequence : bowtie2 build & align 

### Nanopore workflow 
- Quality assessment of input reads : PycoQC
- Align input reads to large Orbi RefSeq : minimap2 align 
- Choose best 10 segments from initial alignment
- Align input reads to best10 RefSeq : minimap2 align 
- Generate consensus sequences: ViralConsensus
- Align input reads to final consensus sequence : minimap2 align
	
These workflows take advantage of nf-core [modules](https://nf-co.re/modules) for many of these components and the overall [nf-core](https://nf-co.re/) design philosophy.

The illumina workflow takes advantage of the [Stenglein Lab Read Preprocessing Pipeline](https://github.com/stenglein-lab/read_preprocessing) and the iVar package[iVar](https://github.com/andersen-lab/ivar?tab=readme-ov-file). 

the nanopore workflow takes advantage of the ViralConsensus tool [ViralConsensus](https://github.com/niemasd/ViralConsensus). 


## Running the Pipeline

1. Clone the pipeline from github and move into the directory
```
git clone https://github.com/tdunham19/OrbiSeq.git
cd OrbiSeq
```

2. Run the pipeline: must specify sequencing platform and where the reference file is located.  
```
nextflow run main.nf --platform ['illumina' or 'nanopore'] --fastq_dir /path/to/fastq/directory --reference /path/to/reference/directory  -resume
```


### Optional Deduplication (Illumina only)

```
nextflow run main.nf --platform illumina --fastq_dir /path/to/fastq/directory --reference /path/to/reference/directory --collapse_duplicate_reads -resume
```


## Testing

To test if the pipeline is working properly the pipeline is provided with small test fastq files. 

   To test the illumina workflow: 
```
nextflow run main.nf -profile test_illumina -resume
```

   To test the nanopore workflow: 
```
nextflow run main.nf -profile test_nanopore -resume
```


## Stopping and Resuming 
- To stop the run
```
control(^) C
```
- To resume the run add the -resume option
```
Nextflow run main.nf --platform --fastq_dir --reference -resume
```


## Dependencies
To run the pipeline the user will need to be working on a computer that has nextflow and singularity installed.

This pipeline requires nextflow version > 24.10.5

There is no specified version of Singularity for this pipeline. The pipeline has been tested with singularity-ce v3.9.9-bionic.