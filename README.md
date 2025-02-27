# OrbiSeq
OrbiSeq is a Nextflow pipeline to analyze short and long sequences for Orbiviruses.

## Contents
- [Pipeline Overview](#Pipeline-Overview)
- [Workflow Steps](#Workflow-Steps)
- [Running the Pipeline](#Running-the-Pipeline)
	- [Optional Deduplication](#Optional-Deduplication)
- [Stopping and Resuming](#Stopping-and-Resuming)
- [Run Test Profile](#Run-Test-Profile)

## Pipeline Overview
OrbiSeq is a nextflow pipeline that will create consensus sequences for any Orbivirus with 10 segments from sequencing data produced with Illumina or Nanopore technologies. This pipeline has the capability to run premade reference files included with the pipeline or custom reference files uploaded by the user. 

### Platform:

This pipeline has the capability to run either Illumina or Nanopore sequencing data. When running this pipeline the user must specify which platform their data came from. 

### Reference:

This pipeline can utlitize any reference genome from Orbiviruses with 10 segments. 
- Premade Reference Sequences 
	- There are premade reference files for Bluetongue virus (BTV) and Epizootic hemorrhagic disease virus (EHDV). More information on the creation of these references can be found in [cite publication]. 
		- The reference files for BTV and EHDV can be found in ./reference/BTV and ./reference/EHDV respectively. 
- Custom Reference Sequences	
	- The user is also able to upload their own custom reference file by placing it in the ./reference directory


## Workflow Steps

### Illumina workflow 
- Preprocessing of input reads : Stenglein Lab read_preprocessing pipeline
	- Optional collapse duplicate reads 
- Align input reads to large Orbi RefSeq : bowtie2 build & align 
- Process files : samtools 
- Choose best 10 segments from initial alignment
- Align input reads to best10 RefSeq : bowtie2 build & align 
- Process files : samtools & bcftools
- Create consensus sequence : bcftools consensus
- Align input reads to final consensu sequence : bowtie2 build & align 

### Nanopore workflow 
- Align input reads to large Orbi RefSeq : minimap2 align 
- Process files : samtools 
- Choose best 10 segments from initial alignment
- Align input reads to best10 RefSeq : minimap2 align 
- Process files : samtools & bcftools
- Create consensus sequence : bcftools consensus
- Align input reads to final consensu sequence : minimap2 align 
	
These workflows take advantage of nf-core [modules](https://nf-co.re/modules) for many of these components and the overall [nf-core](https://nf-co.re/) design philosophy.

Additionally, the illumina workflow takes advantage of the [Stenglein Lab Read Preprocessing Pipeline](https://github.com/stenglein-lab/read_preprocessing).

## Running the Pipeline

1. Clone the pipeline from github and move into the directory
```
git clone https://github.com/tdunham19/OrbiSeq.git
cd OrbiSeq
```
2. Create results directory
```
mkdir results
```
3. Run the pipeline: must specify sequencing platform and where the reference file is located.  
```
nextflow run main.nf --platform ['illumina' or 'nanopore'] --fastq_dir /path/to/fastq/directory --reference /path/to/reference/directory  -resume
```

### Optional Deduplication 

```
nextflow run main.nf --platform ['illumina' or 'nanopore'] --fastq_dir /path/to/fastq/directory --reference /path/to/reference/directory --collapse_duplicate_reads -resume
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

## Run Test Profile

To run the test profile utilize the following code
```
nextflow run main.nf --platform ['illumina' or 'nanopore'] --reference ./reference/BTV/BTV_premade_refseq.fasta --profile test -resume
```