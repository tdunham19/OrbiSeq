# OrbiSeq
OrbiSeq is a Nextflow pipeline to analyze short and long sequences for Orbiviruses.

## Contents
- [Workflow Steps](#Workflow-Steps)
- [Running the Pipeline](#Running-the-Pipeline)
- [Stopping and Resuming](#Stopping-and-Resuming)

## Workflow Steps
- Illumina workflow 
	- Preprocessing of input reads : read_preprocessing MDS pipeline
	- Align input reads to large Orbi RefSeq : bowtie2 build & align 
	- Process files : samtools 
	- Choose best 10 segments from initial alignment
	- Align input reads to best10 RefSeq : bowtie2 build & align 
	- Process files : samtools & bcftools
	- Create consensus sequence : bcftools consensus
	- Align input reads to final consensu sequence : bowtie2 build & align 

- Nanopore workflow 
	- Align input reads to large Orbi RefSeq : minimap2 align 
	- Process files : samtools 
	- Choose best 10 segments from initial alignment
	- Align input reads to best10 RefSeq : minimap2 align 
	- Process files : samtools & bcftools
	- Create consensus sequence : bcftools consensus
	- Align input reads to final consensu sequence : minimap2 align 
	
This takes advantage of nf-core [modules](https://nf-co.re/modules) for many of these components and the overall [nf-core](https://nf-co.re/) design philosophy.

## Running the Pipeline

1. Clone the pipeline from github and move into the directory
```
git clone https://github.com/tdunham19/OrbiSeq.git
cd read_preprocessing
```
2. Create fastq, reference, and results directories 
```
mkdir fastq
mkdir reference
mkdir results
```
3. Run the pipeline with sequencing platform specified
```
# Illumina workflow
nextflow run main.nf --platform illumina
# Nanopore workflow
nextflow run main.nf --platform nanopore
```

## Stopping and Resuming 
- To stop the run
```
control(^) C
```
- To resume the run add the -resume option
```
Nextflow run main.nf --platform specify -resume
```
