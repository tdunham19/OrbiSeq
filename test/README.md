## Test datasets

This directory contains small test datasets from BTV17 samples sequenced side-by-side on Illumina and Nanopore platforms. 

Additionally metadata for the samples is included here with accession numbers from when this was originally sequenced in 2018: [test metadata](./test_metadata.xlsx)

These datasets include:

| Dataset | Platform | Description |
| ------- | -------- | ------- |
| BTV17_illumina_highcov_subset_[R1-2]      | Illumina | A short read dataset from BTV17 CO2018 Cell Culture Isolate, ~900,000 total reads |
| BTV17_illumina_0.1_highcov_subset_[R1-2]  | Illumina | A short read dataset from BTV17 CO2018 Cell Culture Isolate, 1/10th the amount of reads from highcov_subset  |
| BTV17_illumina_0.01_highcov_subset_[R1-2] | Illumina | A short read dataset from BTV17 CO2018 Cell Culture Isolate, 1/100th the amount of reads from highcov_subset |
| BTV17_illumina_lowcov_subset_[R1-2]       | Illumina | A short read low coverage dataset from BTV17 CO2018 Cell Culture Isolate, pipeline should not fail |
| BTV17_nanopore_highcov_subset             | Nanopore | A long read dataset from BTV17 CO2018 Cell Culture Isolate, 5000 total reads  |
| BTV17_nanopore_10th_highcov_subset        | Nanopore | A long read dataset from BTV17 CO2018 Cell Culture Isolate, 500 total reads   |
| BTV17_nanopore_100th_highcov_subset       | Nanopore | A long read dataset from BTV17 CO2018 Cell Culture Isolate, 50 total reads    |

These datasets were created using seqtk and subseting for 1000 reads. 

These datasets are used as input when using the test profile:

```
nextflow run main.nf -resume -profile singularity,test_illumina --platform illumina 
nextflow run main.nf -resume -profile singularity,test_nanopore --platform nanopore 
```