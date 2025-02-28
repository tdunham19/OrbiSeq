## Test datasets

This directory contains small test datasets from BTV17 samples sequenced side-by-side on Illumina and Nanopore platforms. 

Additionally metadata for the samples is included here with accession numbers from when this was originally sequenced in 2018: [minimal metadata](./test_metadata.xlsx)

These datasets include:

| Dataset | Platform | Description |
| ------- | -------- | ------- |
| BTV17_illumina_subset | Illumina | BTV17 CO2018 |
| BTV17_nanopore_subset | Illumina | BTV17 CO2018 |

These datasets were created using seqtk and subseting for 1000 reads. 

These datasets are used as input when using the test profile:

```
nextflow run main.nf -profile test_illumina -resume
nextflow run main.nf -profile test_nanopore -resume
```