# CNV Project Pipeline

This repository contains a Nextflow pipeline to recreate the CNV analyses in the NSPHS.

## CNV Calling
CNV calling is realized in the folder `cnv_calling`.
The workflow specification is `cnv_calling.nf`.

Run the workflow directory using
```
nextflow run -c nextflow.conf [-resume] [-profile {standard,slurm}] [--bams <alignment-files>] [--reference <reference-folder>] cnv_calling/cnv_calling.nf
```
The configuration and directives are specific to Bianca and SLURM.
The `standard` profile runs the pipeline on the machine you are logged in to.
The `slurm` profile submits each process as a SLURM job.

## CNV Matrix Assembly
After calling, sample-level CNVs and associated copy numbers are aligned to non-overlapping 200-bp windows and assembled into a matrix with rows representing harmonized CNV regions and columns samples.
This matrix is then collapsed by merging adjacent bins with consistent copy numbers.

Run the workflow using
```
nextflow run -c nextflow.conf [-resume] [-profile {standard,slurm}] [--raw_variants <files with raw cnv calls>] [--qc-variants <files with QC'ed cnv calls>] cnv_matrix/cnv_matrix.nf
```