# CNV Project Pipeline

This repository contains a Nextflow pipeline to recreate the CNV analyses in the NSPHS.

## CNV Calling
CNV calling is realized in the folder `cnv_calling`.
The workflow specification is `cnv_calling.nf`.

Run the workflow directory using
```
nextflow run -c nextflow.conf [-profile {standard,slurm}] [--bams <alignment-files>] [--reference <reference-folder>] cnv_calling/cnv_calling.nf
```
The configuration and directives are specific to Bianca and SLURM.
The `standard` profile runs the pipeline on the machine you are logged in to.
The `slurm` profile submits each process as a SLURM job.