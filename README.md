# Allele specific RNAseq pipeline
[![Docker Build Status](https://img.shields.io/docker/automated/biocorecrg/asrnaseq.svg)](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/asrnaseq/builds)
[![Nextflow version](https://img.shields.io/badge/Nextflow-20.10.0-brightgreen)](https://www.nextflow.io/)
[![Singularity version](https://img.shields.io/badge/Singularity-v3.2.1-green.svg)](https://www.sylabs.io/)
[![Docker version](https://img.shields.io/badge/Docker-v19.03-blue)](https://www.docker.com/)
<img align="right" href="https://biocore.crg.eu/" src="https://raw.githubusercontent.com/CRG-CNAG/BioCoreMiscOpen/master/logo/biocore-logo_small.png" />

This **Nextflow**[1] based workflow allows alignment of single ends reads to a reference genome or transcriptome (described as a genome in fasta file plus an annotation in GTF format) using either **BWA** or **STAR**.

## Install
You need to install either [Docker](https://docs.docker.com/install/) or [Singularity](https://sylabs.io/guides/3.1/user-guide/installation.html) and [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html). 


Then you can clone the repository:

```bash
git clone --depth 1 git@github.com:biocorecrg/RIPseq.git
```

or

```
git clone --depth 1 https://github.com/biocorecrg/RIPseq.git 
```

depending on your GitHub configuration.


## Input data
The input data required by the pipeline are:
1) (Gzipped) Fastq reads
2) (Gzipped) Reference genome in fasta format
3) (Gzipped) GTF file with annotations

## Configuration files
1) A tsv file containing three columnes indicating
```
SAMPLEID  CONTROLID PEAKCALLER_PROGRAM_NAME
```
See **peakconf.tsv** as an example.

2) The **tool_opt.tsv** file containing information about custom options to be used by mappers and peak callers:
```
#tool	chip_pars	rip_pars
bwa     "-M"    ""
star	""	"--quantMode GeneCounts"
macs2	""	"--nomodel --extsize 100"
epic2	""	""
moaims	""	"strand_specific=2"
```

## Run the pipeline
```bash 
nextflow nextflow run ripseq.nf -bg > log
```

optionally you might want to check your pipeline on Nextflow's [Tower](https://tower.nf/) website. You need to register using an istitutional mail and set the token provided in a variable as described:

```bash
export TOWER_ACCESS_TOKEN=<<<<<TOKEN NUMBER>>>>>>
nextflow run as_rnaseq.nf -bg -with-tower > log 
```

you can check the status of your pipeline live on the tower website.

The parameters for running the pipeline are defined in the file **params.config** that can be changed accordingly.

|parameter|value|
|:---:|:---:|
|reads |$baseDir/test/*_{1,2}.fastq.gz|
|genome |$baseDir/test/GRCm38_68_19.masked.fa.gz|
|annotation |$baseDir/test/Mus_musculus.GRCm38.68_19.gtf|
|tool_opts |$baseDir/tool_opt.tsv|
|peakconfig |$baseDir/peakconf.tsv| 
|type |ripseq|
|output |$baseDir/output|


## Results
The following folders will contain the final outputs:
* QC: containing fastQC results
* Alignments: containing sorted BAM files
* multiQC: A detailed report for the pipeline run.
* Peaks: containing three subfolders with the results from each single peak caller together with the intersection **NAME_multiinter.bed** and the venn diagram **NAME_venn.png**   

