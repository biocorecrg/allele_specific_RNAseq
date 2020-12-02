# Allele specific RNAseq pipeline
[![Docker Build Status](https://img.shields.io/docker/automated/biocorecrg/asrnaseq.svg)](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/asrnaseq/builds)
[![Nextflow version](https://img.shields.io/badge/Nextflow-20.01.0-brightgreen)](https://www.nextflow.io/)
[![Singularity version](https://img.shields.io/badge/Singularity-v3.2.1-green.svg)](https://www.sylabs.io/)
[![Docker version](https://img.shields.io/badge/Docker-v19.03-blue)](https://www.docker.com/)
<img align="right" href="https://biocore.crg.eu/" src="https://raw.githubusercontent.com/CRG-CNAG/BioCoreMiscOpen/master/logo/biocore-logo_small.png" />

This **Nextflow**[1] based workflow allows alignment of either single or paired ends reads to a reference transcriptome described as a genome in fasta file plus an annotation in GTF format using **STAR**[2] aligner considering the variant information provided as a VCF file.

**STAR** implements the **WASP**[3] method for filtering of allele specific alignments and reports within the resulting alignments which variant is found and from which allele.
A python script based on **pysam** [4] efficiently separates the alignment generated by STAR into four different files:
- Allele A
- Allele B
- Reference
- Ambiguous (found both variants on one or different pairs)

**HTseq-count**[5] tool is then used for counting tags separated in those three categories mapping to the genes.
Finally a report is generated using **multiQC**[6] that summarize the results of each step together with an initial QC evaluation fo raw reads done using **FastQC**[7]

## Install
You need to install either [Docker](https://docs.docker.com/install/) or [Singularity](https://sylabs.io/guides/3.1/user-guide/installation.html) and [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html). Then you can clone the repository:

```bash
git clone --depth 1 git@github.com:biocorecrg/allele_specific_RNAseq.git
```

or

```
git clone --depth 1 https://github.com/biocorecrg/allele_specific_RNAseq.git 
```

depending on your GitHub configuration.


## Prepare the input data
You need to extract the SNP information from the global VCF file, so first of all you need to download the VCF file with the index:

```bash
wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.normed.vcf.gz
wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.normed.vcf.gz.tbi
```

then the reference genome:
```
wget ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa
```

and the annotation from Ensembl. We used the version Mus_musculus.GRCm38.68 not available in Ensembl archive.

The module **makeAnno** can be used for generating a VCF file with SNP for the interesting species and a genome with SNP position masked with Ns.

For using the module:

```bash
cd makeAnno
nextflow run make_anno.nf -bg --vcffile mgp.v5.merged.snps_all.dbSNP142.vcf.gz --speciesA CAST_EiJ --speciesB 129S1_SvImJ --genome GRCm38_68.fa --outvcf CAST_EiJ-129S1_SvImJ.vcf > log
```

This can take some memory for masking the genome.

## Run the pipeline
```bash
cd allele_specific_RNAseq; 
nextflow run as_rnaseq.nf -bg > log
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
|strandness |reverse|
|variants |$baseDir/test/19_filt.vcf.gz| 
|output |$baseDir/output_test|
|single |NO|
|varcut |1 |
|title |Allele specific RNAseq project|
|subtitle |This is my wonderful RNA experiment|
|PI |Luca Cozzuto|
|User |Luca Cozzuto|
|UCSCgenomeID |mm10|
|email |mymail@mydomain.eu|

providing a real email address will deliver a mail with the multiqc report when the analysis is finished.


### Fastq reads
Fastq paired ends reads can be either plain or gzipped. 
### Genome
Gzipped masked fasta file of the genome obtained running makeAnno
### annotation
GTF file
### variants
Gzipped VCF file obtained running makeAnno
### single
YES: single end reads. NO: paired ends
### varcut
Number of SNP needed for assigning a read to a variant. 



## Results
The following folder will contain the final outputs:
* Index: the indexed genome
* QC: containing fastQC results
* Alignments: containing sorted BAM files as they were from normal bulk RNAseq
* Report: A detailed report for the pipeline run, that will be sent via email.
* Allele_alignments: containing BAM files containing reads with variants marked as: *alleleA*, *alleleB*, *ref* and *ambiguous*
* cut_N *  A folder that is generated for each SNP cut off chosen containing the following sub-folders:
 * Allele_single_counts: containing the count per gene per sample and single alleles
 * Allele_merged_Counts:  containing the count per gene per allele in a single file for each sample
 * Counts: composite counts per gene per sample. No distincion between alleles
 * Proportions: for each sample the read count per allele are divided for the sum of counts of the variants and multiplied for the composite count.

```
propA = alleleA/(alleleA+alleleB) * composite 
propB = alleleB/(alleleA+alleleB) * composite 
```


## References
1. Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319.
1. Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. 
1. van de Geijn B, McVicker G, Gilad Y, Pritchard JK. WASP: allele-specific software for robust molecular quantitative trait locus discovery. Nat Methods. 2015 Nov;12(11):1061-3.
1. https://github.com/pysam-developers/pysam
1. Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8.
1. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
