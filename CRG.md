# Specific README for CRG

## Install
Singularity is already installed. You need to add to log to nextflow node:
```
ssh -Y nextflow.linux.crg.es
```

Then change your .bashrc file for adding the right version of singularity:

```
vi $HOME/.bashrc
```

type ```i``` for inserting and copy paste this:


```
module use /software/as/el7.2/EasyBuild/CRG/modules/all
module load Singularity/3.2.1
```
and exit typing:

```
:wq
```
and return.

Then you can clone the repository:

```bash
git clone --depth 1 git@github.com:biocorecrg/allele_specific_RNAseq.git
```

or

```
git clone --depth 1 https://github.com/biocorecrg/allele_specific_RNAseq.git 
```

depending on your GitHub configuration.


## Annotation data
Annotation data are available at:
```
/users/bpayer/sequencing_analysis/pipe_data
```
You have:

* SNPs annotations: CAST_EiJ-129S1_SvImJ.vcf.gz 
* Mouse reference genome: Mus_musculus.GRCm38.68.dna.chrom.fa.gz
* Mouse annotation: Mus_musculus.GRCm38.68.gtf.gz

The parameters for running the pipeline are defined in the file **params.crg.config**. You can copy this file and renaming it params.config and change it accordingly. 

In particular you need to change only the following parameters:  

```
	reads        = "/test/*_{1,2}.fastq.gz"
	strandness   = "reverse"
  output       = "$baseDir/output_test"
	single       = "NO"
	varcut       = 1
  title	     = "Allele specific RNAseq project"	
	subtitle     = "This is my wonderful RNA experiment"	
	PI           = "Luca Cozzuto"	
	User	     = "Luca Cozzuto"
	email	     = "mymail@mydomain.eu"	
```


