# Specific README for CRG

Clone the repository:

```bash
git clone --depth 1 git@github.com:biocorecrg/allele_specific_RNAseq.git
```

or

```
git clone --depth 1 https://github.com/biocorecrg/allele_specific_RNAseq.git
```

depending on your GitHub configuration.

## Install New CRG cluster
Singularity is already installed. You need to install nextflow locally and add to the path. 

```
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mkdir -p $HOME/.local/bin/
mv nextflow $HOME/.local/bin/
```

Then change your .bashrc file for adding the right version of JAVA and nextflow

```
vi $HOME/.bashrc
```

type ```i``` for inserting and copy paste this:


``` 
  
  export EASYBUILD_PREFIX=/software/sit/EasyBuild
  module use $EASYBUILD_PREFIX/modules/all
  module load Java/11.0.20
  export PATH=$HOME/.local/bin/:$PATH

```

and exit typing:

```
:wq
```
and return.

For launching the pipeline you need to create a folder in /scratch or /nfs/scratch02 and then pass this via command line:

```
mkdir /scratch/YOURGROUP/workdir

sbatch launch_nf.sh nextflow run as_rnaseq.nf -with-singularity -ansi-log false -profile slurm -w /scratch/YOURGROUP/workdir
```


## Install (OLD CLUSTER)
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

For launching the pipeline you need to create a folder in /nfs/scratch

```
mkdir /scratch/YOURGROUP/workdir

nextflow run as_rnaseq.nf -with-singularity -bg -w /scratch/YOURGROUP/workdir > log.txt
```


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


