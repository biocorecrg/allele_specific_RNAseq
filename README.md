# Allele specific RNAseq pipeline
[![Docker Build Status](https://img.shields.io/docker/automated/biocorecrg/asrnaseq.svg)](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/asrnaseq/builds)
[![Nextflow version](https://img.shields.io/badge/Nextflow-19.10.0-brightgreen)](https://www.nextflow.io/)
[![Singularity version](https://img.shields.io/badge/Singularity-v2.6.1-green.svg)](https://www.sylabs.io/)
[![Docker version](https://img.shields.io/badge/Docker-v19.03-blue)](https://www.docker.com/)
<img align="right" href="https://biocore.crg.eu/" src="https://raw.githubusercontent.com/CRG-CNAG/BioCoreMiscOpen/master/logo/biocore-logo_small.png" />

The pipeline is based on Wasp implementation[1] in STAR[2] aligner. 
In brief the input paired ends reads are inspected per quality by using FastQC[3] and aligned to the reference transcriptome (i.e the genome in fasta file and the annotation in GTF file) using STAR with the option **--waspOutputMode**.
This option will activate the WASP filtering of allele specific alignments[1]. The current version of the aligner is able to use the variant information stored in a VCF file for adding new tags to the output with information about which allele is detected in the read. 

We used this information for assigning a read to:


- Allele A (only variants from allele A or undetermined)
- Allele B (only variants from allele B or undetermined)
- Ambiguous (found both variants from allele A and B on one or both pairs)
- Undetermined (only undetermined variants)

## References
1. van de Geijn B, McVicker G, Gilad Y, Pritchard JK. WASP: allele-specific software for robust molecular quantitative trait locus discovery. Nat Methods. 2015 Nov;12(11):1061-3.
1. Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. 
1. 
