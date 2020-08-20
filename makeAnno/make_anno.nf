#!/usr/bin/env nextflow

/*
 * Copyright (c) 2020, Centre for Genomic Regulation (CRG) and the authors.
 *
 *
 */

/* 
 * Allele specific RNASeq pipeline from Bioinformatics Core @ CRG
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 *
 * 
 */

// Pipeline version
version = '1.0'

params.help            = false
params.resume          = false

log.info """
BIOCORE@CRG RNAseq - N F  ~  version ${version}

╔═╗┬  ┬  ┌─┐┬  ┌─┐  ╔═╗┌─┐┌─┐┌─┐┬┌─┐┬┌─┐  ╦═╗╔╗╔╔═╗┌─┐┌─┐┌─┐ 
╠═╣│  │  ├┤ │  ├┤   ╚═╗├─┘├┤ │  │├┤ ││    ╠╦╝║║║╠═╣└─┐├┤ │─┼┐
╩ ╩┴─┘┴─┘└─┘┴─┘└─┘  ╚═╝┴  └─┘└─┘┴└  ┴└─┘  ╩╚═╝╚╝╩ ╩└─┘└─┘└─┘└
                                                                                       
====================================================
BIOCORE@CRG Allele Specific RNAseq - N F  ~  version ${version}
====================================================
vcffile                       : ${params.vcffile}
vcfindex                      : ${params.vcfindex}
speciesA                      : ${params.speciesA}
speciesB                      : ${params.speciesB}
genome                        : ${params.genome}
outvcf                        : ${params.outvcf}
email for notification        : ${params.email}

"""

log.info "\n"
if (params.help) {
    log.info 'This is the Biocore\'s Allele Specific RNAseq pipeline'
    log.info '\n'
    exit 1
}

if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

/*
 * Setting the reference genome file and the annotation file (validation)
 */
 
genome_file = file(params.genome)
variants_file = file(params.vcffile)
variants_index = file(params.vcfindex)

if( !genome_file.exists() ) exit 1, "Missing genome file: ${genome_file}"
if( !variants_file.exists() ) exit 1, "Missing variations file: ${variants_file}"
if( !variants_index.exists() ) exit 1, "Missing variations index: ${variants_index}"
if( !params.speciesA ) exit 1, "Missing speciesA parameter"
if( !params.speciesB ) exit 1, "Missing speciesB parameter"

/*
 * Run FastQC on raw data
*/
process parseVCF {
    tag { variants_file }
    publishDir "."

    input:
    file(genome_file)
    file(variants_file)
    file(variants_index)
    
    output:
    file(params.outvcf) 

    script:
    """
    parseVCF.py -i ${variants_file} -o ${params.outvcf} -1 ${params.speciesA} -2 ${params.speciesB} -g ${genome_file}
    """

}



/*
 * Mail notification
*/
 
workflow.onComplete {
    def subject = 'Make annotation pipeline execution'
    def recipient = "${params.email}"
    def attachment = "${outputMultiQC}/multiqc_report.html"

    ['mail', '-s', subject, '-a', attachment, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}




workflow.onComplete {
    println "Pipeline BIOCORE@CRG Master of Pore completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
