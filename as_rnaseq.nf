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
reads                         : ${params.reads}
genomeA                       : ${params.genomeA}
genomeB                       : ${params.genomeB}
annotation                    : ${params.annotation}
single (YES or NO)            : ${params.single}
output (output folder)        : ${params.output}
UCSCgenomeID genome ID 
for UCSC hub (optional)       : ${params.UCSCgenomeID}
title                         : ${params.title}
subtitle                      : ${params.subtitle}
PI principal investigator     : ${params.PI}
User                          : ${params.User}
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
 
genome_fileA = file(params.genomeA)
genome_fileB = file(params.genomeB)
annotation_file = file(params.annotation)
multiconfig = file("config.yaml")

if( !genome_fileA.exists() ) exit 1, "Missing genome file: ${genome_fileA}"
if( !genome_fileB.exists() ) exit 1, "Missing genome file: ${genome_fileB}"
if( !annotation_file.exists() ) exit 1, "Missing annotation file: ${annotation_file}"

subsize          = 5000000
outputfolder    = "${params.output}"
outputQC        = "${outputfolder}/QC"
outputMultiQC   = "${params.output}/multiQC"
outputMapping   = "${outputfolder}/Alignments"
outputIndex     = "${outputfolder}/Index"
trimmedReads    = "${outputfolder}/Trimmed"
outputReport    = file("${outputMultiQC}/multiqc_report.html")
outputCounts    = "${outputfolder}/Counts"
UCSCgenomeID    = "${params.UCSCgenomeID}"
tooldb          = file("conf_tools.txt")
if( UCSCgenomeID == "" ) {        UCSCgenomeID = "custom"    }
rootProfiles = "${outputfolder}/Profiles"
outputProfiles = "${rootProfiles}/${UCSCgenomeID}"

/*
* move old multiQCreport
*/

if( outputReport.exists() ) {
  log.info "Moving old report to multiqc_report.html multiqc_report.html.old"
  outputReport.moveTo("${outputMultiQC}/multiqc_report.html.old")
}

Channel
    .fromPath( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_files_for_size; reads_for_fastqc}    


/*
 * For paired-ends experiments
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 * For single-end 
 * Creates the 
 */
Channel
    .fromFilePairs( params.reads , size: (params.single != "NO") ? 1 : 2)
    .into { raw_reads_for_mapping; raw_reads_for_trimming }


/*
 * Mail notification
*/
 
workflow.onComplete {
    def subject = 'Allele specific RNAseq execution'
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



/*
 * Run FastQC on raw data
*/
process QConRawReads {
    tag { read }
    publishDir outputQC
	afterScript 'mv *_fastqc.zip `basename *_fastqc.zip _fastqc.zip`_raw_fastqc.zip'

    input:
    file(read) from reads_for_fastqc.flatten()

    output:
    file("*_fastqc.*") into raw_fastqc_files

    script:
    """
    fastqc -t ${task.cpus} $read
    mv *_fastqc.zip `basename *_fastqc.zip _fastqc.zip`_raw_fastqc.zip
    """

}

/*
 * Extract read length for being used for indexing
*/
process getReadLength {
    input:
    file(single_read_pairs) from read_files_for_size.first()

    output:
    stdout into (read_length_for_trimming, read_length_for_index, read_length_for_profile)

    script:
    """
    	if [ `echo ${single_read_pairs} | grep "gz"` ]; then cat="zcat"; else cat="cat"; fi
    	\$cat ${single_read_pairs} | awk '{num++}{if (num%4==2){line++; sum+=length(\$0)} if (line==100) {printf "%.0f", sum/100; exit} }'
    """
}


Channel.from( ["${genome_fileA.getSimpleName()}", genome_fileA], ["${genome_fileB.getSimpleName()}", genome_fileB] ).set{genomes}

/*
 * Builds the genome index required by the mapping process by using STAR aligner
 * Use read size from the estimated one
 */
 
process buildIndex {
    publishDir outputIndex
    label 'big_comp'
    tag { "${genome_file} with ${annotation_file}" }
    
    input:
    set genome_id, file(genome_file) from genomes
    file annotation_file
    // trick for converting a stdout value in integer...
    val read_size from read_length_for_index.map { it.trim().toInteger() }

    output:
    file genome_id into STARgenomeIndex, STARgenomeIndexForCoverage

    script:
        """
        mkdir ${genome_id}
        if [ `echo ${genome_file} | grep ".gz"` ]; then 
            zcat ${genome_file} > `basename ${genome_file} .gz` 
            STAR --runMode genomeGenerate --genomeDir ${genome_id} --runThreadN ${task.cpus} \
            --genomeFastaFiles `basename ${genome_file} .gz` --sjdbOverhang ${read_size} --sjdbGTFfile ${annotation_file} \
            --outFileNamePrefix `basename ${genome_file} .gz` 
            rm `basename ${genome_file} .gz`
        else 
            STAR --runMode genomeGenerate --genomeDir ${genome_id} --runThreadN ${task.cpus} \
            --genomeFastaFiles ${genome_file} --sjdbOverhang ${read_size} --sjdbGTFfile ${annotation_file} \
            --outFileNamePrefix ${genome_file} 
        fi
        """
}




/*
 * Align reads by using STAR mapper. 
*/

process mapping {
    label 'big_comp'
    tag { pair_id }
    publishDir outputCounts, pattern: "STAR_${pair_id}/*ReadsPerGene.out.tab",  mode: 'copy'
    publishDir outputQC, pattern: "STAR_${pair_id}/*Log.final.out", mode: 'copy'

        input:
        file STARgenome from STARgenomeIndex
        set pair_id, file(reads) from reads_for_mapping_clean

        output:
        set pair_id, file("STAR_${pair_id}/${pair_id}Aligned.sortedByCoord.out.bam") into STARmappedBam_for_qualimap, STARmappedBam_for_indexing
        file("STAR_${pair_id}") into Aln_folders_for_multiqc
        file("STAR_${pair_id}/${pair_id}ReadsPerGene.out.tab")

        script:
        def aligner = new NGSaligner(id:pair_id, reads:reads, index:STARgenome, cpus:task.cpus, output:"STAR_${pair_id}") 
        aligner.doAlignment("STAR")  
        
}



 * Step 7. QualiMap QC.

process qualimap {
    label 'big_comp'
    tag { pair_id }

    input:
    file annotation_file

    set pair_id, file(bamfile) from STARmappedBam_for_qualimap

    output:
    file("QUALIMAP_${pair_id}") into QualiMap_for_multiQC 

    script:
    def qualimp_memory = Math.round(0.6*task.memory.giga ) + "G"
    def qualimode = (params.single != "NO") ? 'pe' : 'se'  

    def qc = new QualityChecker(input:bamfile, annotation_file:annotation_file, output:"QUALIMAP_${pair_id}", strand:"strand-specific-reverse", memory:qualimp_memory, mode:qualimode, extrapars:"-s")
    qc.qualimapRNAseq()  
 }


 * MultiQC QC. 


process tool_report {

        input:
        file(tooldb)

        output:
        file("tools_mqc.txt") into tool_report_for_multiQC

        script:
        """
        make_tool_desc_for_multiqc.pl -l fastqc,star,skewer,qualimap,bedtools,samtools > tools_mqc.txt
        """
}


process multiQC_report {
    publishDir outputMultiQC,  mode: 'copy'

    input:
    file 'pre_config.yaml.txt' from multiconfig
    file 'tools_mqc.txt' from tool_report_for_multiQC
    
    file '*' from Aln_folders_for_multiqc.mix(raw_fastqc_files, trimmed_fastqc_files,logTrimming_for_QC, QualiMap_for_multiQC).flatten().collect()
    
    output:
    file("multiqc_report.html") into multiQC 

    script:
    def reporter = new Reporter(title:params.title, application:"RNA-seq", subtitle:params.subtitle, PI:params.PI, user:params.User, id:UCSCgenomeID, email:params.email,config_file:multiconfig)
    reporter.makeMultiQCreport()
}



 * Index Bam files 

process indexBam {
    publishDir outputMapping, mode: 'copy'
    tag { pair_id }

    input:
    set pair_id, file(bamfile) from STARmappedBam_for_indexing

    output:
    set pair_id, file(bamfile), file("${bamfile}.bai") into bamFiles 
    
    script:
    def misc = new Misc(input:bamfile)
    misc.st_indexBam()
}


 * Convert to BedGraph


process convertToBigWig {
    tag { pair_id }
    publishDir outputProfiles, mode: 'copy' 

    input:
    set pair_id, file(bamfile), file(indexfile) from bamFiles
    file(starindex) from STARgenomeIndexForCoverage
    val read_size from read_length_for_profile.map { it.trim().toInteger() }

    output:
    set pair_id, file("${pair_id}.bw") into bigWig_profiles
    
    script:
    def misc = new Misc(input:bamfile, read_size:read_size, chr_size_file:"${starindex}/chrNameLength.txt", output:"${pair_id}.bw")
    misc.makeAlnProfiles()
    
}


 * make UCSC genome Hub

    process makeGenomeHub {
        publishDir rootProfiles, mode: 'copy'

        when:
        UCSCgenomeID != ""
    
        output:
        file "genomes.txt"
        file "hub.txt"
        
        script: 
        def reporter = new Reporter(id:UCSCgenomeID, title:params.title, subtitle:params.subtitle, email:params.email)
        reporter.makeGenomeUcscHub()

    }

    process makeTrackDB {
        publishDir outputProfiles, mode: 'copy'

        when:
        UCSCgenomeID != ""
        
        input:
        file "*" from bigWig_profiles.collect()

        output:
        file "trackDb.txt"
        
        script: 
        def reporter = new Reporter(id:"composite1", type:"bigWig", extension:"bw", , title:params.title, subtitle:params.subtitle, email:params.email)
        reporter.makeBigWigTrackDB()        
    }
 */
 
workflow.onComplete {
    println "Pipeline BIOCORE@CRG Master of Pore completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}