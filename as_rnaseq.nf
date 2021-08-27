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
genome                        : ${params.genome}
annotation                    : ${params.annotation}
strandness                    : ${params.strandness}
indexfolder                   : ${params.indexfolder}
variants                      : ${params.variants}
single (YES or NO)            : ${params.single}
varcut						  : ${params.varcut}
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
if (params.strandness != "unstranded" && params.strandness != "forward" && params.strandness != "reverse" ) exit 1, "Please define either forward, reverse or unstranded as strandness parameter\n"

/*
 * Setting the reference genome file and the annotation file (validation)
 */
 
genome_file = file(params.genome)
variants_file = file(params.variants)
annotation_file = file(params.annotation)
multiconfig = file("pre_config.yaml")

if( !genome_file.exists() ) exit 1, "Missing genome file: ${genome_file}"
if( !annotation_file.exists() ) exit 1, "Missing annotation file: ${annotation_file}"
if( !variants_file.exists() ) exit 1, "Missing variations file: ${variants_file}"

subsize          = 5000000
outputfolder    = "${params.output}"
outputQC        = "${outputfolder}/QC"
outputMultiQC   = "${params.output}/Report"
outputMapping   = "${outputfolder}/Alignments"
outputvMapping  = "${outputfolder}/cut_${params.varcut}/Allele_alignments"
trimmedReads    = "${outputfolder}/Trimmed"
outputProp		= "${outputfolder}/cut_${params.varcut}/Proportions"
outputReport    = file("${outputMultiQC}/multiqc_report.html")
outputCounts    = "${outputfolder}/Counts"
outputsCounts   = "${outputfolder}/cut_${params.varcut}/Allele_single_counts"
outputmCounts   = "${outputfolder}/cut_${params.varcut}/Allele_merged_Counts"
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
    .set { reads_for_mapping }


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

    input:
    file(read) from reads_for_fastqc.flatten()

    output:
    file("*_fastqc.*") into raw_fastqc_files

    script:
    """
    fastqc -t ${task.cpus} $read
    """

}

/*
 * Extract read length for being used for indexing
*/
process getReadLength {
    input:
    file(single_read_pairs) from read_files_for_size.first()

    output:
    stdout into (read_length_for_index, read_length_for_profile)

    script:
    """
    	if [ `echo ${single_read_pairs} | grep "gz"` ]; then cat="zcat"; else cat="cat"; fi
    	\$cat ${single_read_pairs} | awk '{num++}{if (num%4==2){line++; sum+=length(\$0)} if (line==100) {printf "%.0f", sum/100; exit} }'
    """
}



/*
 * Builds the genome index required by the mapping process by using STAR aligner
 * Use read size from the estimated one
 */
 
process buildIndex {
    publishDir "${params.indexfolder}"
    label 'big_comp'
    errorStrategy 'retry'
    maxRetries 1
    tag { "${genome_file} with ${annotation_file}" }
    
    input:
    file(genome_file)
    file annotation_file
    // trick for converting a stdout value in integer...
    val read_size from read_length_for_index.map { it.trim().toInteger() }

    output:
    path "*", type: 'dir' into STARgenomeIndex, STARgenomeIndexForCoverage

    script:
    def genome_id = genome_file.simpleName
    def over = read_size-1
    """
    mkdir ${genome_id}
    if [ `echo ${genome_file} | grep ".gz"` ]; then 
        zcat ${genome_file} > `basename ${genome_file} .gz` 
        STAR --runMode genomeGenerate --genomeDir ${genome_id} --runThreadN ${task.cpus} \
        --genomeFastaFiles `basename ${genome_file} .gz` --sjdbOverhang ${over} --sjdbGTFfile ${annotation_file} \
        --outFileNamePrefix `basename ${genome_file} .gz` 
        rm `basename ${genome_file} .gz`
    else 
        STAR --runMode genomeGenerate --genomeDir ${genome_id} --runThreadN ${task.cpus} \
        --genomeFastaFiles ${genome_file} --sjdbOverhang ${over} --sjdbGTFfile ${annotation_file} \
        --outFileNamePrefix ${genome_file} 
    fi
    """
}


/*
 * Align reads by using STAR mapper. 
 * To preserve the order of aligned reads we have to use one thread 
 * and sorting order none
*/

process mapping {
    tag { "${pair_id}" }
    label 'big_comp'
    	publishDir outputCounts, pattern: "${pair_id}/${pair_id}ReadsPerGene.out.tab",  mode: 'copy'
    	publishDir outputMapping, pattern: "${pair_id}/*.out.bam*",  mode: 'copy'
    	//publishDir outputQC, pattern: "STAR_${pair_id}/*Log.final.out", mode: 'copy'

        input:
        set pair_id, file(reads) from reads_for_mapping
        file(STARgenome) from STARgenomeIndex
        file(variants_file)
        
        output:
        set pair_id, file("${pair_id}/${pair_id}ReadsPerGene.out.tab") into STAR_counts
        set pair_id, file("${pair_id}/*.out.bam") into STARmappedBam_for_filtering
        file("${pair_id}/*.out.bam.bai")
        file("${pair_id}") into Aln_folders_for_multiqc

        script:
        def output = "${pair_id}"
        def variants = unzipBash("${variants_file}") 
        """
        if [ `echo "${reads}"| cut -f 1 -d " " | grep ".gz"` ]; then gzipped=" --readFilesCommand zcat "; else gzipped=""; fi
            STAR --genomeDir ${STARgenome} \
                 --readFilesIn ${reads} \
                  \$gzipped \
                  --waspOutputMode SAMtag \
                  --outSAMunmapped Within \
                  --outSAMtype BAM SortedByCoordinate \
                  --runThreadN ${task.cpus} \
                  --outFileNamePrefix ${pair_id} \
                  --quantMode GeneCounts \
                  --outSAMattributes NH HI AS nM NM MD jM jI XS MC ch vA vW vG \
                  --varVCFfile ${variants};
                  mkdir ${output}
                  mv *.out.tab ${output}/
                  mv *Aligned* ${output}/
                  mv *Log* ${output}/
 
        	samtools index ${pair_id}/*.out.bam
        """
        
}

/*
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


/*
 * MultiQC QC. 
*/

process tool_report {

        input:
        file(tooldb)

        output:
        file("tools_mqc.txt") into tool_report_for_multiQC

        script:
        """
        make_tool_desc_for_multiqc.pl -l fastqc,star,samtools,htseq > tools_mqc.txt
        """
}



/*
 * Filter Bam files 
 */
process filterBam {
    publishDir outputvMapping, mode: 'copy'
    tag { pair_id }
    label 'big_comp'

    input:
    set pair_id, file(bamfile) from STARmappedBam_for_filtering

    output:
    set pair_id, file("${pair_id}_*_s.bam") into allele_bams
    set pair_id, file("${pair_id}_*_s.bam.bai") 
    
    script:
    """
	splitBamPerAlleles.py -i ${bamfile} -o ${pair_id} -m ${params.varcut}
	for i in ${pair_id}_*.bam; do samtools sort \$i -o `basename \$i .bam`_s.bam; samtools index `basename \$i .bam`_s.bam; done
    """
}

/*
 * Counts tags with HTSEQ
 */

process countTags {
	
	publishDir outputsCounts, mode: 'copy'
    tag { bamfile }

	input:
	file(annotation_file)
	set pair_id, file(bamfile) from allele_bams.transpose()

	output:
	set pair_id, file("*.counts") into tag_counts_for_grouping       
       
	script:
	def strandness = ""
    if (params.strandness=="unstranded") {
    	strandness = "no"
    } else if (params.strandness=="forward") {
        strandness = "yes"
    } else if (params.strandness=="reverse") {
         strandness = "reverse"
    }

    """
    htseq-count -s ${strandness} -f bam ${bamfile} ${annotation_file} > `basename ${bamfile} .bam`.counts
    """
}



/*
 * Group counts tags 
 */

process groupCounts {
	
	publishDir outputmCounts, mode: 'copy', pattern: "*_group.counts"
    tag { pair_id }

	input:
	file(annotation_file)
	set pair_id, file("*") from tag_counts_for_grouping.groupTuple()

	output:
	set pair_id, file("*_group.counts") into count_group_for_proportion
	file("*_group.stats") into counts_stats_for_report
       

    script:
    """
	join_counts.sh ${pair_id}
    """
}

/*
 * Make proportions
 */

process makePropoportions {
	tag { pair_id }
	publishDir outputProp, mode: 'copy'

	input:
	set pair_id, file(group_count), file(read_count) from count_group_for_proportion.join(STAR_counts)
	file(annotation_file)
       
	output:
	file("${pair_id}_composite.txt")

    script:
	def strandness = ""
    if (params.strandness=="unstranded") {
    	strandness = "u"
    } else if (params.strandness=="forward") {
        strandness = "f"
    } else if (params.strandness=="reverse") {
         strandness = "r"
    }
    """
	calc_prop.py -a ${group_count} -c ${read_count} -g ${annotation_file} -s "${strandness}" -o ${pair_id}_composite.txt
    """
}

/*
* 
*/

process multiQC_report {
    publishDir outputMultiQC,  mode: 'copy'

    input:
    file (multiconfig)
    file 'tools_mqc.txt' from tool_report_for_multiQC
    file '*' from Aln_folders_for_multiqc.mix(raw_fastqc_files, counts_stats_for_report).flatten().collect()
    
    output:
    file("multiqc_report.html") into multiQC 

    script:
    """
    export LC_ALL=en_US.utf8
    export LANG=en_US.utf8
    make_conf_multiqc.sh \'${params.title}\' \'${params.subtitle}\' \'${params.PI}\' \'${params.User}\' \'${params.email}\' \'${UCSCgenomeID}\' ${multiconfig} 
    
    cat > counts_mqc.txt << EOL
# id: read_counts
# plot_type: 'bargraph'
# section_name: 'Read counts on alleles'
Sample	Genotype A	Genotype B	Reference	Ambiguous
EOL
    grep -h -v "#" *.stats >> counts_mqc.txt
    multiqc -c config.yaml .
    """

}


// make named pipe 
def unzipBash(filename) { 
    def cmd = filename.toString()
    if (cmd[-3..-1] == ".gz") {
    	cmd = "<(zcat ${filename})"
    }
    return cmd
}

// extract unzip pipe 
def getUnzipName(filename) { 
    def name = filename.toString()
    if (name[-3..-1] == ".gz") {
    	name = name.baseName
    }
    return name
}


workflow.onComplete {
    println "Pipeline BIOCORE@CRG Master of Pore completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
