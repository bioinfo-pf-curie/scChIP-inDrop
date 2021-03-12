#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2020
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/

/*
========================================================================================
<!-- TODO - Pipeline Name -->
========================================================================================
 #### Homepage / Documentation
<!-- TODO - Pipeline code url -->
----------------------------------------------------------------------------------------
*/

// File with text to display when a developement version is used
devMessageFile = file("$baseDir/assets/devMessage.txt")

def helpMessage() {
  if ("${workflow.manifest.version}" =~ /dev/ ){
     log.info devMessageFile.text
  }

  log.info """
  v${workflow.manifest.version}
  ======================================================================

  Usage:
  nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --genome 'hg19' -profile conda
  nextflow run main.nf --samplePlan samplePlan --genome 'hg19' -profile conda

  Mandatory arguments:
    --reads [file]                Path to input data (must be surrounded with quotes)
    --samplePlan [file]           Path to sample plan input file (cannot be used with --reads)
    --genome [str]                Name of genome reference
    -profile [str]                Configuration profile to use. test / conda / multiconda / path / multipath / singularity / docker / cluster (see below)
  
  Inputs:
    --design [file]               Path to design file for extended analysis  
    --singleEnd [bool]            Specifies that the input is single-end reads

  Skip options: All are false by default
    --skipSoftVersion [bool]      Do not report software version
    --skipMultiQC [bool]          Skips MultiQC

  Other options:
    --outDir [file]               The output directory where the results will be saved
    -name [str]                   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
 
  =======================================================
  Available Profiles

    -profile test                Set up the test dataset
    -profile conda               Build a single conda for with all tools used by the different processes before running the pipeline
    -profile multiconda          Build a new conda environment for each tools used by the different processes before running the pipeline
    -profile path                Use the path defined in the configuration for all tools
    -profile multipath           Use the paths defined in the configuration for each tool
    -profile docker              Use the Docker containers for each process
    -profile singularity         Use the singularity images for each process
    -profile cluster             Run the workflow on the cluster, instead of locally

  """.stripIndent()
}

/**********************************
 * SET UP CONFIGURATION VARIABLES *
 **********************************/

// Show help message
if (params.help){
  helpMessage()
  exit 0
}

// Configurable reference genomes

// TODO - ADD HERE ANY ANNOTATION

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes.config file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
customRunName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  customRunName = workflow.runName
}

//Stage config files
Channel
  .fromPath(params.multiqcConfig, checkIfExists: true)
  .set{chMultiqcConfig}
chOutputDocs = file("$baseDir/docs/output.md", checkIfExists: true)
chOutputDocsImages = file("$baseDir/docs/images/", checkIfExists: true)

/************
 * CHANNELS *
 ************/

// Validate inputs
if ((params.reads && params.samplePlan) || (params.readPaths && params.samplePlan)){
  exit 1, "Input reads must be defined using either '--reads' or '--samplePlan' parameters. Please choose one way."
}

//------- Metadata ---------
//--------------------------

if ( params.metadata ){
  Channel
    .fromPath( params.metadata )
    .ifEmpty { exit 1, "Metadata file not found: ${params.metadata}" }
    .set { chMetadata }
}else{
  metadataCh = Channel.empty()
}

//------- Genomes ---------
//-------------------------
genomeRef = params.genome

params.starIndex = genomeRef ? params.genomes[ genomeRef ].starIndex ?: false : false
if (params.starIndex){
  Channel
    .fromPath(params.starIndex, checkIfExists: true)
    .ifEmpty {exit 1, "STAR index file not found: ${params.starIndex}"}
    .set { chStar }
} else {
  exit 1, "STAR index file not found: ${params.starIndex}"
}

params.gtf = genomeRef ? params.genomes[ genomeRef ].gtf ?: false : false
if (params.gtf) {
  Channel
    .fromPath(params.gtf, checkIfExists: true)
    .into { chGtfSTAR; chGtfFC }
}else {
  exit 1, "GTF annotation file not not found: ${params.gtf}"
}

params.bed12 = genomeRef ? params.genomes[ genomeRef ].bed12 ?: false : false
if (params.bed12) {
  Channel  
    .fromPath(params.bed12)
    .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
    .set { chBedGeneCov } 
}else {
  exit 1, "BED12 annotation file not not found: ${params.bed12}"
}

//------- Custom --------
//-----------------------
for ( idx in params.barcodes.keySet() ){
  if ( params.barcodes[ idx ].bwt2 ){
    lastPath = params.barcodes[ idx ].bwt2.lastIndexOf(File.separator)
    bt2Dir = params.barcodes[ idx ].bwt2.substring(0,lastPath+1)
    bt2Base = params.barcodes[ idx ].bwt2.substring(lastPath+1)
    params.barcodes[ idx ].base = bt2Base
    params.barcodes[ idx ].dir = bt2Dir
  }
}

Channel
   .from(params.barcodes)
   .flatMap()
   .map { it -> [ it.key, file(it.value['dir']) ] }
   .ifEmpty { exit 1, "Bowtie2 index not found" }
   .set { chIndexBwt2 } 


//------- Reads ---------
//-----------------------

// Create a channel for input read files
if(params.samplePlan){
  if(params.singleEnd){
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2])]] }
      .into { chRawReadsBowtie2; chRawReadsFastx; chAlignment }
  }else{
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
      .into { chRawReadsBowtie2; chRawReadsFastx; chAlignment}
   }
  params.reads=false
}
else if(params.readPaths){
  if(params.singleEnd){
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied." }
      .into { chRawReadsBowtie2; chRawReadsFastx ; chAlignment}
  } else {
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied." }
      .into { chRawReadsBowtie2 ; chRawReadsFastx; chAlignment}
  }
} else {
  Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { chRawReadsBowtie2 ; chRawReadsFastx; chAlignment}
}

// Make sample plan if not available
if (params.samplePlan){
  chSplan = Channel.fromPath(params.samplePlan)
}else if(params.readPaths){
  if (params.singleEnd){
    Channel
      .from(params.readPaths)
      .collectFile() {
        item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }
      .set{ chSplan }
  }else{
     Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
       .set{ chSplan }
  }
}else{
  if (params.singleEnd){
    Channel
      .fromFilePairs( params.reads, size: 1 )
      .collectFile() {
        item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }
      .set { chSplan }
  }else{
    Channel
      .fromFilePairs( params.reads, size: 2 )
      .collectFile() {
        item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
      }
      .set { chSplan }
   }
}

/***************
 * Design file *
 ***************/

// TODO - UPDATE BASED ON YOUR DESIGN

if (params.design){
  Channel
    .fromPath(params.design)
    .ifEmpty { exit 1, "Design file not found: ${params.design}" }
    .into { chDesignCheck }

  chDesignControl
    .splitCsv(header:true)
    .map { row ->
      if(row.CONTROLID==""){row.CONTROLID='NO_INPUT'}
      return [ row.SAMPLEID, row.CONTROLID, row.SAMPLENAME, row.GROUP, row.PEAKTYPE ]
     }
    .set { chDesignControl }
}else{
  chDesignCheck = Channel.empty()
}

/*******************
 * Header log info *
 *******************/

if ("${workflow.manifest.version}" =~ /dev/ ){
   log.info devMessageFile.text
}

log.info """=======================================================

workflow v${workflow.manifest.version}
======================================================="""
def summary = [:]

summary['Max Memory']     = params.maxMemory
summary['Max CPUs']       = params.maxCpus
summary['Max Time']       = params.maxTime
summary['Container Engine'] = workflow.containerEngine
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


// TODO - ADD YOUR NEXTFLOW PROCESS HERE

process bcMapping {
  tag "${prefix} - ${index}"
  label 'bowtie2'
  label 'highCpu'
  label 'highMem'

  publishDir "${params.outDir}/bcMapping", mode: 'copy'

  input:
  set val(prefix), file(reads), val(index), file(bwt2Idx) from chRawReadsBowtie2.combine(chIndexBwt2)
 
  output:
  set val(prefix), file("*ReadsMatching.txt") into chBarcodeMapped
  set val(prefix), file("*_bowtie2.log") into chBowtie2Logs
  file ("v_bowtie2.txt") into chBowtie2version

  script:
  start = params.barcodes[ index ].start
  size = params.barcodes[ index ].size
  base = params.barcodes[ index ].base
  oprefix = "${prefix}_${index}"

  """
  ##Extract read ids and their index sequences
  gzip -cd  ${reads[1]} | awk -v startIndex=${start} -v sizeIndex=${size} 'NR%4==1{print \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,startIndex,sizeIndex)}' > ${oprefix}Reads.fasta
  #Map indexes (-f) against Index libraries (-x)
  bowtie2 \
    -x ${bwt2Idx}/${base} \
    -f ${oprefix}Reads.fasta \
    -N 1 -L 8 \
    --rdg 0,7 --rfg 0,7 --mp 7,7 \
    --ignore-quals --score-min L,0,-1 \
    -t --no-hd \
    --reorder \
    -p ${task.cpus} > ${oprefix}Bowtie2.sam 2> ${oprefix}_bowtie2.log
 
  ## Keep all alignments, so that ReadsMatching files always have the same size
  # \$1 = read ID and \$3 = barcode number
  awk '/XS/{\$3="*"} {print \$1,\$3}' ${oprefix}Bowtie2.sam > ${oprefix}ReadsMatching.txt

  bowtie2 --version > v_bowtie2.txt

  """
}

/*
process fastxTrimmer {
  tag "${prefix}"
  label 'fastx'
  label 'medCpu'
  label 'medMem'

  publishDir "${params.outDir}/fastxTrimmer", mode: 'copy'

  input:
  set val(prefix), file(reads) from chRawReadsFastx

  output:
  // set val(prefix), file("*_trimmed_G.R2.fastq.gz") into chTrimmedBc
  // set val(prefix), file("*_fastx.log") into chFastxLogs
  // file ("v_fastx.txt") into chFastxVersion
  set val(prefix), file("*_trimmed.R1.fastq"), file("*_trimmed.R2.fastq") into chTrimmedReads
  set val(prefix), file("*_trimmed.log") into chtrimmedReadsLog
  file("v_cutadapt.txt") into chCutadaptVersion

  script:
  linker_length = params.linker_length
  //  fastx_trimmer -i <(gzip -cd ${reads[1]}) -Q 33 -f ${linker_length} -z -o ${prefix}_trimmed_G.R2.fastq > ${prefix}_fastx.log
  // fastx_trimmer -h | grep "FASTX Toolkit" > v_fastx.txt
  """  
  cutadapt -q 33 -G CAACGTGATTGCTTGTGACTTAAA --minimum-length=15 --cores=${task.cpus} -o ${prefix}_trimmed.R1.fastq -p ${prefix}_trimmed.R2.fastq ${reads[0]} ${reads[1]} > ${prefix}_trimmed.log
     # cutadapt -u 84 --minimum-length=15 --cores=6 test.R2.fastq.gz -o test_trimmed.R2.fastq SAME 

  cutadapt -u 84 --minimum-length=15 --cores=6 test.R2.fastq.gz -o test_trimmed.R2.fastq
  """
}
*/

process fastxTrimmer {
  tag "${prefix}"
  label 'fastx'
  label 'medCpu'
  label 'medMem'

  publishDir "${params.outDir}/fastxTrimmer", mode: 'copy'

  input:
  set val(prefix), file(reads) from chRawReadsFastx

  output:
  set val(prefix), file("*_trimmed.R2.fastq.gz") into chTrimmedBc
  //set val(prefix), file("*_fastx.log") into chFastxLogs
  file ("v_fastx.txt") into chFastxVersion

  script:
  linker_length = params.linker_length
  """
  # Trim linker + barcode from R2 reads for genome aligning	
  fastx_trimmer -i <(gzip -cd ${reads[1]}) -f ${linker_length} -z -o ${prefix}_trimmed.R2.fastq.gz 

  fastx_trimmer -h | grep "FASTX Toolkit" > v_fastx.txt
  """
}

process readsAlignment {
  tag "${prefix}"
  label 'star'
  label 'extraCpu'
  label 'extraMem'

  publishDir "${params.outDir}/readsAlignment", mode: 'copy'

  input :
  file genomeIndex from chStar.collect()
  set val(prefix), file(trimmedR2) from chTrimmedBc
  set val(prefix), file(reads) from chAlignment
	
  output :
  set val(prefix), file("*Aligned.sortedByCoord.out.bam") into chAlignedBam
  file "*.out" into chAlignmentLogs
  file("v_star.txt") into chStarVersion

  script:
  """
  gzip -cd ${reads[0]} > R1.fastq
  gzip -cd ${trimmedR2} > R2.fastq
  # Align R2 reads on genome indexes - paired end with R1 - (STAR)
  # Run STAR on barcoded reads
   STAR --alignEndsType EndToEnd --outFilterMultimapScoreRange 2 --winAnchorMultimapNmax 1000 --alignIntronMax 1 --peOverlapNbasesMin 10 --alignMatesGapMax 450 --limitGenomeGenerateRAM 25000000000 --outSAMunmapped Within \
    --runThreadN ${task.cpus} \
    --genomeDir $genomeIndex \
    --readFilesIn R1.fastq R2.fastq \
    --runMode alignReads \
    --outFileNamePrefix ${prefix} \
    --outSAMtype BAM SortedByCoordinate 

  STAR --version &> v_star.txt

  rm *.fastq
  """
}

/*
process  addBarcodes {


  script:
  """
  samtools view -F 256 .bam | awk -v OFS='\t' 'NR%2==1{if(\$12==\"NH:i:1\"){mapped=1;print \$0} else{mapped=0;\$2=4;\$3=\"*\";\$4=0;\$6=\"*\";print \$0}} NR%2==0{if(mapped==1){print \$0} else{\$2=4;\$3=\"*\";\$4=0;\$6=\"*\";print \$0} }' > ${prefix}.sam





  """
}
*/
/***********
 * MultiQC *
 ***********/

/*process getSoftwareVersions{

  # Add the number of barcoded reads as first line of the index count file
  barcoded=`grep "Number of input reads" ${prefix}Log.final.out | cut -d'|' -f2 ` 
  echo "\$(echo 'Barcoded,'\$barcoded | cat - ${bcIndxCounts} )" > ${prefix}_index_mqc.log


  label 'python'
  label 'lowCpu'
  label 'lowMem'
  publishDir path: "${params.outDir}/software_versions", mode: "copy"

  when:
  !params.skipSoftVersions

  input:

  output:
  file 'software_versions_mqc.yaml' into softwareVersionsYaml

  script:
  """
  echo $workflow.manifest.version &> v_pipeline.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  scrape_software_versions.py &> software_versions_mqc.yaml
  """
}


process workflowSummaryMqc {
  when:
  !params.skipMultiQC

  output:
  file 'workflow_summary_mqc.yaml' into workflowSummaryYaml

  exec:
  def yaml_file = task.workDir.resolve('workflow_summary_mqc.yaml')
  yaml_file.text  = """
  id: 'summary'
  description: " - this information is collected when the pipeline is started."
  section_name: 'Workflow Summary'
  section_href: 'https://gitlab.curie.fr/data-analysis/chip-seq'
  plot_type: 'html'
  data: |
        <dl class=\"dl-horizontal\">
  ${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
  """.stripIndent()
}

process multiqc {
  label 'multiqc'
  label 'lowCpu'
  label 'lowMem'
  publishDir "${params.outDir}/MultiQC", mode: 'copy'

  when:
  !params.skipMultiQC

  input:
  file splan from chSplan.collect()
  file multiqcConfig from chMultiqcConfig
  file metadata from chMetadata.ifEmpty([])
  //file ('software_versions/*') from softwareVersionsYaml.collect().ifEmpty([])
  file ('workflow_summary/*') from workflowSummaryYaml.collect()

  output: 
  file splan
  file "*_report.html" into multiqc_report
  file "*_data"

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_report" : "--filename report"
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  //isPE = params.singleEnd ? "" : "-p"
  designOpts= params.design ? "-d ${params.design}" : ""
  modules_list = "-m custom_content"
  """
  mqc_header.py --splan ${splan} --name "PIPELINE" --version ${workflow.manifest.version} ${metadataOpts} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c multiqc-config-header.yaml -c $multiqcConfig $modules_list
  """
}

/****************
 * Sub-routines *
 ****************/

/*
process outputDocumentation {
  label 'python'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.outDir}/pipeline_info", mode: 'copy'

  input:
  file output_docs from chOutputDocs
  file images from chOutputDocsImages

  output:
  file "results_description.html"

  script:
  """
  markdown_to_html.py $output_docs -o results_description.html
  """
}

workflow.onComplete {

  // pipelineReport.html
  def reportFields = [:]
  reportFields['version'] = workflow.manifest.version
  reportFields['runName'] = customRunName ?: workflow.runName
  reportFields['success'] = workflow.success
  reportFields['dateComplete'] = workflow.complete
  reportFields['duration'] = workflow.duration
  reportFields['exitStatus'] = workflow.exitStatus
  reportFields['errorMessage'] = (workflow.errorMessage ?: 'None')
  reportFields['errorReport'] = (workflow.errorReport ?: 'None')
  reportFields['commandLine'] = workflow.commandLine
  reportFields['projectDir'] = workflow.projectDir
  reportFields['summary'] = summary
  reportFields['summary']['Date Started'] = workflow.start
  reportFields['summary']['Date Completed'] = workflow.complete
  reportFields['summary']['Pipeline script file path'] = workflow.scriptFile
  reportFields['summary']['Pipeline script hash ID'] = workflow.scriptId
  if(workflow.repository) reportFields['summary']['Pipeline repository Git URL'] = workflow.repository
  if(workflow.commitId) reportFields['summary']['Pipeline repository Git Commit'] = workflow.commitId
  if(workflow.revision) reportFields['summary']['Pipeline Git branch/tag'] = workflow.revision

  // Render the TXT template
  def engine = new groovy.text.GStringTemplateEngine()
  def tf = new File("$baseDir/assets/onCompleteTemplate.txt")
  def txtTemplate = engine.createTemplate(tf).make(reportFields)
  def reportTxt = txtTemplate.toString()

  // Render the HTML template
  def hf = new File("$baseDir/assets/onCompleteTemplate.html")
  def htmlTemplate = engine.createTemplate(hf).make(reportFields)
  def reportHtml = htmlTemplate.toString()

  // Write summary HTML to a file
  def outputSummaryDir = new File( "${params.summaryDir}/" )
  if( !outputSummaryDir.exists() ) {
    outputSummaryDir.mkdirs()
  }
  def outputHtmlFile = new File( outputSummaryDir, "pipelineReport.html" )
  outputHtmlFile.withWriter { w -> w << reportHtml }
  def outputTxtFile = new File( outputSummaryDir, "pipelineReport.txt" )
  outputTxtFile.withWriter { w -> w << reportTxt }

  // onComplete file
  File woc = new File("${params.outDir}/onComplete.txt")
  Map endSummary = [:]
  endSummary['Completed on'] = workflow.complete
  endSummary['Duration']     = workflow.duration
  endSummary['Success']      = workflow.success
  endSummary['exit status']  = workflow.exitStatus
  endSummary['Error report'] = workflow.errorReport ?: '-'
  String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")
  println endWfSummary
  String execInfo = "Execution summary\n${endWfSummary}\n"
  woc.write(execInfo)

  // final logs
  if(workflow.success){
      log.info "Pipeline Complete"
  }else{
    log.info "FAILED: $workflow.runName"
  }
}
*/