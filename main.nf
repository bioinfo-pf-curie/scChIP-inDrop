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
scChIP-inDrop
========================================================================================
 #### Homepage / Documentation
  scChIP-inDrop
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

  Skip options: All are false by default
    --skipSoftVersion [bool]      Do not report software version. Default is false.
    --skipMultiQC [bool]          Skips MultiQC. Default is false.
    --skipBigWig [bool]           Skips BigWig generation. Default is false.

  Other options:
    --outDir [file]               The output directory where the results will be saved
    -name [str]                   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    --keepRTdup [bool]             Keep RT duplicats. Default is false.
    --distDup [int]                Select the number of bases after gene start sites to detect duplicates. Default is 50.
    --minCounts [int]              Select the minimum count per barcodes after removing duplicates. Default is 500.
    --keepBlacklistRegion [bool]     Keep black region. Default is false.
    --binSize [int]               Bin size to use (in base pairs). Default is 50000.
    --tssWindow [int]              TSS window size to use (in base pairs). Default is 5000.
 
  =======================================================

  Available Profiles:
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
    .set { chGtftoTSSBed }
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

params.blackList = genomeRef ? params.genomes[ genomeRef ].blackList ?: false : false
if (params.blackList) {
  Channel  
    .fromPath(params.blackList)
    .ifEmpty { exit 1, "blackList annotation file not found: ${params.bed12}" }
    .into { chFilterBlackReg; chFilterBlackReg_bamToBigWig } 
}else {
  exit 1, "BED12 annotation file not not found: ${params.bed12}"
}

//------- Custom barcode indexes--------
//--------------------------------------
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
      .into { chRawReadsBowtie2; chRawReadsFastx; chAlignment; chRawReadsWhiteList  }
  }else{
    Channel
      .from(file("${params.samplePlan}"))
      .splitCsv(header: false)
      .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
      .into { chRawReadsBowtie2; chRawReadsFastx; chAlignment; chRawReadsWhiteList }
   }
  params.reads=false
}
else if(params.readPaths){
  if(params.singleEnd){
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied." }
      .into { chRawReadsBowtie2; chRawReadsFastx ; chAlignment; chRawReadsWhiteList }
  } else {
    Channel
      .from(params.readPaths)
      .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
      .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied." }
      .into { chRawReadsBowtie2 ; chRawReadsFastx; chAlignment; chRawReadsWhiteList }
  }
} else {
  Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { chRawReadsBowtie2 ; chRawReadsFastx; chAlignment; chRawReadsWhiteList }
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

/*
process bcMapping {
  tag "${prefix}"
  label 'bowtie2'
  label 'highCpu'
  label 'highMem'

  publishDir "${params.outDir}/bcMapping", mode: 'copy'

  input:
  set val(prefix), file(reads) from chRawReadsBowtie2
 
  output:
  set val(prefix), file("*_read_barcodes.txt") into chReadBcNames
  file ("v_bowtie2.txt") into chBowtie2Version
  // summary of counts
  set val(prefix) , file("*_bowtie2.log") into chBowtie2Log
  // for bowtie2 module in mqc
  set val(prefix) , file("*_index_1_bw2.log") into chIndex1Bowtie2Log
  set val(prefix) , file("*_index_2_bw2.log") into chIndex2Bowtie2Log
  set val(prefix) , file("*_index_3_bw2.log") into chIndex3Bowtie2Log

  script:
  """
  ##Extract three indexes from reads : 1 - 16 = index 1 ; 21 - 36 = index 2; 41 - 56 = index 3
  gzip -cd  ${reads[1]} | awk -v start_index_1=5 -v size_index=16  'NR%4==1{print ">"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_1,size_index)}' > read_indexes_1.fasta
  
  gzip -cd  ${reads[1]} | awk -v start_index_2=25 -v size_index=16 'NR%4==1{print  \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_2,size_index)}' > read_indexes_2.fasta
  
  gzip -cd  ${reads[1]} | awk -v start_index_3=45 -v size_index=16 'NR%4==1{print \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_3,size_index)}' > read_indexes_3.fasta
  
  #Map INDEXES 1 against Index1 library
  bowtie2 -x /data/users/lhadjabe/Gitlab/ChIP-seq_single-cell_LBC/Barcodes/LBC/bowtie_2_index_short/ref_index_1 -f read_indexes_1.fasta \
          -N 1 -L 8 --rdg 0,7 --rfg 0,7 --mp 7,7 --ignore-quals --score-min L,0,-1 -t --no-unal --no-hd -p ${task.cpus} > index_1_bowtie2.sam 2> ${prefix}_index_1_bw2.log
  #Keep only reads that were matched by a unique index 1 + counting matched index1
  awk '/XS/{next} \$2!=4{print \$1,\$3;count++} ;END{print count > \"count_index_1\"}' index_1_bowtie2.sam > reads_matching_index_1.txt
  

  #Map INDEXES 2 against Index2 library
  bowtie2 -x /data/users/lhadjabe/Gitlab/ChIP-seq_single-cell_LBC/Barcodes/LBC/bowtie_2_index_short/ref_index_2 -f read_indexes_2.fasta -N 1 -L 8 --rdg 0,7 --rfg 0,7 --mp 7,7 --ignore-quals --score-min L,0,-1 -t --no-unal --no-hd -p ${task.cpus} > index_2_bowtie2.sam 2> ${prefix}_index_2_bw2.log
  #Keep only reads that were matched by a unique index 2 + counting matched index2
  awk '/XS/{next} \$2!=4{print \$1,\$3;count++} ;END{print count > \"count_index_2\"}' index_2_bowtie2.sam > reads_matching_index_2.txt
  
  
  #Map INDEXES 3 against Index3 library
  bowtie2 -x /data/users/lhadjabe/Gitlab/ChIP-seq_single-cell_LBC/Barcodes/LBC/bowtie_2_index_short/ref_index_3 -f read_indexes_3.fasta -N 1 -L 8 --rdg 0,7 --rfg 0,7 --mp 7,7 --ignore-quals --score-min L,0,-1 -t --no-unal --no-hd -p ${task.cpus} > index_3_bowtie2.sam 2> ${prefix}_index_3_bw2.log
  #Keep only reads that were matched by a unique index 3 + counting matched index3
  awk '/XS/{next} \$2!=4{print \$1,\$3;count++} ;END{print count > \"count_index_3\"}' index_3_bowtie2.sam > reads_matching_index_3.txt
  
  ##Sort indexes by read name: 
  sort -T /scratch/ --parallel=${task.cpus} -k1,1 reads_matching_index_1.txt > reads_matching_index_1_sorted.txt
  
  rm reads_matching_index_1.txt
  
  sort -T /scratch/ --parallel=${task.cpus} -k1,1 reads_matching_index_2.txt > reads_matching_index_2_sorted.txt

  rm reads_matching_index_2.txt
  
  sort -T /scratch/ --parallel=${task.cpus} -k1,1 reads_matching_index_3.txt > reads_matching_index_3_sorted.txt
  
  rm reads_matching_index_3.txt
  
  #Join indexes 1 & 2 together (inner join)
  join -t\$' ' -1 1 -2 1 reads_matching_index_1_sorted.txt reads_matching_index_2_sorted.txt > tmp
  
  #Count matched index 1 & 2
  echo \$(wc -l tmp) | cut -d' ' -f1 > count_index_1_2
  	
  
  #Join indexes (1 & 2) & 3 together to recompose full barcode (inner join)
  join -t\$' ' -1 1 -2 1 tmp reads_matching_index_3_sorted.txt > final
  
  
  #Reformat & count matched index (1 & 2 & 3) <=> barcode
  awk '{print substr(\$1,1)\"\tBC\"substr(\$2,2)substr(\$3,2)substr(\$4,2);count++} ;END{print count > \"count_index_1_2_3\"}' final > ${prefix}_read_barcodes.txt
  
  ##Write logs
  n_index_1=\$(cat count_index_1)
  n_index_2=\$(cat count_index_2)
  n_index_3=\$(cat count_index_3)
  n_index_1_2=\$(cat count_index_1_2)
  n_index_1_2_3=\$(cat count_index_1_2_3)

  ## logs
  echo "## Number of matched indexes 1: \$n_index_1" > ${prefix}_bowtie2.log
  echo "## Number of matched indexes 2: \$n_index_2" >> ${prefix}_bowtie2.log
  echo "## Number of matched indexes 1 and 2: \$n_index_1_2" >> ${prefix}_bowtie2.log
  echo "## Number of matched indexes 3: \$n_index_3" >> ${prefix}_bowtie2.log
  echo "## Number of matched barcodes: \$n_index_1_2_3" >> ${prefix}_bowtie2.log

  ## version
  bowtie2 --version > v_bowtie2.txt
  """
}
*/

process bcMapping {
  tag "${prefix} - ${index}"
  label 'bowtie2'
  label 'extraCpu'
  label 'extraMem'

  publishDir "${params.outDir}/bcMapping", mode: 'copy'

  input:
  set val(prefix), file(reads), val(index), file(bwt2Idx) from chRawReadsBowtie2.combine(chIndexBwt2)
 
  output:
  // read IDs matching each index
  set val(prefix), file("*${index}_ReadsMatchingSorted.txt") into chReadsMatchingIndex
  // counts of the number of reads matching each index
  set val(prefix), file("*${index}_count_index.txt") into chIndexCount
  // for bowtie2 module in mqc
  file ("*Bowtie2.log") into chIndexBowtie2Log 
  // version
  file ("v_bowtie2.txt") into chBowtie2Version

  script:
  // !!! Old design pas ajouté en option => to do 
  start = params.barcodes[ index ].start
  size = params.barcodes[ index ].size
  base = params.barcodes[ index ].base
  oprefix = "${prefix}_${index}"
  """
  ##Extract three indexes from reads (old design): 1 - 16 = index 1 ; 21 - 36 = index 2; 41 - 56 = index 3
  gzip -cd  ${reads[1]} | awk -v start_index_1=${start} -v size_index=${size} 'NR%4==1{print ">"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_1,size_index)}' > ${oprefix}Reads.fasta

  #Map indexes (-f) against Index libraries (-x)
  bowtie2 \
    -x ${bwt2Idx}/${base} \
    -f ${oprefix}Reads.fasta \
    -N 1 -L 8 --rdg 0,7 --rfg 0,7 --mp 7,7 \
    --ignore-quals --score-min L,0,-1 -t \
    --no-unal --no-hd \
    -p ${task.cpus} > ${oprefix}Bowtie2.sam 2> ${oprefix}Bowtie2.log
  #Keep only reads that were matched by a unique index 1 + counting matched index1
  awk '/XS/{next} \$2!=4{print \$1,\$3}' ${oprefix}Bowtie2.sam > ${oprefix}ReadsMatching.txt 
  wc -l < ${oprefix}ReadsMatching.txt  > ${oprefix}_count_index.txt
  
  ##Sort indexes by read name: 
  sort -T /scratch/ --parallel=${task.cpus} -k1,1 ${oprefix}ReadsMatching.txt > ${oprefix}_ReadsMatchingSorted.txt 

  #delete useless files
  rm ${oprefix}ReadsMatching.txt ${oprefix}Bowtie2.sam ${oprefix}Reads.fasta
  
  ## version
  bowtie2 --version > v_bowtie2.txt
  """
}

process bcSubset {
  tag "${prefix}"
  label 'unix'
  label 'extraCpu'
  label 'extraMem'
  publishDir "${params.outDir}/bcSubset", mode: 'copy'

  input:
  // read IDs matching each index
  set val(prefix), file(readsMatchingSorted) from chReadsMatchingIndex.groupTuple()
  // counts of the number of reads matching each index
  set val(prefix), file(count_index) from chIndexCount.groupTuple()
  
  output:
  // correctly barcoded reads
  set val(prefix), file("*_read_barcodes.txt") into chReadBcNames
  // summary of counts
  set val(prefix) , file("*_bowtie2.log") into chBowtie2Log

  script:
  """  
  #Join indexes 1 & 2 together (inner join)
  join -t\$' ' -1 1 -2 1 ${prefix}_indexB_ReadsMatchingSorted.txt ${prefix}_indexC_ReadsMatchingSorted.txt > tmp
  
  #Count matched index 1 & 2
  echo \$(wc -l tmp) | cut -d' ' -f1 > count_index_1_2
  
  #Join indexes (1 & 2) & 3 together to recompose full barcode (inner join)
  join -t\$' ' -1 1 -2 1 tmp ${prefix}_indexD_ReadsMatchingSorted.txt > final
  
  #Reformat & count matched index (1 & 2 & 3) <=> barcode
  awk '{print substr(\$1,1)\"\tBC\"substr(\$2,2)substr(\$3,2)substr(\$4,2);count++} ;END{print count > \"count_index_1_2_3\"}' final > ${prefix}_read_barcodes.txt
  
  ##Write logs
  n_index_1=\$(cat ${prefix}_indexB_count_index.txt)
  n_index_2=\$(cat ${prefix}_indexC_count_index.txt)
  n_index_3=\$(cat ${prefix}_indexD_count_index.txt)
  n_index_1_2=\$(cat count_index_1_2)
  n_index_1_2_3=\$(cat count_index_1_2_3)

  ## logs
  echo "## Number of matched indexes 1: \$n_index_1" > ${prefix}_bowtie2.log
  echo "## Number of matched indexes 2: \$n_index_2" >> ${prefix}_bowtie2.log
  echo "## Number of matched indexes 1 and 2: \$n_index_1_2" >> ${prefix}_bowtie2.log
  echo "## Number of matched indexes 3: \$n_index_3" >> ${prefix}_bowtie2.log
  echo "## Number of matched barcodes: \$n_index_1_2_3" >> ${prefix}_bowtie2.log

  rm count_index_1_2 count_index_1_2_3 tmp final
  """
}


//2-  Trim R2 reads for genome aligning	
process trimReads {
  tag "${prefix}"
  label 'cutadapt'
  label 'medCpu'
  label 'medMem'

  publishDir "${params.outDir}/trimR2", mode: 'copy'

  input:
  //set val(prefix), file(barcodedR1), file(barcodedR2) from chBarcodedReads_Fastx
  set val(prefix), file(reads) from chRawReadsFastx

  output:
  set val(prefix), file("*_trimmed.R2.fastq") into chTrimmedReads
  set val(prefix), file("*_trimmed.log") into chTrimmedReadsLog
  file("v_cutadapt.txt") into chCutadaptVersion

  script:
  barcode_linker_len = params.barcode_linker_length
  """
  # Trim linker + barcode from R2 reads for genome aligning	
  cutadapt -u ${barcode_linker_len} --cores=${task.cpus} ${reads[1]} -o ${prefix}_trimmed.R2.fastq > ${prefix}_trimmed.log
  cutadapt --version &> v_cutadapt.txt
  """
}

//3- Align R2 reads on genome indexes - paired end with R1 - (STAR)
process readsAlignment {
  tag "${prefix}"
  label 'star'
  label 'extraCpu'
  label 'extraMem'

  publishDir "${params.outDir}/readsAlignment", mode: 'copy'

  input :
  file genomeIndex from chStar.collect()
  set val(prefix), file(reads) from chAlignment
  set val(prefix), file(trimmedR2) from chTrimmedReads
  //set val(prefix), file(barcodedR1), file(barcodedR2) from chBarcodedReads_Star

  output :
  set val(prefix), file("*_aligned.bam") into chAlignedBam
  file "*.out" into chAlignmentLogs
  file("v_star.txt") into chStarVersion
  file("v_samtools.txt") into chSamtoolsVersion

  script:
  """
  gzip -cd ${reads[0]} > R1.fastq
  # Align R2 reads on genome indexes - paired end with R1 - (STAR)
  # Run STAR on barcoded reads
   STAR --alignEndsType EndToEnd --outFilterMultimapScoreRange 2 --winAnchorMultimapNmax 1000 \
    --alignIntronMax 1 --peOverlapNbasesMin 10 --alignMatesGapMax 450 \
    --limitGenomeGenerateRAM 25000000000 --outSAMunmapped Within \
    --runThreadN ${task.cpus} \
    --genomeDir $genomeIndex \
    --readFilesIn R1.fastq ${trimmedR2} \
    --runMode alignReads \
    --outFileNamePrefix ${prefix} 

    # Seeds = exactly mapped portions of a read 
    # Anchors = genomic windows (loci) where a read mapping is attempted 
    # --winAnchorMultimapNmax 1000 =  reads having up to 1000 anchor seeds

  STAR --version &> v_star.txt

  samtools view -@ ${task.cpus} -bS ${prefix}Aligned.out.sam > ${prefix}_aligned.bam
  samtools sort -n -@ ${task.cpus} ${prefix}_aligned.bam -o ${prefix}_nsorted.bam && mv ${prefix}_nsorted.bam ${prefix}_aligned.bam
  rm ${prefix}Aligned.out.sam
  samtools --version > v_samtools.txt
  """
}

//4- Add cellular Barcode - (STAR)
process  addBarcodes {
  tag "${prefix}"
  label 'samtools'
  label 'highCpu'
  label 'highMem'

  publishDir "${params.outDir}/addBarcodes", mode: 'copy'

  input:
  set val(prefix), file(alignedBam) from chAlignedBam
  set val(prefix), file(readBarcodes) from chReadBcNames

  output:
  set (prefix), file("*_flagged.bam") into chAddedBarcodes

  script:
  """
  # Remove secondary aligned reads (256 <=> "not primary alignment") & If R1 is unmapped or multimapped (NH != 1), tag R1 & R2 with flag "4" <=> "unmapped" & "chr" = '*'
  samtools view -F 256 ${alignedBam} | awk -v OFS='\t' 'NR%2==1{if(\$12==\"NH:i:1\"){mapped=1;print \$0} else{mapped=0;\$2=4;\$3=\"*\";\$4=0;\$6=\"*\";print \$0}} NR%2==0{if(mapped==1){print \$0} else{\$2=4;\$3=\"*\";\$4=0;\$6=\"*\";print \$0} }' > ${prefix}.sam

  # If read is mapped R1 & unmapped R2 -> set R2 position as '2147483647'
  cat ${prefix}.sam | awk -v OFS='\t' 'NR%2==1{print \$0} NR%2==0{if(\$3==\"*\"){\$4=2147483647;print \$0} else{print \$0} }' > ${prefix}_2.sam

  # Remove comments from the header that produce bugs in the count phase
  samtools view -H ${alignedBam} | sed '/^@CO/ d' > ${prefix}_header.sam

  cat ${prefix}_2.sam >> ${prefix}_header.sam && mv ${prefix}_header.sam ${prefix}.sam && samtools view -b -@ ${task.cpus} ${prefix}.sam > ${prefix}_unique.bam

  rm -f ${prefix}_2.sam ${prefix}.sam
  
  #Keeping R1 aligned + R2 start as tag 'XS' (Switch from Paired End Bam to Single End Bam)
  samtools view ${prefix}_unique.bam | awk '{OFS = \"\t\" ; if(NR%2==1 && !(\$3==\"*\")) {R1=\$0} else if(NR%2==1){R1=0}; if(NR%2==0 && !(R1==0)){tagR2Seq=\"XD:Z:\"\$10; tagR2Pos=\"XS:i:\"\$4;print R1,tagR2Pos,tagR2Seq}}' > ${prefix}_unique.sam
  
  #Sort and join on read names reads barcoded and reads mapped to genome (barcode as tag 'XB') --> filter out unbarcoded OR unmapped reads
  sort -T /scratch/ --parallel=${task.cpus} -k1,1 ${prefix}_unique.sam > ${prefix}_unique_sorted.sam
  
  join -1 1  -2 1  ${prefix}_unique_sorted.sam <(awk -v OFS=\"\t\" '{print \$1,\"XB:Z:\"\$2}' ${readBarcodes}) > ${prefix}_flagged.sam
  
  sed -i 's/ /\t/g' ${prefix}_flagged.sam
  
  #Remove comments from the header that produce bugs in the count phase
  samtools view -H ${prefix}_unique.bam | sed '/^@CO/ d' > ${prefix}_header.sam
  
  cat ${prefix}_flagged.sam >> ${prefix}_header.sam && mv ${prefix}_header.sam ${prefix}_flagged.sam && samtools view -@ ${task.cpus} -b ${prefix}_flagged.sam > ${prefix}_flagged.bam
  
  #Cleaning
  rm -f ${prefix}_unique.bam ${prefix}_flagged.sam ${prefix}_unique_sorted.sam
  """
}


//5- Remove PCR and Reverse Transcription duplicate  - (STAR)
process  removePcrRtDup {
  tag "${prefix}"
  label 'samtools'
  label 'highCpu'
  label 'highMem'

  publishDir "${params.outDir}/removePcrRtDup", mode: 'copy'

  input:
  set (prefix), file(flaggedBam) from chAddedBarcodes

  output:
  set (prefix), file("*_flagged_rmPCR_RT.bam") into chRTremoved
  set (prefix), file("*_removePcrRtDup.log") into chPcrRtCountsLog
  //set (prefix), file("*_flagged_rmPCR_RT.count") into chPcrRtCountTable

  script:
  """
  ##Sort by barcode then chromosome then position R2
  #Find the column containing the barcode tag XB
  barcode_field=\$(samtools view ${flaggedBam} | sed -n \"1 s/XB.*//p\" | sed 's/[^\t]//g' | wc -c)
  #Find the column containing the position R2 tag XS
  posR2_field=\$(samtools view ${flaggedBam} | sed -n \"1 s/XS.*//p\" | sed 's/[^\t]//g' | wc -c)

  printf '@HD\tVN:1.4\tSO:unsorted\n' > ${prefix}_header.sam
  samtools view -H ${prefix}_flagged.bam | sed '/^@HD/ d' >> ${prefix}_header.sam

  #Sort by barcode then chromosome then position R2 then Position R1 (for the PCR removal) 
  #It is important to sort by R1 pos also	because the removal is done by comparing consecutive lin
  samtools view ${prefix}_flagged.bam | LC_ALL=C sort -T /scratch/ --parallel=${task.cpus} -t \$'\t' -k \"\$barcode_field.8,\$barcode_field\"n -k 3.4,3g -k \"\$posR2_field.6,\$posR2_field\"n -k 4,4n >> ${prefix}_header.sam && samtools view -@ ${task.cpus} -b ${prefix}_header.sam > ${prefix}_flagged.sorted.bam

  samtools view ${prefix}_flagged.sorted.bam | awk -v bc_field=\$barcode_field '{print substr(\$bc_field,6)}' |  uniq -c > ${prefix}_flagged.count
  samtools view ${prefix}_flagged.sorted.bam | awk -v bc_field=\$barcode_field -v R2_field=\$posR2_field 'BEGIN {countR1unmappedR2=0;countPCR=0};NR==1{print \$0;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\"); lastR1Pos=\$4} ; NR>=2{split( \$R2_field,R2Pos,\":\");R1Pos=\$4; if(R2Pos[3]==2147483647){print \$0;countR1unmappedR2++; next}; if( (R1Pos==lastR1Pos) && (R2Pos[3]==lastR2Pos[3]) && ( \$3==lastChrom ) && (\$bc_field==lastBarcode) ){countPCR++;next} {print \$0;lastR1Pos=\$4;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\") }} END {print countPCR > \"count_PCR_duplicates\";print countR1unmappedR2 > \"countR1unmappedR2\"}' > ${prefix}_flagged_rmPCR.sam
  
  samtools view -H ${prefix}_flagged.sorted.bam  | sed '/^@CO/ d' > ${prefix}_header.sam
  
  cat ${prefix}_flagged_rmPCR.sam >> ${prefix}_header.sam && samtools view -@ ${task.cpus} -b ${prefix}_header.sam > ${prefix}_flagged_rmPCR.bam
  
  #Create count Table from flagged - PCR dups (already sorted by barcode)
  samtools view ${prefix}_flagged_rmPCR.bam | awk -v bc_field=\$barcode_field '{print substr(\$bc_field,6)}' |  uniq -c > ${prefix}_flagged_rmPCR.count
  
  ## Sort flagged_rmPCR file
  samtools sort -@ ${task.cpus} ${prefix}_flagged_rmPCR.bam > ${prefix}_flagged_rmPCR_sorted.bam
  
  ## Rename flagged_rmPCR file
  mv ${prefix}_flagged_rmPCR_sorted.bam ${prefix}_flagged_rmPCR.bam
  
  ## Index flagged_rmPCR file
  samtools index ${prefix}_flagged_rmPCR.bam
  
  if [ ${params.keepRTdup} == 'false' ] 
  then

    #Remove RT duplicates (if two consecutive reads have the same barcode and same R2 chr&start) but not same R1 
    cat ${prefix}_flagged_rmPCR.sam | awk -v bc_field=\$barcode_field -v R2_field=\$posR2_field 'BEGIN{count=0};NR==1{print \$0;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\")} ; NR>=2{split( \$R2_field,R2Pos,\":\");if((R2Pos[3]==lastR2Pos[3]) && (R2Pos[3]!=2147483647) && (lastR2Pos[3]!=2147483647)  && ( \$3==lastChrom ) && (\$bc_field==lastBarcode) ){count++;next} {print \$0;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\") }} END {print count > \"count_RT_duplicates\"}' > ${prefix}_flagged_rmPCR_RT.sam
    
    samtools view -H ${prefix}_flagged.bam  | sed '/^@CO/ d' > ${prefix}_header.sam

    cat ${prefix}_flagged_rmPCR_RT.sam >> ${prefix}_header.sam && samtools view -@ ${task.cpus} -b ${prefix}_header.sam > ${prefix}_flagged_rmPCR_RT.bam 
    
    #Create count Table from flagged - PCR dups - RT dups  (already sorted by barcode)
    samtools view ${prefix}_flagged_rmPCR_RT.bam | awk -v bc_field=\$barcode_field '{print substr(\$bc_field,6)}' | uniq -c > ${prefix}_flagged_rmPCR_RT.count
    
    ## Sort flagged_rmPCR_RT file
    samtools sort -@ ${task.cpus} ${prefix}_flagged_rmPCR_RT.bam > ${prefix}_flagged_rmPCR_RT_sorted.bam
    
    ## Rename flagged_rmPCR_RT file
    mv ${prefix}_flagged_rmPCR_RT_sorted.bam ${prefix}_flagged_rmPCR_RT.bam
  
  else

    ## Copy flagged_rmPCR to flagged_rmPCR_RT
    cp ${prefix}_flagged_rmPCR.bam ${prefix}_flagged_rmPCR_RT.bam
    
    cp ${prefix}_flagged_rmPCR.count ${prefix}_flagged_rmPCR_RT.count
    
    ## Set RT duplicate count to 0
    echo 0 > count_RT_duplicates
  fi
  
  ## Index flagged_rmPCR_RT file
  samtools index ${prefix}_flagged_rmPCR_RT.bam
  
  ## Rename flagged.sorted -> flagged
  mv ${prefix}_flagged.sorted.bam ${prefix}_flagged.bam
  
  ## Logs
  n_mapped_barcoded=\$(samtools view -c  ${prefix}_flagged.bam)
  n_pcr_duplicates=\$(cat count_PCR_duplicates)
  n_rt_duplicates=\$(cat count_RT_duplicates)
  n_R1_mapped_R2_unmapped=\$(cat countR1unmappedR2)
  
  n_unique_except_R1_unmapped_R2=\$((\$n_mapped_barcoded - \$n_pcr_duplicates - \$n_rt_duplicates))

  echo "## Number of reads mapped and barcoded: \$n_mapped_barcoded" >> ${prefix}_removePcrRtDup.log
  echo "## Number of pcr duplicates: \$n_pcr_duplicates" >> ${prefix}_removePcrRtDup.log
  echo "## Number of rt duplicates: \$n_rt_duplicates" >> ${prefix}_removePcrRtDup.log
  echo "## Number of R1 mapped but R2 unmapped: \$n_R1_mapped_R2_unmapped" >> ${prefix}_removePcrRtDup.log
  echo "## Number of reads after PCR and RT removal (not R1 unmapped R2): \$n_unique_except_R1_unmapped_R2" >> ${prefix}_removePcrRtDup.log

  ## Remove all non used files
  rm -f count* *.sam
  """
}

// 6-Remove duplicates by window (if R2 is unmapped)
process  removeDup {
  tag "${prefix}"
  label 'bedtools'
  label 'medCpu'
  label 'medMem'

  publishDir "${params.outDir}/removeDup", mode: 'copy'

  input:
  set (prefix), file(noPcrRtBam) from chRTremoved
  file blackListBed from chFilterBlackReg  

  output:
  set (prefix), file("*_rmDup.bam"),  file("*_rmDup.bam.bai") into chNoDup_ScBed, chNoDup_bigWig, chNoDup_countMatrices
  set (prefix), file("*_rmDup.count") into chDupCounts, chRemoveDupBarcodeLog, chDistribUMIs
  set (prefix), file("*_rmDup.log") into chRemoveDupLog
  file("v_bedtools.txt") into chBedtoolsVersion
  
  script:
  """
  # window param
  if [ ! -z ${params.distDup} ]; then
	  rmDup.py -v -i ${noPcrRtBam} -o ${prefix}_rmDup.bam -d ${params.distDup} > ${prefix}_rmDup.log
  else
	  rmDup.py -v -i ${noPcrRtBam} -o ${prefix}_rmDup.bam > ${prefix}_rmDup.log
  fi
    
  #Create count Table from flagged - PCR dups - RT dups and window-based rmDup (need to sort by barcode)
  barcode_field=\$(samtools view ${prefix}_rmDup.bam  | sed -n \"1 s/XB.*//p\" | sed 's/[^\t]//g' | wc -c)
  
  samtools view ${prefix}_rmDup.bam | awk -v bc_field=\$barcode_field '{print substr(\$bc_field,6)}' | sort | uniq -c > ${prefix}_rmDup.count	

  # Removing encode black regions
  if [[${params.keepBlacklistRegion} == "false"]]
  then
    bedtools intersect -v -abam ${prefix}_rmDup.bam -b ${blackListBed} > ${prefix}_rmDup_rmBlackReg.bam && mv ${prefix}_rmDup_rmBlackReg.bam ${prefix}_rmDup.bam
  fi

  bedtools --version &> v_bedtools.txt

  ## Index BAM file
  samtools index ${prefix}_rmDup.bam
  """
}

// Plot read distributions accross cells (used to filtre out low read cells)
process distribUMIs{
  tag "${prefix}"
  label 'R'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/distribUMIs", mode: 'copy'

  input:
  
  set val(prefix), file(countedReadsPerCell_matrix) from chDistribUMIs

  output:
  set val(prefix), file("*distDF.mqc") into mqcDistribUMI
  set val(prefix), file("*distribution.pdf") into pdfDist
  file ("v_R.txt") into chRversion

  script:
  """
  umisDistribution.r ${countedReadsPerCell_matrix} ${prefix}
  R --version &> v_R.txt
  """ 
}

// 7-Generate bigWig fileq - Create pseudo bulk
process bamToBigWig{
  tag "${prefix}"
  label 'deeptools'
  label 'extraCpu'
  label 'extraMem'

  publishDir "${params.outDir}/bamToBigWig", mode: 'copy'

  when:
  !params.skipBigWig

  input:
  set val(prefix), file (rmDupBam), file (rmDupBai) from chNoDup_bigWig
  file blackListBed from chFilterBlackReg_bamToBigWig    
      
  output:
  set val(prefix), file("*.bw") into chBigWig
  file("v_deeptools.txt") into chBamCoverageVersion
  set (prefix), file("*_bamToBigWig.log") into chBamToBigLog
  
  script:
  """
  if [[${params.keepBlacklistRegion} == "false"]]
  then
      bamCoverage --bam ${rmDupBam} --outFileName ${prefix}.bw --numberOfProcessors ${task.cpus} --normalizeUsing CPM --ignoreForNormalization chrX --binSize 50 --smoothLength 500 --extendReads 150 --blackListFileName $blackListBed &> ${prefix}_bamToBigWig.log
  else
      bamCoverage --bam ${rmDupBam} --outFileName ${prefix}.bw --numberOfProcessors ${task.cpus} --normalizeUsing CPM --ignoreForNormalization chrX --binSize 50 --smoothLength 500 --extendReads 150 &> ${prefix}_bamToBigWig.log
  fi

  deeptools --version &> v_deeptools.txt
  """
}

//7bis-Generate scBed file
process bamToScBed{
  tag "${prefix}"
  label 'samtools'
  label 'highCpu'
  label 'highMem'

  publishDir "${params.outDir}/bamToScBed", mode: 'copy'

  input:
  set (prefix), file (rmDupBam), file (rmDupBai) from chNoDup_ScBed
  
  output:
  file (prefix) into chScBed
  
  script:
  """
  mkdir -p ${prefix}
  for i in \$(echo ${params.minCounts} | sed 's/,/ /g'); do mkdir -p ${prefix}/scBed_\$i/ ; done
  
  #Get barcode field & read length
  barcode_field=\$(samtools view ${rmDupBam}  | sed -n "1 s/XB.*//p" | sed 's/[^\t]//g' | wc -c)
  
  #Create header
  samtools view -H ${rmDupBam} | sed '/^@HD/ d' > ${prefix}_tmp_header.sam
  
  #Sort by Barcode, Chr, Pos R1 :
  samtools view ${rmDupBam} | LC_ALL=C sort -T /scratch/ --parallel=${task.cpus} -t \$'\t' -k \"\$barcode_field.8,\$barcode_field\"n -k 3.4,3g -k 4,4n >> ${prefix}_tmp_header.sam
  
  samtools view -@ ${task.cpus} -b ${prefix}_tmp_header.sam > ${prefix}_tmp.sorted.bam
  
  #Convert to bedgraph: Input must be sorted by barcode, chr, position R1
  samtools view ${prefix}_tmp.sorted.bam | awk -v odir=${prefix}/scBed -v bc_field=\$barcode_field -v OFS="\t" -v count=${params.minCounts} '
  BEGIN{
    split(count,min_counts,",")
  }
  NR==1{
    lastBC=substr(\$bc_field,6,15);
    i=1
    chr[i] = ${prefix}
    start[i] = ${params.minCounts}
    end[i] = ${params.minCounts} +1
  }
  NR>1{
  if(lastBC==substr(\$bc_field,6,15)){
    i = i +1
    chr[i] = ${prefix}
    start[i] = ${params.minCounts}
    end[i] = ${params.minCounts} +1
    }
    else{
    for(m=1; m<=length(min_counts);m++){
      if(i > min_counts[m]){
        for (x=1; x<=i; x++){
          out = odir"_"min_counts[m]"/"lastBC".bed"
          print chr[x],start[x],end[x] >> out
        }
      }
    }
    i=0
    }
     lastBC=substr(\$bc_field,6,15);
}
'

#Gzip
if [ -f scBed*/*.bed ];then
  for i in scBed*/*.bed; do gzip -9 \$i; done
fi

rm -f ${prefix}_tmp_header.sam ${prefix}_tmp.sorted.bam
"""
}

// 8 - Generate TSS bed file
process gtfToTSSBed {
  tag "${prefix}"
  label 'unix'
  label 'medCpu'
  label 'medMem'

  publishDir "${params.outDir}/gtfToTSSBed", mode: 'copy'

  input:
  file gtf from chGtftoTSSBed
  
  output:
  file("*_TSS*.bed") into tssBedFile

  script:
  """ 
  create_transcript_annotation.sh $gtf ${params.tssWindow}
  """
}


// 8 - Generate count matrix
process countMatrices {
  tag "${prefix}"
  label 'samtools'
  label 'highCpu'
  label 'highMem'

  publishDir "${params.outDir}/countMatrices", mode: 'copy'

  input:
  file tssBed from tssBedFile
  set (prefix), file (rmDupBam), file (rmDupBai) from chNoDup_countMatrices
  set (prefix), file(countTable) from chDupCounts

  output:
  set val(prefix), file ("*.tsv.gz") into chCountMatrices
  set val(prefix), file ("*_counts.log") into chCountMatricesLog
  file("v_python.txt") into chPythonVersion
  
  script:
  """
  #Counting unique BCs
  bc_prefix=\$(basename ${rmDupBam} | sed -e 's/.bam\$//')
  barcodes=\$(wc -l ${countTable} | awk '{print \$1}')
  echo "Barcodes found = \$barcodes" > ${prefix}_counts.log
  
  # Counts are generated either per bin (--bin) and per genomics features (--bed)

  # 50000 & 5000
  for bsize in ${params.binSize}
  do
    opts="-b \$bsize "
    if [ ! -z ${params.minCounts} ]; then
        opts="\$opts -f ${params.minCounts} "
	  fi
        sc2counts.py -i ${rmDupBam} -o ${prefix}_counts_bin_${params.binSize}.tsv \$opts -s \$barcodes -v
  done
 
  # 5000
  for bed in $tssBed
  do
    opts="-B \$bed"
    if [ ! -z ${params.minCounts} ]; then
        opts="\$opts -f ${params.minCounts} "
    fi
    sc2counts.py -i ${rmDupBam} -o ${prefix}_counts_TSS_${params.tssWindow}.tsv \$opts -s \$barcodes -v
  done
  
  for i in ${prefix}*.tsv; do gzip -9 \$i; done
  
  python --version &> v_python.txt
  """
}

process create10Xoutput{
  tag "${prefix}"
  label 'R'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/create10Xoutput", mode: 'copy'

  input:
  set val(prefix), file(binMatx), file(tssMatx) from chCountMatrices

  output:
  file (prefix) into chOut10X
  
  script:
  """
  mkdir ${prefix}

  create10Xoutput.r ${binMatx} binMatx/
  mv binMatx/ ${prefix}/

  create10Xoutput.r ${tssMatx} tssMatx/
  mv tssMatx/ ${prefix}/

  """ 
}


/***********
 * MultiQC *
 ***********/

process getSoftwareVersions{
  label 'python'
  label 'lowCpu'
  label 'lowMem'
  publishDir path: "${params.outDir}/software_versions", mode: "copy"

  when:
  !params.skipSoftVersions

  input:
  file ("v_bowtie2.txt") from chBowtie2Version.first().ifEmpty([])
  file("v_cutadapt.txt") from chCutadaptVersion.first().ifEmpty([])
  file("v_star.txt") from chStarVersion.first().ifEmpty([])
  file("v_samtools.txt") from chSamtoolsVersion.first().ifEmpty([])
  file("v_deeptools.txt") from chBamCoverageVersion.first().ifEmpty([])
  file("v_bedtools.txt") from chBedtoolsVersion.first().ifEmpty([])
  file("v_python.txt") from chPythonVersion.first().ifEmpty([])
  file("v_R.txt") from chRversion.first().ifEmpty([])

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
  label 'unix'
  label 'lowCpu'
  label 'lowMem'

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
  file ('software_versions/*') from softwareVersionsYaml.collect().ifEmpty([])
  file ('workflow_summary/*') from workflowSummaryYaml.collect()
  //Modules
  file ('star/*') from chAlignmentLogs.collect().ifEmpty([])
  //file ('trimming/*') from chTrimmedReadsLog.collect().ifEmpty([])
  file("bamToBigWig/*") from chBamToBigLog.collect().ifEmpty([])
  file ('index/*') from chIndexBowtie2Log.collect().ifEmpty([])
  //Logs
  file("bowtie2/*") from chBowtie2Log.collect().ifEmpty([])
  file("removeRtPcr/*") from chPcrRtCountsLog.collect().ifEmpty([])
  file("cellThresholds/*") from chRemoveDupBarcodeLog.collect().ifEmpty([])
  file("rmDup/*") from chRemoveDupLog.collect().ifEmpty([])
  // Weighted histogram
  file ('countUMI/*') from mqcDistribUMI.collect().ifEmpty([])


  output: 
  file splan
  file "*report.html" into multiqc_report
  file "*_data"

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_report" : "--filename report"
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  //isPE = params.singleEnd ? "" : "-p"
  //designOpts= params.design ? "-d ${params.design}" : ""
  modules_list = "-m custom_content -m bowtie2 -m star"
  """
  stat2mqc.sh ${splan} 
  mqc_header.py --splan ${splan} --name "scChIP-seq" --version ${workflow.manifest.version} ${metadataOpts} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c multiqc-config-header.yaml -c $multiqcConfig $modules_list
  """
}

/****************
 * Sub-routines *
 ****************/
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
