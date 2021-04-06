#!/bin/bash

splan=$1

## Catch sample names
all_samples=$(awk -F, '{print $1}' $splan)

## Table headers
## Barcodes
echo "Sample_id,Sample_name,Barcoded,Index 1 and 2 found not 3,Index 1 found not 2 and 3,Index 2 found not 1 and 3,Index 3 found not 1 and 2,No Index Found ~ genomic DNA" > scChIPseq_barcode.csv
## Mapping
echo "Sample_id,Sample_name,Deduplicated reads, Window duplicates,RT duplicates,PCR duplicates,Uniquely mapped not barcoded,Mapped to multiple loci,Unmapped" > scChIPseq_alignments.csv
## Summary table
echo "Sample_id,Sample_name,%Aligned,%Aligned_Barcoded,%Unique_Reads" > scChIPseq_table.csv

for sample in $all_samples
do
    ## sample name
    sname=$(awk -F, -v sname=$sample '$1==sname{print $2}' $splan)

    # BOWTIE2
    match_index_1=$(grep -e "## Number of matched indexes 1:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    match_index_2=$(grep -e "## Number of matched indexes 2:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    match_index_1_2=$(grep -e "## Number of matched indexes 1 and 2:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    match_index_3=$(grep -e "## Number of matched indexes 3:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    match_barcode=$(grep -e "## Number of matched barcodes:" $scChIPseq_logfile | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')

    # Remove RT & PCR duplicats
    uniquely_mapped_and_barcoded=$(grep -e "## Number of reads mapped and barcoded:" removeRtPcr/${sample}_removePcrRtDup.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    pcr_duplicates=$(grep -e "## Number of pcr duplicates:" removeRtPcr/${sample}_removePcrRtDup.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    rt_duplicates=$(grep -e "## Number of rt duplicates:" removeRtPcr/${sample}_removePcrRtDup.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    R1_mapped_R2_unmapped=$(grep -e "## Number of R1 mapped but R2 unmapped:" removeRtPcr/${sample}_removePcrRtDup.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    reads_after_pcr_rt_rm=$(grep -e "## Number of reads after PCR and RT removal (not R1 unmapped R2):" removeRtPcr/${sample}_removePcrRtDup.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    
    # rmDup
    R2_unmapped_duplicates=$(grep -e "## Number of duplicates:" rmDup/${sample}_rmDup.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    unique_reads=$(grep -e "## Number of reads after duplicates removal:" rmDup/${sample}_rmDup.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')

    # STAR 
    total_reads=`grep -e "Number of input reads " star/${sample}Log.final.out | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*' `
    uniquely_mapped=`grep "Uniquely mapped reads number" star/${sample}Log.final.out | awk '{print $NF}'`
    uniquely_mapped_percent=`grep "Uniquely mapped reads %" star/${sample}Log.final.out | awk '{print $NF}' | sed -e 's/%//'`


    ## Data for the barcode matching graph
    reads_after_pcr_rt_rm=$( calc $reads_after_pcr_rt_rm - $R1_mapped_R2_unmapped)
    index_1_2_not_3=$( calc $match_index_1_2 - $match_barcode)
    index_1_not_2_not_3=$( calc $match_index_1 - $index_1_2_not_3 - $match_barcode)
    index_2_not_1_3=$( calc $match_index_2 - $match_index_1_2)
    index_3_not_1_2=$( calc $match_index_3 - $match_barcode)
    no_index_found=$( calc $total_reads - $match_barcode - $index_1_2_not_3 - $index_1_not_2_not_3 - $index_2_not_1_3 - $index_3_not_1_2)
    uniquely_mapped_and_barcoded_percent=$( calc 100*$uniquely_mapped_and_barcoded / $total_reads)
    unique_reads_percent=$( calc 100*$unique_reads/ $total_reads)
    # Table
    echo "${sample},$sname,$match_barcode,$index_1_2_not_3,$index_1_not_2_not_3,$index_2_not_1_3,$index_3_not_1_2,$no_index_found" >> scChIPseq_barcode.csv

    ## Data for mapping 
    uniquely_mapped_unbarcoded=$( calc $uniquely_mapped - $uniquely_mapped_and_barcoded)
    multimapped=$( calc $multimapped + $multimapped_toomany)
    unmapped=$unmapped_count
    # Table
    echo "${sample},$sname,$unique_reads,$R2_unmapped_duplicates,$rt_duplicates,$pcr_duplicates,$uniquely_mapped_unbarcoded,$multimapped,$unmapped" >> scChIPseq_alignments.csv

    ## Data for cell thresholds
    n100=$( sed 's/^\s*//g' cellThresholds/${sample}_rmDup.count | awk -v limit=100 '$1>=limit && NR>1{c++} END{print c+0}')
    n500=$( sed 's/^\s*//g' cellThresholds/${sample}_rmDup.count | awk -v limit=500 '$1>=limit && NR>1{c++} END{print c+0}' )
    n1000=$( sed 's/^\s*//g' cellThresholds/${sample}_rmDup.count | awk -v limit=1000 '$1>=limit && NR>1{c++} END{print c+0}')
    n1500=$( sed 's/^\s*//g' cellThresholds/${sample}_rmDup.count | awk -v limit=1500 '$1>=limit && NR>1{c++} END{print c+0}')
    desc="Number of barcodes with more than 100 unique reads = $n100 <br>Number of barcodes with more than 500 unique reads = $n500 <br>Number of barcodes with more than 1000 unique reads = $n1000 <br>Number of barcodes with more than 1500 unique reads = $n1500" 
    sed -i "s|{desc}|$desc|g" ../assets/multiqcConfig.yaml

    ## Summary table
    echo "${sample},$sname,$uniquely_mapped_percent,$uniquely_mapped_and_barcoded_percent,$unique_reads_percent" >> scChIPseq_table.csv

done

