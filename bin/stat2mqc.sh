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
# The column names have to be the same as the ID column in the multiqcConfig.yaml !!!!! 
echo -e "Sample_id,Sample_name,Cells>1000reads,Reads(median)/cell,Aligned,Aligned & Barcoded,Unique Reads" > scChIPseq_table.csv

for sample in $all_samples
do
    ## sample name
    sname=$(awk -F, -v sname=$sample '$1==sname{print $2}' $splan)

    # BOWTIE2
    match_index_1=$(grep -e "## Number of matched indexes 1:" bowtie2/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    match_index_2=$(grep -e "## Number of matched indexes 2:" bowtie2/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    match_index_1_2=$(grep -e "## Number of matched indexes 1 and 2:" bowtie2/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    match_index_3=$(grep -e "## Number of matched indexes 3:" bowtie2/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')
    match_barcode=$(grep -e "## Number of matched barcodes:" bowtie2/${sample}_bowtie2.log | sed 's/.*://g' | grep -o -e '[0-9]*\.*[0-9]*')

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
    multimapped=$(grep -e "Number of reads mapped to multiple loci " star/${sample}Log.final.out | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')
    multimapped_toomany=$(grep -e "Number of reads mapped to too many loci " star/${sample}Log.final.out | sed 's/.*|//g' | grep -o -e '[0-9]*\.*[0-9]*')    

    ## Data for the barcode matching graph
    #reads_after_pcr_rt_rm=$( echo "scale=2; ($reads_after_pcr_rt_rm - $R1_mapped_R2_unmapped)" | bc -l ) 
    reads_after_pcr_rt_rm=$(echo "$reads_after_pcr_rt_rm $R1_mapped_R2_unmapped" | awk ' { printf "%.2f", $1-$2 } ')
    #index_1_2_not_3=$( echo "scale=2; ($match_index_1_2 - $match_barcode)" | bc -l ) 
    index_1_2_not_3=$(echo "$match_index_1_2 $match_barcode" | awk ' { printf "%.2f", $1-$2 } ')
    #index_1_not_2_not_3=$( echo "scale=2; ( $match_index_1 - $index_1_2_not_3 - $match_barcode)" | bc -l ) 
    index_1_not_2_not_3=$(echo "$match_index_1 $index_1_2_not_3 $match_barcode" | awk ' { printf "%.2f", $1-$2-$3 } ')
    #index_2_not_1_3=$( echo "scale=2; ( $match_index_2 - $match_index_1_2)" | bc -l )
    index_2_not_1_3=$(echo "$match_index_2 $match_index_1_2" | awk ' { printf "%.2f", $1-$2 } ')
    #index_3_not_1_2=$( echo "scale=2; ( $match_index_3 - $match_barcode)" | bc -l ) 
    index_3_not_1_2=$(echo "$match_index_3 $match_barcode" | awk ' { printf "%.2f", $1-$2 } ')
    #no_index_found=$( echo "scale=2; ( $total_reads - $match_barcode - $index_1_2_not_3 - $index_1_not_2_not_3 - $index_2_not_1_3 - $index_3_not_1_2)" | bc -l ) 
    no_index_found=$(echo "$total_reads $match_barcode $index_1_2_not_3 $index_1_not_2_not_3 $index_2_not_1_3 $index_3_not_1_2" | awk ' { printf "%.2f", $1-$2-$3-$4-$5-$6 } ')
    #uniquely_mapped_and_barcoded_percent=$( echo "scale=2; ( 100*$uniquely_mapped_and_barcoded / $total_reads)" | bc -l )
    uniquely_mapped_and_barcoded_percent=$(echo "$uniquely_mapped_and_barcoded $total_reads" | awk ' { printf "%.2f", 100*$1/$2 } ')
    #unique_reads_percent=$( echo "scale=2; ( 100*$unique_reads/ $total_reads)" | bc -l )
    unique_reads_percent=$(echo "$unique_reads $total_reads" | awk ' { printf "%.2f", 100*$1/$2 } ')

    # Table
    echo "${sample},$sname,$match_barcode,$index_1_2_not_3,$index_1_not_2_not_3,$index_2_not_1_3,$index_3_not_1_2,$no_index_found" >> scChIPseq_barcode.csv

    ## Data for mapping - STAR
    #total_mapped=$( echo "scale=2; ($uniquely_mapped + $multimapped + $multimapped_toomany)" | bc -l )
    total_mapped=$(echo "$uniquely_mapped $multimapped $multimapped_toomany" | awk ' { printf "%.2f", $1+$2+$3 } ')

    #unmapped_count=$( echo "scale=2; ($total_reads - $total_mapped)" | bc -l )
    unmapped_count=$(echo "$total_reads $total_mapped" | awk ' { printf "%.2f", $1-$2 } ')

    #total_unmapped_percent=$( echo "scale=2; ($unmapped_mismatches_percent + $unmapped_tooshort_percent + $unmapped_other_percent)" | bc -l )
    total_unmapped_percent=$(echo "$unmapped_mismatches_percent $unmapped_tooshort_percent $unmapped_other_percent" | awk ' { printf "%.2f", $1+$2+$3 } ')

    #uniquely_mapped_unbarcoded=$( echo "scale=2; ( $uniquely_mapped - $uniquely_mapped_and_barcoded)" | bc -l ) 
    uniquely_mapped_unbarcoded=$(echo "$uniquely_mapped $uniquely_mapped_and_barcoded" | awk ' { printf "%.2f", $1-$2 } ')

    #multimapped=$( echo "scale=2; ( $multimapped + $multimapped_toomany)" | bc -l )
    multimapped=$(echo "$multimapped $multimapped_toomany" | awk ' { printf "%.2f", $1+$2 } ')

    unmapped=$unmapped_count
    # Table
    echo "${sample},$sname,$unique_reads,$R2_unmapped_duplicates,$rt_duplicates,$pcr_duplicates,$uniquely_mapped_unbarcoded,$multimapped,$unmapped" >> scChIPseq_alignments.csv

    ## Data for cell thresholds
    # total cells 
    nbCell=$(wc -l < cellThresholds/${sample}_rmDup.count)
    # nb cells with more than 1000 reads
    n1000=$( sed 's/^\s*//g' cellThresholds/${sample}_rmDup.count | awk -v limit=1000 '$1>=limit && NR>1{c++} END{print c+0}')

    # Median reads per cell with more than 1000 reads
    if (( $n1000>1 ))
    then 
        awk -v limit=1000 '$1>=limit && NR>1 {print $1}' cellThresholds/${sample}_rmDup.count | sort -n > list_nbReads_over1000
        mod=$(($n1000%2))
        if (( $mod == 0 ))
        then
            # get the first number that cut in two the number of read range
            line_first=$(( $n1000/2 ))
            first_num=$(sed "${line_first}q;d" list_nbReads_over1000)
            # get the second number that cut in two the number of read range
            line_sec=$(( $line_first+1 ))
            sec_num=$(sed "${line_sec}q;d" list_nbReads_over1000)
            
            #median=$( echo "scale=0; (($first_num+$sec_num)/2)" | bc -l )
            median=$(echo "$first_num $sec_num" | awk ' { printf "%.1f", ($1+$2)/2 } ')
        else
            #median=$( echo "scale=0; (($nbcell+1)/2)" | bc -l )
            line_median=$(( $n1000/2 ))
            median=$(sed "${line_median}q;d" list_nbReads_over1000)
        fi
    else
    # if there is one cell take the first line == number of reads within this cell
    # if there is no cell, it will write 0
        median=$(sed "1q;d" cellThresholds/${sample}_rmDup.count)
    fi
    
    ## Summary table
    echo -e "${sample},$sname,$n1000,$median,$uniquely_mapped_percent,$uniquely_mapped_and_barcoded_percent,$unique_reads_percent" >> scChIPseq_table.csv

done

