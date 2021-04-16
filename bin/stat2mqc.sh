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
echo "Sample_id,Sample_name,#Cells,#Cell>1000reads,Median (>1000reads), %Aligned,%Aligned_Barcoded,%Unique_Reads" > scChIPseq_table.csv

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
    reads_after_pcr_rt_rm=$( echo "scale=2; ($reads_after_pcr_rt_rm - $R1_mapped_R2_unmapped)" | bc -l ) 
    index_1_2_not_3=$( echo "scale=2; ($match_index_1_2 - $match_barcode)" | bc -l ) 
    index_1_not_2_not_3=$( echo "scale=2; ( $match_index_1 - $index_1_2_not_3 - $match_barcode)" | bc -l ) 
    index_2_not_1_3=$( echo "scale=2; ( $match_index_2 - $match_index_1_2)" | bc -l ) 
    index_3_not_1_2=$( echo "scale=2; ( $match_index_3 - $match_barcode)" | bc -l ) 
    no_index_found=$( echo "scale=2; ( $total_reads - $match_barcode - $index_1_2_not_3 - $index_1_not_2_not_3 - $index_2_not_1_3 - $index_3_not_1_2)" | bc -l ) 
    uniquely_mapped_and_barcoded_percent=$( echo "scale=2; ( 100*$uniquely_mapped_and_barcoded / $total_reads)" | bc -l ) 
    unique_reads_percent=$( echo "scale=2; ( 100*$unique_reads/ $total_reads)" | bc -l ) 
    # Table
    echo "${sample},$sname,$match_barcode,$index_1_2_not_3,$index_1_not_2_not_3,$index_2_not_1_3,$index_3_not_1_2,$no_index_found" >> scChIPseq_barcode.csv

    ## Data for mapping - STAR
    total_mapped=$( echo "scale=2; ($uniquely_mapped + $multimapped + $multimapped_toomany)" | bc -l ) 
    unmapped_count=$( echo "scale=2; ($total_reads - $total_mapped)" | bc -l ) 
    total_unmapped_percent=$( echo "scale=2; ($unmapped_mismatches_percent + $unmapped_tooshort_percent + $unmapped_other_percent)" | bc -l ) 
    uniquely_mapped_unbarcoded=$( echo "scale=2; ( $uniquely_mapped - $uniquely_mapped_and_barcoded)" | bc -l ) 
    multimapped=$( echo "scale=2; ( $multimapped + $multimapped_toomany)" | bc -l ) 
    unmapped=$unmapped_count
    # Table
    echo "${sample},$sname,$unique_reads,$R2_unmapped_duplicates,$rt_duplicates,$pcr_duplicates,$uniquely_mapped_unbarcoded,$multimapped,$unmapped" >> scChIPseq_alignments.csv

    ## Data for cell thresholds
    nbCell=$(wc -l < cellThresholds/${sample}_rmDup.count)
    n1000=$( sed 's/^\s*//g' cellThresholds/${sample}_rmDup.count | awk -v limit=1000 '$1>=limit && NR>1{c++} END{print c+0}')

    if (( $n1000>1 ))
    then 
        awk -v limit=1000 '$1>=limit && NR>1 {print $1}' cellThresholds/${sample}_rmDup.count | sort -n | uniq > list
        nb_lines=$(wc -l < list)
        mod=$(($nb_lines%2))
        if (( $mod == 0 ))
        then
            line_first=$(( $nb_lines/2 ))
            line_sec=$(( $line_first+1 ))
            first_num=$(sed "${line_first}q;d" list)
            sec_num=$(sed "${line_sec}q;d" list)
            median=$( echo "scale=0; (($first_num+$sec_num)/2)" | bc -l )
        else
            median=$( echo "scale=0; (($nbcell+1)/2)" | bc -l )
        fi
    else
        median=$n1000
    fi
    
    ## Summary table
    echo "${sample},$sname,$nbCell,$n1000,$median,$uniquely_mapped_percent,$uniquely_mapped_and_barcoded_percent,$unique_reads_percent" >> scChIPseq_table.csv
done

