  removeRTduplicates.sh ${prefix} ${task.cpus} ${flaggedBam} ######### les mettre en arguments dans le script 


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