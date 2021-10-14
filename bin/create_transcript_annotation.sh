#!/bin/bash

# Create an annotation based on gencode gtf 
# Take all transcripts for protein coding & long non coding RNA

annotation_gtf=$1

name=$(basename $annotation_gtf .gtf)

#Keep only transcripts
awk -v FS="\t" -v OFS="\t" '$3=="transcript"{print $0}' $annotation_gtf > tmp1

#Format to BED
awk -v FS="\t" -v OFS="\t" '{
  split($9,array,";")
  gsub(/ gene_type "/,"",array[3])
  gsub(/"/,"",array[3])
  gsub(/ gene_name "/,"",array[4])
  gsub(/"/,"",array[4]) 
  gsub(/ transcript_name "/,"",array[6])
  gsub(/"/,"",array[6])
  print $1,$4,$5,array[3],array[6],array[4],$7}'  tmp1 > tmp2
 

#Keep only protein coding, lncRNA  & remove column
awk -v FS="\t" -v OFS="\t" '$4=="lncRNA"||$4=="protein_coding"{print $1,$2,$3,$5,$6,$7}' tmp2 > gencode.v34.annotation.bed

# add +- 5000 around geneTSS
awk -v FS="\t" -v OFS="\t" '{
  if($6=="+"){
    print $1,$2-5000,$2+5000,$4,$5,$6
    
  } else{
    print $1,$3-5000,$3+5000,$4,$5,$6
  }
}' gencode.v34.annotation.bed > tmp3.bed

awk -v FS="\t" -v OFS="\t" '
BEGIN{
  transcript=""
  first_start=$2
  last_end=""
  last_gene_name=""
  last_chr=""
  last_strand=""
}
{
  if($2<last_end && $5==last_gene_name){
    if(transcript!="") {transcript=transcript","$4}
    else {
    transcript=$4
    if($3>$last_end){last_end=$3}
    }
  } else {
    if(transcript==""){
      if(NR>1) {print $0}
    }
    else {
      print last_chr,first_start,last_end,transcript,last_gene_name,last_strand
    }
    transcript=$4
    first_start=$2
    last_end=$3
  }
  last_gene_name=$5
  last_chr=$1
  last_strand=$6
}' tmp3.bed > tmp4

awk -v OFS="\t" '{if($2<0){$2=0}; print $0}' tmp4 $name.bed
rm tmp*
