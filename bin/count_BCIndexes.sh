#!/bin/bash

for file in *NOTjoinedBarcodes.txt
do 
    sample=`basename $file _NOTjoinedBarcodes.txt`

    # Two indexes only
    awk -v count=0 '$2!="*" && $3!="*" && $4=="*" {letterFirst=substr($2,1,1); letterSec=substr($3,1,1); count++}; END{print letterFirst letterSec "," count }' $file > $sample"double_counts.txt"
    awk -v count=0 '$2=="*" && $3!="*" && $4!="*" {letterFirst=substr($3,1,1); letterSec=substr($4,1,1); count++}; END{print letterFirst letterSec "," count }' $file >> $sample"double_counts.txt"   
    awk -v count=0 '$2!="*" && $3=="*" && $4!="*" {letterFirst=substr($2,1,1); letterSec=substr($4,1,1); count++}; END{print letterFirst letterSec "," count }' $file >> $sample"double_counts.txt"    
    # Reorganize
    awk '{gsub("CB","BC",$1); gsub("DC","CD",$1); gsub("BD","DB",$1) } 1' $sample"double_counts.txt" | sort > $sample"_indexes_mqc.log"

    # One index only
    awk -v count=0 '$2=="*" && $3=="*" && $4 != "*" {letter=substr($4,1,1); count++}; END{print letter "," count }' $file > one_counts.txt  
    awk -v count=0 '$2=="*" && $3!="*" && $4 == "*" {letter=substr($3,1,1);count++}; END{print letter "," count }' $file >> one_counts.txt   
    awk -v count=0 '$2!="*" && $3=="*" && $4 == "*" {letter=substr($2,1,1);count++}; END{print letter "," count }' $file >> one_counts.txt  

    sort one_counts.txt >> $sample"_indexes_mqc.log"
    
    # None
    awk -v count=0 '$2=="*" && $3=="*" && $4=="*" {count++}; END{print "None," count }' $file >> $sample"_indexes_mqc.log"

done



