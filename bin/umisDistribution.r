#!/usr/bin/env Rscript
.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))

##############################################################################
#  This script read a full read/gene/barcode count matrix and plot its weighted distribution of reads per cell
#
#  Copyright (c) 2020 - Institut Curie
#
#  File author(s):
#      Louisa Hadj Abed <louisa.hadj-abed@curie.fr>
#
#  Distributed under the terms of the CeCILL-B license.
#  The full license is in the LICENSE file, distributed with this software.
#
##############################################################################
library(plotrix)
library(reshape2)
library(dplyr)
# Aim weighted hist : better vizualisation
#  cellules avec peu de reads/umis >> cellules avec bcp reads/umis 
#  Pour voir les 2 distributions de manière égale on donne plus de poids aux cellules avec bcp reads/umis. 

# Arguments
countMatrix<- as.character(commandArgs(TRUE)[1])
prefix = as.character(commandArgs(TRUE)[2])

# Load data
matrix<-read.table(countMatrix, header=FALSE, sep = "")

# get a matrix with barcode names in the first column and the number of reads in the second
# Only if more than one cell
if(nrow(matrix)>1){
    #longMatrix<-data.frame(Barcodes=colnames(matrix)[-1], nbReads=colSums(matrix[,-1]))
    colnames(matrix)<-c("nbReads", "Barcodes")
    
    # Histogram and pdf export
    pdf(paste0(as.character(prefix), 'distribution.pdf'))
    # y = sum des reads par bin, x= #readsTot/cell (on ne sait cb de cellules ont entre x1 et x2 reads)
    p<-weighted.hist(log10(matrix$nbReads), w=log10(matrix$nbReads), ylab ="Number of barcodes", xlab="log10(UMIs)")
    # somme des reads par bin <=> somme des reads de toutes les cellules qui ont entre x1 et x2 
    # Pour les pics les plus à gauche: pic haut = bcp de cellules qui ont peu de reads
    # ==> Le weighted hist permet de montrer la proportion de cellules ayant un faible nombre de reads
    dev.off()
    
    # Export table to create the plot with multiqc
    yy<-data.frame(p["counts"])
    # double all y values to create bin limits 
    y<-yy %>% slice(rep(1:n(), each = 2))
    
    xx<-data.frame(x=p[["breaks"]])
    # remove first and last value in x 
    x<-data.frame(x=p[["breaks"]][-c(1,length(p[["breaks"]]))])
    x<-rbind(xx,x)
    # increasing order
    x<-arrange(x, as.numeric(x))
    
    weightedHist_DF<-cbind(x, y)
    write.table(weightedHist_DF, paste0(as.character(prefix), "_distDF.mqc"),
                sep=',', row.names=FALSE, col.names=FALSE)
}else{ # If only one cell, create empty files
    pdf(paste0(as.character("prefix"), 'distribution.pdf'))
    dev.off()
    weightedHist_DF<-data.frame()
    write.table(weightedHist_DF, paste0(as.character(prefix), "_distDF.mqc"),
                sep=',', row.names=FALSE, col.names=FALSE)
}


