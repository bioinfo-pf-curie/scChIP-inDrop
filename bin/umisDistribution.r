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

# Aim weighted hist : better vizualisation
#  cellules avec peu de reads/umis >> cellules avec bcp reads/umis 
#  Pour voir les 2 distributions de manière égale on donne plus de poids aux cellules avec bcp reads/umis. 

# Arguments
countMatrix<- as.character(commandArgs(TRUE)[1])
prefix = as.character(commandArgs(TRUE)[2])

# Load data
matrix<-read.table(countMatrix, header=TRUE, sep = "\t")

# get a matrix with barcode names in the first column and the number of reads in the second
longMatrix<-data.frame(Barcodes=colnames(matrix[,-1]), nbReads=colSums(matrix[,-1]))

# Histogram and pdf export
pdf(paste0(as.character(prefix), 'distribution.pdf'))
# y = sum des reads par bin, x= #readsTot/cell (on ne sait cb de cellules ont entre x1 et x2 reads)
p<-weighted.hist(log10(longMatrix$nbReads), w=longMatrix$nbReads, ylab ="sum #reads/bin", xlab="log10(#reads/cell)")
# somme des reads par bin <=> somme des reads de toutes les cellules qui ont entre x1 et x2 
# Pour les pics les plus à gauche: pic haut = bcp de cellules qui ont peu de reads
# ==> Le weighted hist permet de montrer la proportion de cellules ayant un faible nombre de reads
dev.off()

# Export table to create the plot with multiqc

# 1st - calculate the center (== mean) of each bin because only a curve can be drawn with mqc
# To do so, the center of each bin is calculated
# Get first list of mean bins
breaks<-list(p[["breaks"]])
spitNum<-length(breaks[[1]])/2
splitedBreaks1 <- lapply(breaks, function(x) split(unlist(x), cut(seq_along(unlist(x)), spitNum, labels = F)))
splitedBreaks1<-unlist(splitedBreaks1, recursive = F)
x1<-lapply(splitedBreaks1, mean)
# Do the same but begining=2nd value to have means of missing bins
breaks<-list(p[["breaks"]][-c(1,length(breaks[[1]]))])
spitNum<-length(breaks[[1]])/2
splitedBreaks2 <- lapply(breaks, function(x) split(unlist(x), cut(seq_along(unlist(x)), spitNum, labels = F)))
splitedBreaks2<-unlist(splitedBreaks2, recursive = F)
x2<-lapply(splitedBreaks2, mean)
# Merge both x and sort result to obtain all the point of the x abscisse
xUnsorted<-unlist(append(x1,x2))
x<-sort(xUnsorted)

# Get y values == counts
y<-p["counts"]

weightedHist_DF<-data.frame(y=y,x=x)
write.table(weightedHist_DF, paste0(as.character(prefix), "_distDF.mqc"),
          sep=',', row.names=FALSE, col.names=FALSE)

