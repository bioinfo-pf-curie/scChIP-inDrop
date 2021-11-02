#!/usr/bin/env Rscript
.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))

umiMatrix<- as.character(commandArgs(TRUE)[1])
dir_res= as.character(commandArgs(TRUE)[2])


info = file.info(umiMatrix)

if(info$size <=58){ # if empty file
    
    
}else{
    #####  Long to sparse matrix:
    ####-----------------------------------------
    library(Matrix)
    data<-read.csv(umiMatrix, sep = '\t', header = T)
    library(reshape2)
    rm(data)
    data1<-melt(data)
    data1$X=as.factor(data1$X)
    data1$variable=as.factor(data1$variable)
    data1$value=as.factor(data1$value)
    X <- with(data1, sparseMatrix(i=as.numeric(X),
                                  j=as.numeric(variable),
                                  x=as.numeric(value),
                                  dimnames=list(levels(X), levels(variable))
    )
    )
    
    #### Save results like 10X 
    ####-----------------------------------------
    library("DropletUtils")
    cell.ids <- levels(data1[,2])
    ngenes <- nrow(X)
    gene.ids <- paste0("ID", seq_len(ngenes))
    gene.symb <-levels(data1[,1])
    
    write10xCounts(dir_res, X, gene.id=gene.ids, 
                   gene.symbol=gene.symb, barcodes=cell.ids)
}







