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
    data<-read.csv(umiMatrix, sep = '\t', header = F)
    data$V1=as.factor(data$V1)
    data$V2=as.factor(data$V2)
    data$V3=as.factor(data$V3)
    X <- with(data, sparseMatrix(i=as.numeric(V1),
                                 j=as.numeric(V2),
                                 x=as.numeric(V3),
                                 dimnames=list(levels(V1), levels(V2))
    )
    )
    
    #### Save results like 10X 
    ####-----------------------------------------
    library("DropletUtils")
    cell.ids <- levels(data[,2])
    ngenes <- nrow(X)
    gene.ids <- paste0("ID", seq_len(ngenes))
    gene.symb <-levels(data[,1])
    
    write10xCounts(dir_res, X, gene.id=gene.ids, 
                   gene.symbol=gene.symb, barcodes=cell.ids)
}







