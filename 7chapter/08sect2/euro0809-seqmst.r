# sequential minimum spanning tree for euro0809 GARCH-filtered returns
#  transformed to normal scores.

# using gausstrvine.mst and igraph package

library(CopulaModel)
library(igraph) 

# normal scores correlation matrix for GARCH-filtered log return data
rmat=read.table("euro0809gfR.txt",skip=1,nrow=7,header=F)
rmat=as.matrix(rmat)
print(rmat)

#source("../../R/gausstrvineMST.R")
seqbest=gausstrvine.mst(rmat,ntrunc=6,iprint=T)

cat("\n============================================================\n")
cat("\nfinal result for sequential MST\n")
print(seqbest)
#save(seqbest,file="euro0809-seqbestvines.RData")
