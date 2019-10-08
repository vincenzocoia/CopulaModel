# DAX GARCH-filtered log return data; sequential MST for best truncated vines

library(CopulaModel)
library(igraph) 

load("dax1112gf.RData") # 
r1112=cor(dax1112gf$zscore)
nn=nrow(dax1112gf$zscore)
d=ncol(dax1112gf$zscore)

#label
# [1] "ALV.DE"  "CBK.DE"  "DB1.DE"  "DBK.DE"  "MUV2.DE" "BAS.DE"  "BAYN.DE"
# [8] "LIN.DE"  "SDF.DE"  "MRK.DE"  "SIE.DE"  "TKA.DE"  "HEI.DE"  "BMW.DE" 
#[15] "DAI.DE"  "VOW3.DE" "ADS.DE"  "BEI.DE"  "HEN3.DE" "MEO.DE"  "DPW.DE" 
#[22] "LHA.DE"  "EOAN.DE" "RWE.DE"  "FME.DE"  "FRE.DE"  "DTE.DE"  "IFX.DE" 
#[29] "SAP.DE"

bestmst=gausstrvine.mst(r1112,ntrunc=20,iprint=F)
save(file="dax1112-seqbestvines.RData",bestmst)
#names(bestmst)
#[1] "RVM"        "mst"        "treeweight" "trunclevel" "truncval" 
#names(bestmst$RVM)
# "VineA"  "pc"     "Matrix" "Cor"

absmax=function(x) { max(abs(x),na.rm=T) }

cat("\ndiagonal of vine array\n")
print(diag(bestmst$RVM$VineA))
cat("\nFirst 5 rows of vine array\n")
print(bestmst$RVM$VineA[1:5,])
cat("\nmax partial correlation by tree\n")
print(apply(bestmst$RVM$pc,1,absmax))


#diag(bestmst$RVM$VineA)
# [1] 23 27  1 24 20 22  6 12 15 11 13  9 28 29 14 16 17  3  4  2  5 21  7 10 8
#[26] 19 18 26 25

#apply(daxmst$RVM$pc,1,absmax)
# [1] 0.84867684 0.35548681 0.29550766 0.20796784 0.19805493 0.15500492
# [7] 0.14779944 0.13060986 0.12138667 0.17383070 0.11048147 0.14888916
#[13] 0.14961065 0.11334717 0.09976289 0.11136982 0.14776438 0.12733927
#[19] 0.16461645 0.07032140 0.13076549 0.09721258 0.10641595 0.14056496
#[25] 0.19917042 0.07916924 0.06275588 0.00000000       -Inf
