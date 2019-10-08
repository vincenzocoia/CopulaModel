# KL divergence dcop1 vs dcop2 
library(CopulaModel)
library(cubature)

tau=0.5
th.bvn=bvn.b2cpar(tau)
th.pla=depmeas2cpar(tau,"tau","plackett")  # 11.40486
th.frk=depmeas2cpar(tau,"tau","frank")  # 5.736287

cat("\nBVN vs Plackett, tau=0.5\n")
KLcopvsbvn(rh=th.bvn,dcop2=dpla,param2=th.pla,copname2="Plackett",UB=7,iprint=T)
cat("\nBVN vs Frank, tau=0.5\n")
KLcopvsbvn(rh=th.bvn,dcop2=dfrk,param2=th.frk,copname2="Frank",UB=7,iprint=T)
cat("\nPlackett vs Frank, tau=0.5\n")
KLcopvscop("Plackett",th.pla,dcop1=dpla,"Frank",th.frk,dcop2=dfrk,UB=7,iprint=T)


