# check of function bb1.tau2eqlm
# this is useful for MLE with BB1 copulas in vines
library(CopulaModel)

for(tau in seq(.1,.9,.1))
{ cat("\ntau=", tau,"\n")
  out=bb1.tau2eqlm(tau,iprint=T)
  cat("(th,de)=", out,"\n")
  cat("backcheck\n")
  cat(bb1.cpar2lm(out[1:2]),"\n")
  cat(bb1.cpar2tau(out[1:2]),"\n")
}

