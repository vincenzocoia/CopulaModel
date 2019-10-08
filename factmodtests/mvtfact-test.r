# multivariate t with factor structure

library(CopulaModel)

data(euro07gf);
udat=euro07gf$uscore;
d=ncol(udat);
st1=rep(0.4,d);
st2=rep(0.4,2*d);

# Gaussian if df>300
for(df in c(3,5,10,15,20))
{ tdata=qt(udat,df);
  cat("\ndf=", df,"\n")
  cat("1-factor MVt\n")
  out1t=mvtpfact(tdata,st1,pfact=1,df=df,prlevel=1)
  cat("\n2-factor MVt\n")
  out2t=mvtpfact(tdata,st2,pfact=2,df=df,prlevel=1)
  cat("============================================================\n")
  st1=out1t$estimate
  st2=out2t$estimate
}

