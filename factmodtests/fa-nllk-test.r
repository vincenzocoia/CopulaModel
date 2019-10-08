# run factanal and then get nllk, 
# compare with mvtpfact with 1 and 2 factors
library(CopulaModel) 

data(euro07gf);
zdat=euro07gf$zscore
rr=cor(zdat)
d=ncol(zdat)
n=nrow(zdat)

start=rep(0.4,d)
st2=rep(0.4,2*d)

out1=mvtpfact(zdat,start,df=350,pfact=1,prlevel=0)
print(out1$estimate)
# 0.7508538 0.9357568 0.9457988 0.9703055 0.8975374 0.9501660 0.8146135
# convert to nllk on N(0,1) and not U(0,1) scale
adj=.5*sum(zdat^2)+.5*n*d*log(2*pi)
out1$znllk=out1$min+adj

#out2=ml2mvtfact(st2,df=350,ifixed=rep(F,2*d),zdat,prlevel=0)
out2=mvtpfact(zdat,st2,df=350,pfact=2,prlevel=0)
print(out2$estimate)
# [1] 0.6631025 0.7521069 0.7577258 0.7598880 0.7369290 0.7841513 0.9990000
# [8] 0.4804267 0.8424189 0.8650239 0.9392054 0.7577941 0.8621264 0.4188591
out2$znllk=out2$min+adj
print(c(out1$znllk,out2$znllk,out1$znllk-out2$znllk))
# 1278.992168 1271.333944    7.658225

# these are a bit different from result from factanal because 
# diag of covariance matrix is not 1
# but out1-out2 should compare with nllk difference of factanal

out=factanal2nllk(rr,3,n,iprint=T)
# max abs diff  0.05753447  for  1 factors
# max abs diff  0.02526525  for  2 factors
# max abs diff  0.008528348 for  3 factors

# max is 3 factors for dimension d=7
print(out)
# 1278.992 1271.391 1266.753
print(-diff(out)) # matches results of mvtpfact 
# [1] 7.600734 4.638273

