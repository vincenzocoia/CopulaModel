# check acf for AR(2) 

library(CopulaModel)

rh1=.6
rh2=.5
pc=(rh2-rh1^2)/(1-rh1^2)
ph1=rh1*(1-pc)
ph2=pc

rhv=rep(0,10)
rhv[1]=rh1; rhv[2]=rh2
for(i in 3:10)
{ rhv[i]=ph1*rhv[i-1]+ph2*rhv[i-2] }
print(rhv)

d=10
{ pcmat=matrix(0,d,d)
  cat("\nd=",d,"\n")
  for(j in 1:(d-1)) pcmat[j,j+1]=rh1
  for(j in 1:(d-2)) pcmat[j,j+2]=pc
  rmat=pcor2cor.dvine(pcmat)
  print(rmat)
}

# OK matches
