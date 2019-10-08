# best gaussian truncated vines with 1,2,3,...d-2 trees

library(CopulaModel)

# polychoric correlation matrix for Science item response data
rmat=matrix(c(
1.0, 0.09947529, 0.201211361, 0.34626600, 0.08950463, 0.182484548, 0.40804986,
0.09947529, 1.0,-0.082776320,-0.02804555, 0.46381814, 0.410792782,-0.03652279,
0.20121136,-0.08277632, 1.0, 0.47861777,-0.10388833,-0.007641919, 0.20864149,
0.34626600,-0.02804555, 0.478617767, 1.0,-0.03596023, 0.102731469, 0.37689932,
0.08950463, 0.46381814,-0.103888330,-0.03596023, 1.0, 0.434794093,-0.01434073,
0.18248455, 0.41079278,-0.007641919, 0.10273147, 0.43479409, 1.0, 0.1179876,
0.40804986,-0.03652279, 0.208641491, 0.37689932,-0.01434073, 0.117987630, 1.0),
7,7)

out=gausstrvine(rmat,iprint=T)
d=nrow(rmat)

for(ell in 1:(d-2))
{ A=vnum2array(d,out$bnum[ell])
  cat("\ntruncation level ", ell,"\n")
  print(A)
  perm=out$permmat[,ell]
  Aperm=varrayperm(A,perm)
  print(Aperm)
  cat("check on log determinant\n")
  pcmat=out$pcarr[,,ell]
  logdet=sum(log(1-pcmat[1:ell,]^2))
  print(logdet)
}

