# compare two versions of Gauss-Legendre quadrature
library(CopulaModel)
library(statmod)

# nq = number of quadrature points for Gauss-Legendre
# Output: comparison of quadrature points and weights with function in statmod
chk=function(nq)
{ cat("\nChecking for nq=", nq,"\n")
  out=gauss.quad.prob(nq,dist="uniform")
  out2=gausslegendre(nq)
  print(out$nodes)
  print(out2$nodes)
  print(out$weights)
  print(out2$weights)
  print(max(abs(out$nodes-out2$nodes)))
  print(max(abs(out$weights-out2$weights)))
  print(sum(out2$weights))
  invisible(0)
}
  
# nq=7
#$nodes
# 0.02544604 0.12923441 0.29707742 0.50000000 0.70292258 0.87076559 0.97455396
#$weights
# 0.06474248 0.13985270 0.19091503 0.20897959 0.19091503 0.13985270 0.06474248

chk(7)
chk(8)
chk(15)
chk(21)
chk(25)
chk(31)
