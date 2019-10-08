# compare multivariate cdf for Gumbel and Galambos copulas

library(CopulaModel)

# wrapper function for comparisons
# tau = Kendall tau value in (0,1)
# d = dimension
# Output within function: output of multivariate Gumbel and Galambos
#    cdfs at some random d-tuples of uniform scores
wrapcmp=function(tau,d)
{ cat("tau=", tau, " d=", d,"\n")
  cpargum=gum.tau2cpar(tau)
  cpargal=depmeas2cpar(tau,"tau","galambos")
  cat(cpargum,cpargal,"\n")
  set.seed(123)
  cat("uvec pmgum pmgal\n")
  for(i in 1:5)
  { uu=runif(d,0,1)
    tem1=pmgum(uu,cpargum)
    tem2=pmgal(uu,cpargal)
    cat(uu," : ", tem1,tem2,"\n")
  }
  cat("\n============================================================\n\n")
  invisible(0)
}

wrapcmp(.2,2)
wrapcmp(.2,3)
wrapcmp(.2,4)
wrapcmp(.2,5)
wrapcmp(.5,2)
wrapcmp(.5,3)
wrapcmp(.5,4)
wrapcmp(.5,5)
wrapcmp(.7,2)
wrapcmp(.7,3)
wrapcmp(.7,4)
wrapcmp(.7,5)

# Galambos cdf a little larger when tau is fixed 
