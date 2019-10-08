# compare taucor (fast Kendall tau) with cor( method="kendall")
# in the case of ties

library(CopulaModel)
ktcor=function(x,y) { cor(x,y,method = "kendall") }
spcor=function(x,y) { cor(x,y,method = "spearman") }

# simulation to compare computation of Kendall's tau
# n = sample size
# nsim = number of simulation replications
# seed = seed for set.seed()
# fnname = function name, ktau or taucor
# rho = correlation for bivariate normal, used for simulated data
# iround = number of decimal places for rounding
# iprint = print flag for the generated data sets
# Output: Kendall's tau is printed in the function
ktsimul2=function(n,nsim,seed,fnname,rho=0.5,iround=1,iprint=F)
{ set.seed(seed)
  r1=sqrt(1-rho^2)
  for(isim in 1:nsim)
  { x = rnorm(n,0,1)
    y = rho*x+ r1*rnorm(n,0,1)
    x = round(x,iround)
    y = round(y,iround)
    tau=fnname(x,y)
    if(iprint) print(cbind(x,y))
    cat(isim, tau,"\n")
  }
  tau
}


set.seed(123)

n=1000; nsim=5;
n=40
rho=0.5
rho=0.7
rho=-0.4
rho=-0.8
seed=123
cat("\ntrue rho=", rho,"\n")
cat("method=kendall O(n^2)\n")
ktsimul2(n,nsim,seed,ktcor,rho,iprint=T)
cat("Knight's O(n*log(n)) implementation\n")
ktsimul2(n,nsim,seed,taucor,rho,iprint=T)

# matches with "method=kendall O(n^2)\n"
