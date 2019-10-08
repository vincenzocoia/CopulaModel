# compare output and timing for Kendall tau computation: cor and taucor
library(CopulaModel)

# compare taucor (fast Kendall tau) with cor( method="kendall")
ktcor=function(x,y) { cor(x,y,method = "kendall") }
spcor=function(x,y) { cor(x,y,method = "spearman") }

# simulation to compare computation of Kendall's tau
# n = sample size
# nsim = number of simulation replications
# seed = seed for set.seed()
# fnname = function name, ktau or taucor
# rho = correlation for bivariate normal, used for simulated data
# Output: Kendall's tau is printed in the function
ktsimul=function(n,nsim,seed,fnname,rho=0.5)
{ set.seed(seed)
  r1=sqrt(1-rho^2)
  for(isim in 1:nsim)
  { x = rnorm(n,0,1)
    y = rho*x+ r1*rnorm(n,0,1)
    tau=fnname(x,y)
    cat(isim, tau,"\n")
  }
  tau
}

set.seed(123)
for(isim in 1:10)
{ n = sample(1:1000,1)
  print(n)
  x = rnorm(n,0,1)
  y = rnorm(n,0,1)
  tau1=ktcor(x,y); tau2=taucor(x,y)
  print(c(tau1,tau2))
}

n=1000; nsim=5;
n=10000; nsim=5;
rho=0.5
seed=123
for(rho in seq(-.9,.9,.2))
{ cat("\ntrue rho=", rho,"\n")
  cat("method=kendall O(n^2)\n")
  time0=proc.time()
  ktsimul(n,nsim,seed,ktcor,rho)
  time1 = proc.time() - time0
  time0=proc.time()
  cat("Knight's O(n*log(n)) implementation\n")
  ktsimul(n,nsim,seed,taucor,rho)
  time2 = proc.time() - time0
  # user  system elapsed
  cat(time1,time2,"\n")
}
