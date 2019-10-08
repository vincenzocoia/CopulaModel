# acf function for AR2
# rho1 = lag1 correlation
# rho2 = lag2 correlation
# d = dimension of desired Toeplitz matrix
# Output : acf function to lag d-1
ar2acf=function(rho1,rho2,d)
{ pc=(rho2-rho1^2)/(1-rho1^2)
  ph1=rho1*(1-pc)  # first regression coefficient 
  ph2=pc
  rhv=rep(0,d-1)
  rhv[1]=rho1; rhv[2]=rho2
  for(i in 3:(d-1))
  { rhv[i]=ph1*rhv[i-1]+ph2*rhv[i-2] }
  c(1,rhv)
}

