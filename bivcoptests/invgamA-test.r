# Archimedean copula with inverse gamma LT:
# LT involves Bessel function besselK
library(CopulaModel)

cpar=0.5
cpar=1.5
cat("cpar=", cpar, "tau =", invgamA.cpar2tau(cpar),"\n")
u=seq(.1,.6,.1)
v=seq(.4,.9,.1)
print(pinvgamA(u,v,cpar))
print(pcondinvgamA(v,u,cpar))
print(dinvgamA(u,v,cpar))

cat("\ncheck pcop, pcond, dcop\n")
u=.3
#u=.8
v=seq(.4,.9,.1)
chkcopderiv(u,v,cpar,bcdf=pinvgamA,pcond=pcondinvgamA,bpdf=dinvgamA,str="invgamA",eps=1.e-5)

