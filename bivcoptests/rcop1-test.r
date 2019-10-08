# random generation for 1-parameter copula families

library(CopulaModel)

sp=.5

n=100
#n=100000

set.seed(123)
pla.par=pla.rhoS2cpar(sp)
out=rpla(n,pla.par,icheck=T)

set.seed(123)
frk.par=depmeas2cpar(sp,"rhoS","frank")
out=rfrk(n,frk.par,icheck=T)

set.seed(123)
mtcj.par=depmeas2cpar(sp,"rhoS","mtcj")
out=rmtcj(n,mtcj.par,icheck=T)

set.seed(123)
joe.par=depmeas2cpar(sp,"rhoS","joe")
out=rjoe(n,joe.par,icheck=T)

set.seed(123)
gum.par=depmeas2cpar(sp,"rhoS","gumbel")
out=rgum(n,gum.par,icheck=T)

set.seed(123)
gal.par=depmeas2cpar(sp,"rhoS","galambos")
out=rgal(n,gal.par,icheck=T)

set.seed(123)
hr.par=depmeas2cpar(sp,"rhoS","huesler-reiss") 
hr.par=depmeas2cpar(sp,"rhoS","hr") 
out=rhr(n,hr.par,icheck=T)

