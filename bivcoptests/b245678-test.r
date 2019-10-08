# checks of 1-parameter copula families
# Plackett, MTCJ, Joe/B5, Gumbel, Galambos

library(CopulaModel)
#source("../R/copcond.R")

be=.5
#be=.8 # causes problem with HR
# test of beta to copula parameter for 1-parameter copula families
par.pla=pla.b2cpar(be)
par.mtcj=mtcj.b2cpar(be)
par.joe=joe.b2cpar(be)
par.gum=gum.b2cpar(be)
par.gal=gal.b2cpar(be)
par.hr=hr.b2cpar(be)
par.fgm=0.5
cat("beta=",be,"\n")

u=.3
v=seq(.4,.9,.1)
cat("\nChecking pcop, pcond, dcop\n")
chkcopderiv(u,v,par.pla,bcdf=ppla,pcond=pcondpla,bpdf=dpla,str="pla",eps=1.e-5)
chkcopderiv(u,v,par.mtcj,bcdf=pmtcj,pcond=pcondmtcj,bpdf=dmtcj,str="mtcj",eps=1.e-5)
chkcopderiv(u,v,par.joe,bcdf=pjoe,pcond=pcondjoe,bpdf=djoe,str="joe",eps=1.e-5)
chkcopderiv(u,v,par.gum,bcdf=pgum,pcond=pcondgum,bpdf=dgum,str="gum",eps=1.e-5)
chkcopderiv(u,v,par.gal,bcdf=pgal,pcond=pcondgal,bpdf=dgal,str="gal",eps=1.e-5)
chkcopderiv(u,v,par.hr,bcdf=phr,pcond=pcondhr,bpdf=dhr,str="hr",eps=1.e-5)
chkcopderiv(u,v,par.fgm,bcdf=pfgm,pcond=pcondfgm,bpdf=dfgm,str="fgm",eps=1.e-5)

set.seed(123)
u=runif(1)
v=runif(5)
chkcopderiv(u,v,par.pla,bcdf=ppla,pcond=pcondpla,bpdf=dpla,str="pla",eps=1.e-5)
chkcopderiv(u,v,par.mtcj,bcdf=pmtcj,pcond=pcondmtcj,bpdf=dmtcj,str="mtcj",eps=1.e-5)
chkcopderiv(u,v,par.joe,bcdf=pjoe,pcond=pcondjoe,bpdf=djoe,str="joe",eps=1.e-5)
chkcopderiv(u,v,par.gum,bcdf=pgum,pcond=pcondgum,bpdf=dgum,str="gum",eps=1.e-5)
chkcopderiv(u,v,par.gal,bcdf=pgal,pcond=pcondgal,bpdf=dgal,str="gal",eps=1.e-5)
chkcopderiv(u,v,par.hr,bcdf=phr,pcond=pcondhr,bpdf=dhr,str="hr",eps=1.e-5)
chkcopderiv(u,v,par.fgm,bcdf=pfgm,pcond=pcondfgm,bpdf=dfgm,str="fgm",eps=1.e-5)

cat("\nchecking pcond/qcond\n")
u=seq(.1,.6,.1)
v=seq(.4,.9,.1)
chkcopcond(u,v,par.pla,pcondpla,qcondpla,"pla")
chkcopcond(u,v,par.mtcj,pcondmtcj,qcondmtcj,"mtcj")
chkcopcond(u,v,par.joe,pcondjoe,qcondjoe,"joe")
chkcopcond(u,v,par.gum,pcondgum,qcondgum,"gum")
chkcopcond(u,v,par.gal,pcondgal,qcondgal,"gal")
# error here to check for large beta
chkcopcond(u,v,par.hr,pcondhr,qcondhr,"hr")
chkcopcond(u,v,par.fgm,pcondfgm,qcondfgm,"fgm")

set.seed(1234)
u=runif(3)
v=runif(3)
chkcopcond(u,v,par.pla,pcondpla,qcondpla,"pla")
chkcopcond(u,v,par.mtcj,pcondmtcj,qcondmtcj,"mtcj")
chkcopcond(u,v,par.joe,pcondjoe,qcondjoe,"joe")
chkcopcond(u,v,par.gum,pcondgum,qcondgum,"gum")
chkcopcond(u,v,par.gal,pcondgal,qcondgal,"gal")
chkcopcond(u,v,par.hr,pcondhr,qcondhr,"fgm")
chkcopcond(u,v,par.fgm,pcondfgm,qcondfgm,"fgm")

cat("\n============================================================\n")
