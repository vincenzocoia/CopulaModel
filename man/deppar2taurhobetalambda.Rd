\name{deppar2taurhobetalambda}
\Rdversion{1.1}
\alias{deppar2taurhobetalambda}
\alias{bvn.b2cpar}
\alias{bvn.cpar2b}
\alias{bvn.cpar2rhoS}
\alias{bvn.cpar2tau}
\alias{bvn.rhoS2cpar}
\alias{bvn.tau2cpar}
\alias{bvt.b2cpar}
\alias{bvt.cpar2b}
\alias{bvt.cpar2lm}
\alias{bvt.cpar2tau}
\alias{bvt.lm2cpar}
\alias{bvt.tau2cpar}
\alias{pla.b2cpar}
\alias{pla.cpar2rhoS}
\alias{pla.rhoS2cpar}
\alias{frk.b2cpar}
\alias{frk.cpar2rhoS}
\alias{frk.cpar2tau}
\alias{mtcj.b2cpar}
\alias{mtcj.cpar2lm}
\alias{mtcj.cpar2tau}
\alias{mtcj.lm2cpar}
\alias{mtcj.tau2cpar}
\alias{joe.b2cpar}
\alias{joe.cpar2lm}
\alias{joe.cpar2tau}
\alias{joe.lm2cpar}
\alias{gum.b2cpar}
\alias{gum.cpar2lm}
\alias{gum.cpar2rhoS}
\alias{gum.cpar2tau}
\alias{gum.lm2cpar}
\alias{gum.tau2cpar}
\alias{gal.b2cpar}
\alias{gal.cpar2lm}
\alias{gal.cpar2rhoS}
\alias{gal.cpar2tau}
\alias{gal.lm2cpar}
\alias{hr.b2cpar}
\alias{hr.cpar2lm}
\alias{hr.cpar2rhoS}
\alias{hr.cpar2tau}
\alias{hr.lm2cpar}
\alias{bb1.b2cpar}
\alias{bb1.cpar2lm}
\alias{bb1.cpar2tau}
\alias{bb1.lm2cpar}
\alias{bb2.cpar2tau}
\alias{bb3.cpar2tau}
\alias{bb4.cpar2lm}
\alias{bb4.lm2cpar}
\alias{bb5.cpar2rhoS}
\alias{bb5.cpar2tau}
\alias{bb6.cpar2tau}
\alias{bb7.cpar2lm}
\alias{bb7.cpar2tau}
\alias{bb7.lm2cpar}
\alias{bb8.cpar2tau}
\alias{bb9.cpar2tau}
\alias{bb10.cpar2tau}
\alias{tev.cpar2b}
\alias{tev.cpar2lm}
\alias{tev.cpar2rhoS}
\alias{tev.cpar2tau}
\title{
Bivariate copulas: mappings of copula parameter and dependence measures
}
\description{
Bivariate copulas: mappings of copula parameter (cpar) and
Kendall's tau, Spearman's rho, Blomqvist's beta, tail dependence lambda
where these are simple.
}
\usage{
bvn.b2cpar(b)  # b=beta
bvn.cpar2b(rho) 
bvn.cpar2rhoS(rho) 
bvn.cpar2tau(rho) 
bvt.cpar2b(cpar) 
bvt.cpar2lm(cpar)
bvt.lm2cpar(lm,nu)  # lm=lambda
pla.b2cpar(b) 
pla.cpar2rhoS(cpar)
pla.rhoS2cpar(rhoS, cpar0=1.5,mxiter=25,eps=1.e-6,iprint=F,mxstep=10)
frk.b2cpar(b, cpar0=0,mxiter=20,eps=1.e-8,iprint=F)
frk.cpar2rhoS(cpar)
frk.cpar2tau(cpar)
mtcj.b2cpar(b, cpar0=0,mxiter=20,eps=1.e-8,iprint=F)
mtcj.cpar2lm(cpar) 
mtcj.cpar2tau(cpar) 
mtcj.lm2cpar(lm) 
mtcj.tau2cpar(tau) 
joe.b2cpar(b, cpar0=0,mxiter=20,eps=1.e-8,iprint=F)
joe.cpar2lm(cpar) 
joe.cpar2tau(cpar)
joe.lm2cpar(lm) 
gum.b2cpar(b)
gum.cpar2lm(cpar) 
gum.cpar2rhoS(cpar) 
gum.cpar2tau(cpar) 
gum.lm2cpar(lm) 
gum.tau2cpar(tau) 
gal.b2cpar(b)
gal.cpar2lm(cpar) 
gal.cpar2rhoS(cpar) 
gal.cpar2tau(cpar) 
gal.lm2cpar(lm) 
hr.b2cpar(b)
hr.cpar2lm(cpar) 
hr.cpar2rhoS(cpar) 
hr.cpar2tau(cpar) 
hr.lm2cpar(lm) 
bb1.b2cpar(b,de,thstart=1,mxiter=30,eps=1.e-6)
bb1.cpar2lm(cpar)
bb1.cpar2tau(cpar) 
bb1.lm2cpar(lmpar)  
#and other similar functions (see Details below)
}
\arguments{
\item{cpar}{copula parameter (scalar or vector), 
rho and df for bvt.cpar2lm}
\item{rho}{parameter for bivariate normal or t }
\item{tau}{Kendall's tau for the copula }
\item{rhoS}{Spearman's rho for the copula }
\item{b}{Blomqvist's beta for the copula }
\item{de}{second parameter for BB1 copula }
\item{lm}{tail dependence for the copula if just upper (lower) tail dependence}
\item{nu}{degree of freedom parameter for bivariate t}
\item{lmpar}{vector with (lml,lmu) for lower and upper tail dependence 
parameters }
\item{cpar0}{starting point for Newton-Raphson iterations;
there is a default starting point in some cases if cpar=0 is specified}
\item{mxiter}{maximum number of iterations}
\item{eps}{tolerance for convergence}
\item{iprint}{print flag for iterations}
\item{mxstep}{bound on step size for Newton-Raphson iterations}
\item{thstart}{starting point for Newton-Raphson iterations for bb1.b2cpar}
}
\details{
For abbreviations of names of copula families (before the dot in function
names), see pcop help page.

Function names after the dot are abbreviations:

b2cpar for Blomqvist's beta to copula parameter
(also for bvt)

lm2cpar for tail dependence lambda to copula parameter
(also for  bb4, bb7) 

rhoS2cpar for Spearman's rho to copula parameter

tau2cpar for Kendall's tau to copula parameter
(also for ipsA) 

cpar2b for copula parameter to Blomqvist's beta
(also for tev)

cpar2lm for copula parameter to tail dependence lambda
(also for tev, bvt, bb4, bb7)

cpar2rhoS for copula parameter to Spearman's rho
(also for tev, bb5)

cpar2tau for copula parameter to Kendall's tau
(also for tev, bvt, bb2, bb3, bb5, bb6, bb7, bb8, bb9, bb10, imitlefA, ipsA)
}
\value{
copula parameter or one of tau, rhoS, beta, lambda
}
\seealso{
\code{\link{pcop}}
}
\examples{
pla.rhoS2cpar(.5,cpar0=1.84)
be=seq(.1,.9,.1)
cpar=pla.b2cpar(be)
beta=4*ppla(.5,.5,cpar)-1
print(rbind(cpar,beta))
bb1.lm2cpar(c(.4,.6))
frk.cpar2tau(2)
}
\keyword{copula}
\keyword{dependence}
