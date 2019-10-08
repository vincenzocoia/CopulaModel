\name{copMarkovts}
\Rdversion{1.1}
\alias{mltscop1}
\alias{mltscop2}
\alias{rmspe.mccop1}
\alias{rmspe.mccop2}
\alias{rmspe.tvn}
\alias{nb1mcnllk}
\alias{nb2mcnllk}
\alias{gp1mcnllk}
\alias{gp2mcnllk}
\alias{pomcnllk}
\alias{nb1mc2nllk}
\alias{nb2mc2nllk}
\alias{gp1mc2nllk}
\alias{gp2mc2nllk}
\alias{pomc2nllk}
\alias{nb1ar2nllk}
\alias{nb2ar2nllk}
\alias{gp1ar2nllk}
\alias{gp2ar2nllk}
\alias{poar2nllk}
\title{
maximum likelihood estimation and
negative log-likelihood for copula-based Markov count time series models
}
\description{
maximum likelihood estimation and
negative log-likelihood for copula-based Markov count time series models
}
\usage{
mltscop1(y,x,pcop=pbvncop,start,family="Po",iprint=F,cparlb=0,cparub=30,prlevel=1)
mltscop2(y,x,pcop3=pmxid3ls,pcop2=pmxid2ls,start,family="Po",iprint=F,cparlb=0,cparub=30,prlevel=1)
rmspe.mccop1(param,yy,xdat=0,pcop=pbvncop,family="Po",iprint=F)
rmspe.mccop2(param,yy,xdat=0,pcop3=pmxid3ls,pcop2=pmxid2ls,family="Po",iprint=F)
rmspe.tvn(param,yy,xdat=0,family="Po",iprint=F)
nb1mcnllk(param,yy,xdat=0,pcop=pbvncop,cparlb=0,cparub=30) # NB1, Markov order1
nb2mcnllk(param,yy,xdat=0,pcop=pbvncop,cparlb=0,cparub=30) # NB2, Markov order1
gp1mcnllk(param,yy,xdat=0,pcop=pbvncop,cparlb=0,cparub=30) # GP1, Markov order1
gp2mcnllk(param,yy,xdat=0,pcop=pbvncop,cparlb=0,cparub=30) # GP2, Markov order1
pomcnllk(param,yy,xdat=0,pcop=pbvncop,cparlb=0,cparub=30) # Poisson, Markov order1
# Markov order 2 count time series below
nb1mc2nllk(param,yy,xdat=0,pcop3=pmxid3ls,pcop2=pmxid2ls,cparlb=0,cparub=30)
nb2mc2nllk(param,yy,xdat=0,pcop3=pmxid3ls,pcop2=pmxid2ls,cparlb=0,cparub=30)
gp1mc2nllk(param,yy,xdat=0,pcop3=pmxid3ls,pcop2=pmxid2ls,cparlb=0,cparub=30)
gp2mc2nllk(param,yy,xdat=0,pcop3=pmxid3ls,pcop2=pmxid2ls,cparlb=0,cparub=30)
pomc2nllk(param,yy,xdat=0,pcop3=pmxid3ls,pcop2=pmxid2ls,cparlb=0,cparub=30)
# Markov order 2 count time series based on trivariate Gaussian below
nb1ar2nllk(param,yy,xdat)
nb2ar2nllk(param,yy,xdat)
gp1ar2nllk(param,yy,xdat)
gp2ar2nllk(param,yy,xdat)
poar2nllk(param,yy,xdat)
}
\arguments{
\item{y}{count times series vector (nx1)}
\item{x}{matrix of covariates (n x ncovar), where ncovar is number
of covariates/predictors}
\item{yy}{count times series vector}
\item{xdat}{matrix of covariates, 0 by default for no covariates}
\item{param}{parameter vector (includes univariate and dependence parameters)}
\item{start}{starting point of param for nlm optimization}
\item{pcop}{function for bivariate copula cdf for Markov order 1
models} 
\item{pcop3}{function for trivariate copula cdf for Markov order 2
models} 
\item{pcop2}{function for bivariate marginal copula cdf for Markov 
order 2 models} 
\item{family}{univariate count regression model: one of "Po", "NB1", "NB2", "GP1", "GP2"}
\item{cparlb}{lower bound for copula parameters}
\item{cparub}{upper bound for copula parameters}
\item{iprint}{print flag for of intermediate results}
\item{prlevel}{print.level for nlm()}
}
\details{
mltscop1 is a front end to rmspe.mccop1 and one of 
nb1mcnllk, nb1mcnllk, gp1mcnllk, gp2mcnllk, pomcnllk for maximum likelihood
of copula-based Markov order 1 count time series model.

mltscop2 is a front end to rmspe.mccop2 and one of 
nb1mc2nllk, nb1mc2nllk, gp1mc2nllk, gp2mc2nllk, pomcnllk for maximum likelihood
of copula-based Markov order 2 count time series model.

}
\value{negative log-likelihood for any of the nllk functions,

$rmse and $pred (root mean square prediction error and 1-step predictions
for any of the rmspe functions,

$nllk, $mle (maximum likelihood estimate), $acov
(asymptotic covariance matrix, $rmpse for the mltscop1, mltscop2 functions
}
\examples{
data(wcb8594)
y=wcb8594$nclaims
cat("Frank NB2, 2 covariates\n")
ucpar2=c(1.7,-.2,-.2,.4,1.5)
frknb2=mltscop1(y,wcb8594[,2:3],pfrk,start=ucpar2,fam="NB2",iprint=TRUE)
cat("pos stable LT + Frank in trivariate mixture of max-id\n")
par1=c(1.7,-.2,-.2,.4,1.5,1.1)
lsfrknb1=mltscop2(y,wcb8594[,2:3],pcop3=pmxid3ls,pcop2=pmxid2ls,start=par1,
  family="NB1", iprint=TRUE,cparlb=c(0,0),cparub=30)
cat("AR(2) Gaussian/ NB1\n")
nb1param=c(1.7,-.2,-.2,.4,.5,.3)
out=nlm(nb1ar2nllk,p=nb1param,hessian=TRUE,steptol=1.e-4,yy=y,xdat=wcb8594[,2:3],
  print.level=1)
acov=solve(out$hess)
SE=sqrt(diag(acov)) 
cat("SEs=",SE,"\n")
outpred=rmspe.tvn(out$estimate,y,wcb8594[,2:3],family="NB1",iprint=TRUE)
cat("RMSPE=",outpred$rmse,"\n")
}
\keyword{copula}
\keyword{maximum likelihood}
\keyword{regression}
\keyword{ts}
