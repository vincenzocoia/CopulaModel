\name{discreteresponse}
\Rdversion{1.1}
\alias{discreteresponse}
\alias{latentBVNnllk1}
\alias{latentBVNnllk2}
\alias{ieenllk}
\alias{MVNlatent1}
\alias{MVNlatent2}
\alias{vinebivEvsO1}
\title{
Copula models for count response data
}
\description{
Copula models for count response data, including discretized MVN
diagnostics from bivariate margins
}
\usage{
latentBVNnllk1(rho,param,ucdf,xdat1,xdat2,y1,y2)
latentBVNnllk2(rho,par1,par2,cdf1,cdf2,xdat1,xdat2,y1,y2)
ieenllk(param,upmf,ydat,xdat,LB,UB) # independent estimating equations
MVNlatent1(ydat,xdat,nrep,upmf,ucdf,upmfcdf,mx,ustart,LB,UB,prlevel=0)
   # discretized multivariate normal/Gaussian with univariate marginal
   # regression model being specified by upmf/ucdf and common
   # regression parameters for each margin
MVNlatent2(ydat,xdat,nrep,unllks,upmfs,ucdfs,upmfcdfs,mx,ustart,LB,UB,prlevel=0)
   # discretized multivariate normal/Gaussian with univariate marginal
   # regression models being specified by upmfs/ucdfs and different
   # regression parameters for each margin
vinebivEvsO1(ydat,xdat,nrep,mx,uparam,upmfcdf,cpar,A,pcop)
}
\arguments{
\item{rho}{latent correlation parameter}
\item{param}{common parameter vector for univariate margin of the latent 
bivariate normal model; assumed estimated in a previous step in latentBVNnllk1}
\item{par1}{univariate parameter vector for univariate margin 1 of the latent 
bivariate normal model; assumed estimated in a previous step}
\item{par2}{univariate parameter vector for univariate margin 2 of the latent 
bivariate normal model; assumed estimated in a previous step}
\item{uparam}{univariate parameter vector, common to all univariate margin}
\item{ucdf}{function for univariate cdf}
\item{upmf}{function for common univariate pmf for ieenllk}
\item{upmfcdf}{function that efficiently computes univariate pmf and cdf
up to an upper limit mx}
\item{unllks}{vector of strings for univariate negative log-likelihoods}
\item{upmfs}{vector of strings for pmf functions for margins 1,...,d}
\item{ucdfs}{vector of strings for cdf functions for margins 1,...,d}
\item{upmfcdfs}{vector of strings for pmf/cdf functions for computing up
to mx}
\item{mx}{bound used for Expected vs Observed tables in univariate/bivariate
margins}
\item{cdf1}{function for univariate cdf for margin 1}
\item{cdf2}{function for univariate cdf for margin 1}
\item{xdat1}{covariates for margin 1 (nxq), q=#covariates}
\item{xdat2}{covariates for margin 2 (nxq), q=#covariates}
\item{xdat}{covariate matrix ((n*d)xq), q=#covariates, d=#repeated measurements}
\item{y1}{count response 1, vector of length n}
\item{y2}{count response 2, vector of length n}
\item{ydat}{count response, vector of length n*d}
\item{nrep}{#repeated measurements per subject}
\item{ustart}{starting parameter point for univariate model}
\item{LB}{lower bound for param}
\item{UB}{upper bound for param}
\item{A}{dxd vine array with 1:d on diagonal}
\item{pcop}{function for copula cdf (with scalar parameter) in tree 1 of vine }
\item{cpar}{parameter values for pair-copulas in tree 1 of vine}
\item{prlevel}{print.level for nlm()}
}
\value{
negative log-likelihood value for latentBVNnllk1, latentBVNnllk2;

negative of sum of univariate log-likelihoods for ieenllk;

list(uparam vector or parmat matrix,rhvec,E1arr,O1arr,E2arr,O2arr) 
[univariate parameter estimates.
latent correlation estimates, Expected and Observed counts for
univariate/bivariate] for MVNlatent1,MVNlatent2;

list(rhosub,E1arr,O1arr,E2arr,O2arr) for
mvnbivEvsO1 where rhosub is based on subset and
rhopar is based on full data set;

list(E2arr,O2arr) for vinebivEvsO1  with Expected and Observed counts for
bivariate margins in first tree of vine.
}
\seealso{
\code{\link{rvinediscrete}}
}
\examples{
# count regression model with common betas for all margins: longitudinal count
data(rwmsubset)
rwm=rwmsubset
rwm$agec=(rwm$age-50)/10
rwm$ageq=(rwm$agec)^2
rwm$handfra=rwm$handper/100
xdat=cbind(rwm$sex,rwm$agec,rwm$ageq,rwm$hsat,rwm$handfra,rwm$univ)
xdat=as.matrix(xdat)
ydat=rwm$docvis
nc=ncol(xdat)
nc1=nc+1
out=MVNlatent1(ydat,xdat,nrep=5,upmf=nb1pmf,ucdf=nb1cdf,upmfcdf=nb1pmfcdf,
  mx=7,ustart=c(1.7,.3,.2,.1,-.2,.7,-.5,1.5),
  LB=c(rep(-20,nc1),0),UB=rep(10,nc1+1),prlevel=1)
print(out$uparam)
print(out$rhvec)
# 
# GP1 count regression models for each margin: uncommon regression coefficients 
data(kzrepmeas)
kz=kzrepmeas
kz$agehun=kz$age/100
xdat=cbind(kz$agehun,kz$sex,kz$msmok)
ydat=kz[,6]
nrep=4
outgp1=MVNlatent2(ydat,xdat,nrep,unllks=rep("gp1nllk",4),ucdfs=rep("gp1cdf",4),
  upmfcdfs=rep("gp1pmfcdf",4),
  mx=3,ustart=c(0,0,0,0,1),LB=c(-20,-20,-20,-20,0),UB=rep(10,5),prlevel=0)
print(outgp1$uparam)
print(outgp1$rhvec)
}
\keyword{copula}
\keyword{regression}
