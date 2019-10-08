\name{ordinal.bivcop}
\Rdversion{1.1}
\alias{ordinal.bivcop}
\alias{bivcopOrdinalnllk}
\title{
Negative log-likelihood for bivariate marginal copula model for discrete variables
}
\description{
Negative log-likelihood for bivariate marginal copula model for discrete 
variables; this function is a counterpart of polychoric() when
the bivariate Gaussian copula is replaced by another bivariate copula.
}
\usage{
ordinal.bivcop(odatfr,ucuts,pcop,cparstart,LB=0,UB=10,iprint=F,prlevel=0)
bivcopOrdinalnllk(cpar,ucuts,bfr,jj1,jj2,pcop,LB=0,UB=10)
}
\arguments{
\item{odatfr}{nrecx(d+1) matrix: d columns of ordinal responses and final
column with frequency of each distinct observed d-vector (it could be
just a vector of 1s)}
\item{ucuts}{cutpoints in U(0,1) scale as (ncateg+1)xd matrix, obtained via 
unifcuts; so ncateg=nrow(ucuts)-1 and d=ncol(ucuts)}
\item{pcop}{function with pair-copula cdf}
\item{cpar}{copula parameter for pcop}
\item{cparstart}{vector of starting points, dimension is d*(d-1)/2
times the dimension of the parameter for pcop }
\item{bfr}{vector of bivariate frequencies}
\item{jj1}{index of first variable when enumerating through pairs}
\item{jj2}{index of second variable}
\item{LB}{lower bound on parameter of pcop}
\item{UB}{upper bound on parameter of pcop}
\item{iprint}{flag for printing of intermediate results}
\item{prlevel}{print.level for nlm for numerical optimization}

}
\value{
List with

$nllkvec = vector of length d*(d-1)/2 with negative log-likelihoods
at the bivariate MLEs for each pair

$cparvec = vector of  length d*(d-1)/2 times the dimension of the
pcop parameter, with the bivariate MLEs for each pair

$summary = ncateg x (2*ncateg) x d*(d-1)/2 array with observed and expected
bivariate marginal counts.

Compare output of function polychoric() with the bivariate normal/Gaussian
copula as a latent model.
}
\seealso{
\code{\link{bivcopnllk}}
\code{\link{ordinal}}
}
\examples{
# convert bivariate t with fixed degree of freedom to use with above
pbvtcop5=function(u,v,rho) { pbvtcop(u,v,c(rho,5)) }
pbvtcop10=function(u,v,rho) { pbvtcop(u,v,c(rho,10)) }
data(ltmconv) # to use data set 'env'
d=3
nc=4  # ncateg
ucuts=unifcuts(env[,1:d])
bd=pnorm(-6)
ucuts=rbind(rep(bd,d),ucuts,rep(1-bd,d))
zcuts=qnorm(ucuts)
odatfr=ordinal2fr(env[,1:d],nc)
polyr=polychoric(odatfr,zcuts,iprint=FALSE,prlevel=0)
rhst=cormat2vec(polyr$polych)
rhstgal=rhst
cparst=depmeas2cpar(rhstgal,"rhoN","galambos")
cat("\nGalambos\n")
outgal=ordinal.bivcop(odatfr,ucuts,pcop=pgal,cparstart=cparst,iprint=TRUE,
   prlevel=0,LB=0,UB=10)
cat("\nt(5)\n")
outbvt5=ordinal.bivcop(odatfr,ucuts,pcop=pbvtcop5,cparstart=rhst,iprint=TRUE,
   prlevel=0,LB=-1,UB=1)
cat("\nt(10)\n")
outbvt10=ordinal.bivcop(odatfr,ucuts,pcop=pbvtcop10,cparstart=rhst,iprint=TRUE,
   prlevel=0,LB=-1,UB=1)
cat("\nBB1\n")
taust=bvn.cpar2tau(rhst)
dd=d*(d-1)/2
cparst=NULL
for(ii in 1:dd)
{ cpar=bb1.tau2eqlm(taust[ii])
  cparst=c(cparst,cpar[1:2])
}
outbb1=ordinal.bivcop(odatfr,ucuts,pcop=pbb1,cparstart=cparst,iprint=TRUE,
   prlevel=0,LB=rep(c(0,1),dd),UB=5)
cat("\nnllk comparison\n")
print(cbind(outgal$nllkvec,outbvt5$nllkvec,outbvt10$nllkvec,outbb1$nllkvec))
}
\keyword{copula}
\keyword{maximum likelihood}
