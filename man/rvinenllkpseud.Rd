\name{rvinenllkpseud}
\Rdversion{1.1}
\alias{rvinenllkpseud}
\alias{rvinenllkpseud2}
\title{
Pseudo observations after fitting truncated R-vine copulas 
}
\description{
Pseudo observations after fitting truncated R-vine copulas 
}
\usage{
rvinenllkpseud(parvec,udat,A,logdcopnames,pcondnames,np)
   # common pair-copula family for each tree
rvinenllkpseud2(parvec,udat,A,ntrunc,logdcopmat,pcondmat,np)
   # can be different pair-copula family for each edge of vine
}
\arguments{
\item{parvec}{parameter vector for the model}
\item{udat}{nxd matrix of uniform scores for rvinenllk functions}
\item{A}{dxd vine array with 1:d on diagonal}
\item{ntrunc}{truncation level between 1 and d-1}
\item{logdcopnames}{string vector with names of log copula pdfs of length 
ntrunc, ntrunc=truncation level}
\item{pcondnames}{string vector with names of copula conditional cdfs of length
ntrunc, ntrunc=truncation level}
\item{logdcopmat}{matrix of names of log copula pdfs for trees 1,...,ntrunc} 
\item{pcondmat}{matrix of names of conditional cdfs for trees 1,...,ntrunc}
\item{np}{dxd where np[ell,j] is size for parameter[ell,j]
for bivariate copula in tree ell, linking variables j and A[ell,j]}
}
\value{
\item{nllk}{negative log-likelihood}
\item{condforw}{nx(d-ntrunc) matrix with \eqn{C_{j|a_{ntrunc,j};S}}
in the forward direction}
\item{condbackw}{nx(d-ntrunc) matrix with \eqn{C_{a_{ntrunc,j}|j;S}} in the 
backward direction}
}
\seealso{
\code{\link{rvinenllk}}
}
\examples{
parmat=matrix(c(0,1.5,2,2.5,2.2,0,0,.5,.4,.6,0,0,0,.2,.2,0,0,0,.2,.2,0,0,0,0,.2),5,5,byrow=TRUE)
parvec=c(parmat[1,2:5],parmat[2,3:5],parmat[3,4:5],parmat[4,5])
pcondnames=c("pcondgum","pcondbvncop","pcondbvncop","pcondbvncop")
qcondnames=c("qcondgum","qcondbvncop","qcondbvncop","qcondbvncop")
np=matrix(1,5,5)
C5=Cvinearray(5)
set.seed(123)
udat=rvinesimvec(300,C5,parvec,np,qcondnames,pcondnames)
# fit 1-truncated
logdcopnames1="logdgum"
mle1=nlm(rvinenllk1.trunc,p=parvec[1:4],
    udat=udat,A=C5,logdcopnames=logdcopnames1,pcondnames=pcondnames[1],
    hessian=TRUE,iterlim=30,print.level=1,LB=1,UB=20)
pseud1=rvinenllkpseud(mle1$estimate,udat,C5,logdcopnames1,pcondnames[1],np)
# should be 2|1 3|1 4|1 5|1 in $condforw
zdat1=nscore(pseud1$condforw)
out=semicortable(zdat1,inscore=TRUE)
out[,1]=out[,1]+1
out[,2]=out[,2]+1
print(out)
# fit 2-truncated
logdcopnames2=c("logdgum","logdbvncop")
mle2=nlm(rvinenllk1.trunc,p=parvec[1:7],
    udat=udat,A=C5,logdcopnames=logdcopnames2,pcondnames=pcondnames[1:2],
    hessian=TRUE,iterlim=30,print.level=1,LB=c(1,1,1,1,-1,-1,-1),
    UB=c(20,20,20,20,1,1,1))
pseud2=rvinenllkpseud(mle2$estimate,udat,C5,logdcopnames2,pcondnames[1:2],np)
# should be 3|12 4|12 5|12 in $condforw
zdat2=nscore(pseud2$condforw)
out=semicortable(zdat2,inscore=TRUE)
out[,1]=out[,1]+2
out[,2]=out[,2]+2
print(out)
}
\keyword{maximum likelihood}
\keyword{vine}
