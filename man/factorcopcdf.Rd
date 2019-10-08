\name{factorcopcdf}
\Rdversion{1.1}
\alias{factorcopcdf}
\alias{pfact1cop}
\alias{pfact2cop}
\alias{pcondfact1}
\alias{pcondfact2}
\alias{dfact1cop}
\alias{dfact2cop}
\alias{pfact1gau}
\alias{pfact1frk}
\alias{pfact1gum}
\alias{pfact1bb1}
\alias{pfact2gau}
\alias{pfact2frk}
\alias{pfact2gum}
\alias{pfact2bb1frk}
\title{
Bivariate marginal copula cdfs for 1-factor and 2-factor copula models
}
\description{
Bivariate marginal copula cdfs for 1-factor and 2-factor copula models
}
\usage{
pfact1cop(u1,u2,pcondcop,param,nq)
pcondfact1(u2,u1,pcondcop,dcop,param,nq)
dfact1cop(u1,u2,dcop,param,nq)
pfact2cop(u1,u2,pcondcop1,pcondcop2,param1,param2,nq)
pcondfact2(u2,u1,pcondcop1,pcondcop2,dcop1,dcop2,param1,param2,nq)
dfact2cop(u1,u2,pcondcop1,dcop1,dcop2,param1,param2,nq)
pfact1frk(u1,u2,param)  # nq defaulted to 35 etc
pfact2frk(u1,u2,pmatrix)  # nq defaulted to 35 etc
}
\arguments{
\item{u1}{vector of values in interval 0,1; }
\item{u2}{vector of values in interval 0,1; same length as u1 }
\item{param}{vector of length 2 or 2xq matrix where q is number of parameters} 
\item{param1}{vector of length 2 or 2xq matrix where q is number of parameters 
for the bivariate copula (e.g. BB1) in pcondcop; 
parameters that link observed variables to latent factor 1 }
\item{param2}{vector of length 2,
parameters that link observed variables to latent factor 2 }
\item{pcondcop}{function for conditional cdf of linking copula for 1-factor model}
\item{pcondcop1}{function for conditional cdf of linking copula for factor 1
for 2-factor model}
\item{pcondcop2}{function for conditional cdf of linking copula for factor 1
for 2-factor model}
\item{dcop}{function for pdf of linking copula for 1-factor model}
\item{dcop1}{function for pdf of linking copula for factor 1 for 2-factor model}
\item{dcop2}{function for pdf of linking copula for factor 2}
\item{nq}{number of quadrature points for Gauss-Legendre quadrature }
\item{pmatrix}{For most choices,
first column has parameters that link observed variables to latent factor 1,
second column has parameters that link observed variables to latent factor 2.
For pfact2bb1frk: first two columns for BB1 parameters, third column
for Frank for latent factor 2.  }
}
\value{
cdf or conditional cdf or pdf value(s)
}
\seealso{
\code{\link{factorcopsim}}
}
\examples{
th1a=frk.b2cpar(.7)
th1b=frk.b2cpar(.6)
th2a=frk.b2cpar(.5)
th2b=frk.b2cpar(.4)
u1=seq(.1,.9,.2)
u2=seq(.3,.7,.1)
pfact1frk(u1,u2,c(th1a,th1b))
pfact2frk(u1,u2,matrix(c(th1a,th1b,th2a,th2b),2,2))
}
\keyword{copula}
\keyword{factor model}
