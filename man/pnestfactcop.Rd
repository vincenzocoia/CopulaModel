\name{pnestfactcop}
\Rdversion{1.1}
\alias{pnestfactcop}
\alias{pnest2cop}
\alias{pnest2frk}
\alias{pnest2gum}
\alias{pnest2t}
\alias{pnest2tgum}
\alias{pnest2tbb1}
\alias{pnest2gumbb1}
\title{
Bivariate marginal copula cdfs for nested-factor copula models
}
\description{
Bivariate marginal copula cdfs for nested-factor copula models
}
\usage{
pnest2cop(u1,u2,dcop1,pcondcop2,param1,param2,nq)
pnest2frk(u1,u2,param)  # nq defaulted to 35 etc
pnest2gum(u1,u2,param)  # nq defaulted to 35 etc
pnest2t(u1,u2,param,df)  # nq defaulted to 35 etc
pnest2tgum(u1,u2,param,df)  # nq defaulted to 35 etc
pnest2tbb1(u1,u2,param,df)  # nq defaulted to 35 etc
pnest2gumbb1(u1,u2,param)  # nq defaulted to 35 etc
}
\arguments{
\item{u1}{vector of values in interval 0,1; }
\item{u2}{vector of values in interval 0,1; same length as u1 }
\item{param1}{vector of length 2 or 2xq matrix where q is number of parameters 
for the bivariate copula (e.g. BB1) in dcop1; 
parameters that link observed variables to common latent}
\item{param2}{vector of length 2 or 2xq matrix where q is number of parameters
for the bivariate copula in pcondcop2;
parameters that link observed variables to nested group latent variable}
\item{param}{column 1 has parameters for global/common latent, column 2 
(and column 3 for pnest2tbb1,pnest2gumbb1) has parameters for group latent}  
\item{df}{shape or df parameter for bivariate t linking copula}
\item{dcop1}{function for pdf of copula family for global/common latent} 
\item{pcondcop2}{function for conditional cdf of copula family for nested
group factor }
\item{nq}{number of quadrature points for Gauss-Legendre quadrature }
}
\seealso{
\code{\link{structcop}}
}
\details{
This function is the bivariate marginal cdf of nested copula for two 
variables in different groups (for within the same group, the
bivariate maeginal cdf is a 1-factor copula margin.
}
\value{
cdf value(s)
}
\examples{
pnest2fgm=function(u1,u2,param)
{ f=pnest2cop(u1,u2,dfgm,pcondfgm,param[,1],param[,2],35);
  f=f*(f<=1)+(f>1)
  f
}
th1a=frk.b2cpar(.7)
th1b=frk.b2cpar(.6)
th2a=frk.b2cpar(.5)
th2b=frk.b2cpar(.4)
u1=seq(.1,.9,.2)
u2=seq(.3,.7,.1)
pnest2frk(u1,u2,matrix(c(th1a,th1b,th2a,th2b),2,2))
pnest2gum(u1,u2,matrix(c(th1a,th1b,th2a,th2b),2,2))
pnest2t(u1,u2,matrix(c(.5,.6,.5,.4),2,2),c(5,5))
pnest2tgum(u1,u2,matrix(c(.5,.6,1.5,1.4),2,2),5)
pnest2tbb1(u1,u2,matrix(c(.5,.6,.4,.5,1.5,1.4),2,3),5)
pnest2gumbb1(u1,u2,matrix(c(1.5,1.6,.4,.5,1.5,1.4),2,3))
pnest2fgm(u1,u2,matrix(c(.5,.4,.6,.7),2,2))
}
\keyword{copula}
