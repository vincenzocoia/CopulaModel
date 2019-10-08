\name{taucor}
\Rdversion{1.1}
\alias{taucor}
\alias{taub}
\title{
Kendall's tau for two variables
}
\description{
Data version of Kendall's tau for two variables, using Knight's algorithm
}
\usage{
taucor(x,y=0,pflag=0)
taub(otab) # Kendall's taub, adjusted for ties for 2-way table
}
\arguments{
\item{x}{first variable, ties allowed; n x d data matrix if y=0}
\item{y}{second variable, ties allowed}
\item{pflag}{print flag,
0 = silent,
1 = print number of ties, exchange count and other summeries
2 = print sorted data
3 = print sorting iterations}
\item{otab}{ordinal 2-way table/matrix}
}
\value{
\item{tau}{Kendall's tau with value in (-1,1); denominator adjusts for ties
in case of ordinal variables}
}
\references{
Knight WR (1966). A computer method for calculating Kendall's tau with
ungrouped data. Journal of the American Statistical Association,
61, 436-439.
}
\examples{
set.seed(1234)
x=rnorm(50)
y=.6*x+.8*rnorm(50)
tau=taucor(x,y)
tau2=cor(x,y,method="kendall")  # this is too slow when length(x) is large
cat(tau,tau2,"\n")
z=.5*x+.5*y+sqrt(.5)*rnorm(50)
tau=taucor(cbind(x,y,z))
tau2=cor(cbind(x,y,z),method="kendall")
print(max(abs(tau-tau2)))
# check for two ordinal variables
n=25
set.seed(123)
x1=floor(runif(n,0,3))
x2=floor(runif(n,0,3))
otab=table(x1,x2)
taub(otab) # should be same as next two lines
taucor(x1,x2)
cor(x1,x2,method="kendall")
}
\keyword{dependence}
