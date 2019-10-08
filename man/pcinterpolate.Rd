\name{pcinterpolate}
\Rdversion{1.1}
\alias{pcinterpolate}
\alias{pcderiv}
\title{
Interpolate a monotone function (via piecewise cubic Hermite)
}
\description{
Interpolate a monotone function from a two-column table
}
\usage{
pcderiv(x,fn)
pcinterpolate(x,fn,deriv,xnew)
}
\arguments{
\item{x}{vector of x values }
\item{fn}{function values corresponding to each value in x} 
\item{deriv}{vector of (estimated) derivative values assuming monotonicity,
output of pcderiv and input to pcinterpolate }
\item{xnew}{vector of new x values to interpolate values of the function. }
}
\details{The monotone piecewise cubic interpolation algorithm is from
the references given below.
}
\references{
Fritsch FN and Carlson RE (1980), Monotone piecewise
cubic interpolation. Siam J Numerical Analysis, 17, 238-246.

Kahaner D, Moler CB and Nash S (1989). Numerical Methods and Software
Prentice Hall.
}
\value{
\item{deriv}{for pcderiv, estimated derivatives at x values }
\item{y}{for pcinterpolate, length(xnew) x 2  matrix of interpolated fn values 
and derivatives}
}

\seealso{
\code{\link{makedeptable}}
}
\examples{
n=21
x=seq(0,pi/2,length=n)
fn=sin(x)
der=pcderiv(x,fn)
print(cbind(der,cos(x)))
xnew=seq(.05,1.,.05)
out=pcinterpolate(x,fn,der,xnew)
fval=sin(xnew)
print(cbind(out[,1],fval,abs(out[,1]-fval)))
cat("max err in deriv : ", max(abs(der-cos(x))), "\n")
cat("max err in interp : ", max(abs(out[,1]-fval)), "\n")
}
\keyword{interpolation}
