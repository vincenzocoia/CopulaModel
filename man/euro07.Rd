\name{euro07}
\alias{euro07}
\alias{euro07lr}
\alias{euro07gf}
\alias{euro07names}
\title{log returns and GARCH-filtered log returns for some Euro markets 2007}
\usage{
data(euro07lr) # objects euro07names and euro07lr
data(euro07gf) # objects euro07names and euro07gf
}
\description{
Log returns and GARCH filtered values of log returns for 
OSEAX FTSE AEX FCHI SSMI GDAXI ATX 
(market indexes in Norway, UK, Holland, France, Switzerland, Germany,
Austria). This is a small data set with n=239 that can be used
for illustration of functions for fitting vine and factor copulas. The
original source is http://quote.yahoo.com

euro07lr has two objects: (i) euro07names has the above labels
for the markets and (ii) euro07lr is a 239x7 matrix of log returns in
2007, one column for each market (the markets were merged to common
dates) before returns were obtained

euro07gf has two objects: (i) euro07names has the above labels
for the markets and (ii) euro07lr is a list with several matrices
given below.
}

\format{
The following are components.
\describe{
\item{filter}{239x7 matrix of GARCH filtered returns}

\item{uscore}{239x7 matrix of empirical U(0,1) scores of GARCH filtered returns}

\item{zscore}{239x7 matrix of empirical normal scores of GARCH filtered returns}

\item{uscmodel}{239x7 matrix of U(0,1) scores of GARCH filtered returns
based on assuming standardized Student t distributions for the innovations}

\item{zscmodel}{239x7 matrix of normal scores of GARCH filtered returns
based on assuming standardized Student t distributions for the innovations}

\item{sigmat}{239x7 matrix of volatilities}

\item{coef}{5x7 matrix of GARCH parameter estimates}
rows are  mu, omega, alpha1, beta1, shape,
where 'shape' is the shape or degree of freedom parameter for the 
Student t innovations.}
}
\keyword{datasets}
