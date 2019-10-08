\name{ordinalex}
\alias{ordinalex}
\title{Multivariate ordinal data set}
\usage{
data(ordinalex) # components b0cut, b1, d, ncateg, ncl, rmat, xx, yy
}
\description{
A simulated data set yy (200x4) with a repeated measures ordinal response 
with 4 measurements and 3 categories, 
a covariate matrix xx (200x4) with values in (-1,1). 
The values used for the simulation are:
b0cut[1]=-0.5, b0cut[2]=0.5, b1=0.4, latent AR(1) correlation=0.5
}
\format{
The following are components.
\describe{
\item{\code{b0cut}}{cutpoints for the intercept b0}
\item{\code{b1}}{slope b1}
\item{\code{d}}{dimension of repeated measures}
\item{\code{ncateg}}{number of ordinal categories}
\item{\code{ncl}}{number of clusters}
\item{\code{rmat}}{latent correlation matrix}
\item{\code{xx}}{200x4 covariate matrix with values in (-1,1)} 
\item{\code{yy}}{200x4 response variable matrix with values in 1,2,3} 
}
}
\keyword{datasets}
