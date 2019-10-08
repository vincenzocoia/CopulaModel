\name{bifct}
\Rdversion{1.1}
\alias{bifct}
\alias{bifct2}
\alias{trifct}
\alias{subgr.consistent}
\title{
bi-factor and tri-factor structure 
}
\description{
bi-factor and tri-factor structure: matrix inverse and determinant 
}
\usage{
bifct(grsize,rh1,rh2)
bifct2(grsize,rh1,rh2) # version 2, with inverse and det of a smaller matrix
trifct(grsize,sbgrsize,rh1,rh2,rh3)
subgr.consistent(grsize,sbgrsize) # check that grsize,sbgrsize are consistent
}
\arguments{
\item{grsize}{vector of group sizes for mgrp groups with sum(grsize)=d} 
\item{sbgrsize}{vector of subgroup sizes by partitioning grsize vector}
\item{rh1}{vector of correlations with global latent variable}
\item{rh2}{vector of partial correlations with group latent variable given 
global; ith variable linked to group g}
\item{rh3}{vector of partial correlations with subgroup latent variable given
global and group latent; ith variable linked to group sg}
}
\value{
$fctmat (correlation matrix), $fctdet (determinant), 
$fctinv (inverse) for bifct, bifct2, trifct.

TRUE or FALSE for subgr.consistent depending on consistency
}
\seealso{
\code{\link{mvtfact}}
}
\examples{ 
d=29
# bi-factor
grsize=c(5,5,3,3,4,2,2,2,3) # sum(grsize)=d
bifpar=c(0.84,0.63,0.58,0.78,0.79,  
 0.87,0.80,0.74,0.71,0.57,   0.83,0.77,0.80, 0.74,0.78,0.71, 
 0.71,0.58,0.63,0.64, 0.77,0.73, 0.72,0.72, 0.45,0.49, 0.66,0.70,0.69, 
 0.67,0.58,0.15,0.70,0.47, 0.32,0.27,0.73,0.19,0.12, 0.35,0.23,0.53,
 0.89,0.63,0.63, 0.28,0.51,0.80,-0.02, 0.37,0.37, 0.82,0.82, 
 0.70,0.70,-0.25,0.55,0.10)
bifobj=bifct(grsize,bifpar[1:d],bifpar[(d+1):(2*d)])
rmat=bifobj$fctmat
print(det(rmat)-bifobj$fctdet)
print(max(abs(solve(rmat)-bifobj$fctinv)))
# tri-factor
grsize=c(11,18) # sum(grsize)=d
sbgrsize = c(5,4,2,5,3,3,2,2,3) # sum(sbgrsize)=d
subgr.consistent(grsize,sbgrsize) 
trifpar=c(0.86,0.77,0.73,0.71,0.54,0.70,
 0.55, 0.60,  0.64, 0.40, 0.44, 0.84, 0.65, 0.58,  0.78, 0.78, 0.83, 0.79,
 0.81, 0.77,  0.80, 0.73, 0.76, 0.73, 0.70, 0.70,  0.62, 0.72, 0.68, 0.25,
 0.37, 0.46,  0.23, 0.42, 0.28, 0.43, 0.45,-0.05,  0.61, 0.62, 0.54, 0.09,
 0.09, 0.34,  0.45, 0.07,-0.06,-0.02,-0.30,-0.12, -0.27, 0.13, 0.06, 0.56,
 0.50, 0.51, -0.31, 0.03, 0.81, 0.18, 0.21,-0.02, -0.16, 0.19, 0.37, 0.81,
 0.05, 0.48,  0.48, 0.50, 0.66, 0.10, 0.66, 0.29,  0.30, 0.11, 0.46, 0.86,
 0.60, 0.56,  0.37, 0.37, 0.76, 0.76, 0.14, 0.31,  0.18)
trifobj=trifct(grsize,sbgrsize,trifpar[1:d],trifpar[(d+1):(2*d)],
   trifpar[(2*d+1):(3*d)])
rmat=trifobj$fctmat
print(det(rmat)-trifobj$fctdet)
print(max(abs(solve(rmat)-trifobj$fctinv)))
}
\keyword{factor model}
