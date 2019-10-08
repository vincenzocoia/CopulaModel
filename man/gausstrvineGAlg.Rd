\name{gausstrvineGAlg}
\Rdversion{1.1}
\alias{gausstrvineGAlg}
\alias{gausstrvine.galg}
\title{
Gaussian truncated d-dimensional vines based on minimum 
spanning trees combined with a genetic algorithm
}
\description{
Gaussian truncated d-dimensional vines based on a genetic algorithm
and spanning trees with weights of one minus squared partial correlation
}
\usage{
gausstrvine.galg(cormat,target=0.99,nneigh=10,nbmst=10,
  selectcrit="prop",selectfact=0.1,elite=0.05,method="nfi",n=0,iprint=T)
}
\arguments{
\item{cormat}{dxd correlation matrix} 
\item{target}{target fit value (near but less than 1; 1-alpha in
publication)}
\item{nneigh}{number of nearest neighbors of minimum spanning tree}
\item{nbmst}{number of best minimum spanning trees}
\item{selectcrit}{selection criterion ("prop" for proportional or "rank")}
\item{selectfact}{selection factor with default 0.1
(total fraction of solutions retained after elitism and random selection)}
\item{elite}{percentage for elitism (default 0.05)}
\item{method}{fit index ("nfi", "cfi" or "ifi");
NFI=normed fit index, CFI=comparative fit index, IFI=incremental fit index}
\item{n}{sample size for the correlation matrix, default 0 if unknown}
\item{ntrunc}{specified upper bound in truncation level to consider}
\item{iprint}{print flag for intermediate results}
}
\value{
\item{RVM}{object with $RVM$VineA = d-dimensional vine array,
$RVM$pc = partial correlations by tree,
$RVM$Matrix = vine array in VineCopula format [d:1,d:1],
$RVM$Cor = partial correlations in VineCopula format [d:1,d:1]
}
\item{mst}{spanning trees 1,2,...d-1: $mst[[1]], $mst[[2]], ...}
\item{treeweight}{vector of length d-1 with
sum_edge log(1-rho[edge]^2) for trees 1,...d-1}
\item{trunclevel}{same as inputted ntrunc}
\item{truncval}{sum treeweight[1:trunclevel] / sum treeweight[1:(d-1)] }
}
\seealso{
\code{\link{gausstrvine}}
\code{\link{gausstrvineMST}}
}
\details{
This function depends on the minimum spanning tree algorithm in the
library igraph.
}
\references{
Brechmann and Joe (2015), Truncation of vine copulas using fit indices.
J Multivariate Analysis, v 138, 19-33.
}
\examples{
\dontrun{
library(igraph) # version 0.7.1 and 1.0.0 work
rrvec=c(
-0.04,-0.08, 0.47,-0.07, 0.67, 0.50,-0.02, 0.52, 0.61, 0.46, 0.00,-0.06,
-0.08,-0.07,-0.10, 0.15,-0.03,-0.14,-0.15,-0.06, 0.19, 0.03,-0.01,-0.04,
-0.03,-0.03, 0.10, 0.22, 0.19, 0.10,-0.09,-0.11, 0.02,-0.02, 0.62, 0.18, 
 0.10, 0.04, 0.02,-0.03, 0.05, 0.07, 0.20, 0.27, 0.16, 0.10, 0.18,-0.09,
-0.10, 0.01, 0.00, 0.32, 0.04, 0.53, 0.05,-0.16, 0.03, 0.12, 0.15, 0.06,
-0.25,-0.83,-0.24,-0.52,-0.23,-0.26,-0.11,-0.49, 0.05, 0.01,-0.14,-0.03,
-0.38,-0.11,-0.57,-0.18,-0.43, 0.35, 0.15,-0.02,-0.14,-0.13,-0.05, 0.17, 
 0.95, 0.21, 0.63, 0.18, 0.33,-0.81,-0.38,-0.05, 0.07, 0.13, 0.10, 0.08,
-0.33,-0.32,-0.15,-0.13,-0.21,-0.06, 0.42, 0.13,-0.31, 0.40,-0.03,-0.11,
-0.08,-0.05,-0.07, 0.13, 0.01, 0.25, 0.06, 0.22,-0.12,-0.17, 0.15, 0.02,
 0.66, 0.00,-0.11,-0.06,-0.02,-0.02, 0.18, 0.04, 0.26, 0.11, 0.16,-0.18,
-0.18, 0.19,-0.03, 0.52, 0.29,-0.01,-0.02,-0.04,-0.01,-0.02, 0.05, 0.01,
 0.09, 0.07, 0.04,-0.05,-0.05, 0.05, 0.02, 0.15, 0.30, 0.35,-0.22,-0.12,
-0.17,-0.09,-0.07, 0.10,-0.02, 0.15, 0.03, 0.14,-0.10, 0.01, 0.12, 0.02,
 0.80, 0.50, 0.19)

rmat=corvec2mat(rrvec)
nn=1187; d=ncol(rmat)
mxtrunc=10
bestmst=gausstrvine.mst(rmat,ntrunc=mxtrunc,iprint=F)
set.seed(123)
kbmst=gausstrvine.galg(rmat,selectcrit="prop",method="nfi",n=nn,iprint=F)
# there is a random component in the genetic algorithm
# so different seeds can lead to different suboptimal truncated vines
set.seed(12345)
kbmst2=gausstrvine.galg(rmat,selectcrit="prop",method="nfi",n=nn,iprint=F)
# Use the different vine arrays for further analysis such as
#  replace bivariate Gaussian copulas on edges with tail asymmetric copulas

# greedy algorithm : sequential MST
cat("\nsequential MST\n")
perm=diag(kbmst$RVM$VineA)
Aori=kbmst$RVM$VineA
iperm=order(perm)
AA=varrayperm(Aori,iperm)
pcmat=kbmst$RVM$pc
rperm=rmat[perm,perm]
npar=0
for(ell in 1:mxtrunc)
{ Rtrun=pcor2cor.truncvine(pcmat,AA,ntrunc=ell)
  npar=npar+(d-ell)
  if(ell>1) Rtrun=Rtrun$rmat
  cat("\nMST",ell,"-truncated R-vine with npar=", npar,"\n")
  cat("max abs difference =", max(abs(rperm-Rtrun)),"\n")
  cat("avg abs difference =", sum(abs(rperm-Rtrun))/d/(d-1),"\n")
  outtr=corDis(Rtrun,rperm,nn,npar=npar)
  cat("corDis:Dfit, 2*nllk, AIC, BIC=", outtr,"\n")
}
cat("\n============================================================\n")

# non-greedy MST with genetic algorithm
cat("\nnon-greedy, seed=123\n")
perm=diag(kbmst$RVM$VineA)
Aori=kbmst$RVM$VineA
iperm=order(perm)
AA=varrayperm(Aori,iperm)
pcmat=kbmst$RVM$pc
rperm=rmat[perm,perm]
npar=0
for(ell in 1:mxtrunc)
{ Rtrun=pcor2cor.truncvine(pcmat,AA,ntrunc=ell)
  npar=npar+(d-ell)
  if(ell>1) Rtrun=Rtrun$rmat
  cat("\nGenetic-alg with MST",ell,"-truncated R-vine with npar=", npar,"\n")
  cat("max abs difference =", max(abs(rperm-Rtrun)),"\n")
  cat("avg abs difference =", sum(abs(rperm-Rtrun))/d/(d-1),"\n")
  outtr=corDis(Rtrun,rperm,nn,npar=npar)
  cat("corDis:Dfit, 2*nllk, AIC, BIC=", outtr,"\n")
}
}
}
\keyword{partial correlation}
\keyword{vine}
