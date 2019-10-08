\name{gausstrvineMST}
\Rdversion{1.1}
\alias{gausstrvineMST}
\alias{gausstrvine.mst}
\title{
Gaussian truncated d-dimensional vines based on sequential minimum 
spanning trees
}
\description{
Gaussian truncated d-dimensional vines based on sequential minimum 
spanning trees with weights of one minus squared partial correlation
}
\usage{
gausstrvine.mst(rmat,ntrunc,iprint=F)
}
\arguments{
\item{rmat}{dxd correlation matrix} 
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
}
\details{
This function depends on the minimum spanning tree algorithm in the
library igraph0.
}
\examples{
\dontrun{
d=5
library(igraph0) # version 0.5.6 works
rmat=matrix(c(
1.00,0.76,0.76,0.74,0.67,
0.76,1.00,0.91,0.93,0.86,
0.76,0.91,1.00,0.94,0.85,
0.74,0.93,0.94,1.00,0.88,
0.67,0.86,0.85,0.88,1.00), d,d)
colnames(rmat) = rownames(rmat) = paste("V",1:d,sep="")
out=gausstrvine.mst(rmat,ntrunc=3,iprint=TRUE)
print(out) }
}
\keyword{partial correlation}
\keyword{vine}
