
# R interface
# nq = number of quadrature points
# Output
#   nodes: nq-vector of quadrature nodes
#   weights: nq-vector of quadrature weights
# link to C code which is translation of jacobi.f in Stroud and Secrest (1966)
gausslegendre=function(nq)
{ # need a shift because C code indexes 1 to nq, not 0 to  nq-1
  if(nq<=0 || nq>70) stop("nq between 1 and 70")
  nq1=nq+1
  out= .C("gauleg", x1=as.double(0), x2=as.double(1),
         as.integer(nq), xq=as.double(rep(0,nq1)), wq=as.double(rep(0,nq1)) )
  list(nodes=out$xq[-1], weights=out$wq[-1])
}

# this might be needed for some functions
#gldefault=gausslegendre(25)

