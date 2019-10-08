# interface to plot contour of bivariate copula density with N(0,1) margins

# cpar = copula parameter
# zvec = grid points of N(0,1) density for computing bivariate density
# dcop = function name of bivariate copula density with arguments u,v,cpar
# irefl = T to take reflection of dcop
# Output: nothing is returned but a contour plot is displayed
contourBivCop=function(cpar,zvec,dcop,irefl=F)
{ f=dnorm(zvec)
  cdf=pnorm(zvec)
  nn=length(zvec)
  dens=matrix(0,nn,nn)
  for(i1 in 1:nn)
  { for(i2 in 1:nn)
    { dens[i1,i2]=dcop(cdf[i1],cdf[i2],cpar)*f[i1]*f[i2] }
  }
  if(irefl) contour(zvec,zvec,dens[nn:1,nn:1])
  else contour(zvec,zvec,dens)
  invisible(dens)
}

