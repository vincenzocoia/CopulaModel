
# independence, comonotonicity, countercomonotonicity copulas

# u = vector or scalar with values in (0,1)
# v = value in (0,1) if u is scalar
#   for dimensions d>=3, use v=-1 default and input a vector
# cpar is not used, but is listed as an argument to match generic form
pindepcop=function(u,v=-1,cpar=0)
{ if(v[1]<0) return(prod(u)) else return(u*v)
}

# u = vector or scalar
# v = value in (0,1) if u is scalar
#   for dimensions d>=3, use v=-1 default and input a vector
# cpar is not used, but is listed as an argument to match generic form
pcomonocop=function(u,v=-1,cpar=0)
{ if(v[1]<0) return(min(u)) else pmin(u,v)
}

# u = 2-vector or scalar
# v = value in (0,1) if u is scalar
# cpar is not used, but is listed as an argument to match generic form
pcountermono=function(u,v=-1,cpar=0)
{ if(v[1]<0) return(max(u[1]+u[2]-1,0)) else pmax(u+v-1,0) 
}


# conditional cdf for C^+, C^\perp, C^-
# for use in ptree2cop
# see above documentation for arguments of functions

pcondindep=function(v,u,cpar=0) { v }

pcondcomono=function(v,u,cpar=0) 
{ ii=(u<=v) 
  as.numeric(ii)
}

pcondcountermono=function(v,u,cpar=0)
{ ii=(v>=1-u) 
  as.numeric(ii)
}

