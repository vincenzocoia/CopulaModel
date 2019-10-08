
# generic functions in the .Rd files: 
# pcond qcond pcop dcop logdcop rcop rmcop 

pcop=function(u,v,cpar) { pfrk(u,v,cpar) } 
dcop=function(u,v,cpar) { dfrk(u,v,cpar) } 
logdcop=function(u,v,cpar) { logdfrk(u,v,cpar) } 
pcond=function(v,u,cpar) { pcondfrk(v,u,cpar) } 
qcond=function(p,u,cpar) { qcondfrk(p,u,cpar) } 
rcop=function(n,cpar,icheck=F) { rfrk(n,cpar) } 
rmcop=function(n,d,cpar) { rmfrk0(n,d,cpar) } 

# list of functions from other packages
# pmvnorm, GenzBretz  (library mvtnorm)
# garchFit, pstd (library fGarch)
# minimum.spanning.tree (library igraph0)
# permn    (library combinat)
# abind (library abind)
