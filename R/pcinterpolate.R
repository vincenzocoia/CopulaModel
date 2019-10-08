
# R interface
# Piecewise Cubic Interpolation, Hermite cubic monotone interpolant
# link to f77 code 
 
# piecewise cubic to get derivative at grid points
# x = vector of values (assumed increasing in value)
# fn = vector with same length of x, monotone in x
# Output:
#   deriv = vector of derivatives at x values (based on interpolant)
pcderiv=function(x,fn)
{ if(length(x) != length(fn)) stop("x and fn should have same length")
  #if(!is.loaded("pchez")) dyn.load("./dpch.so")
  n=length(x)
  out= .Fortran("pchez", as.integer(n), as.double(x), as.double(fn),
   deriv=as.double(rep(0,n)), as.logical(FALSE), as.double(1.), as.integer(1),
   ierr=as.integer(0) )
  if(out$ierr!=0) cat("ierr code in pchez is ", out$ierr,"\n")
  out$d
}

# piecewise cubic to interpolate at new x values 
# x = vector of values (assumed increasing in value)
# fn = vector with same length of x, monotone in x
# deriv = vector of derivatives at x values (based on interpolant)
# xnew = vector of new x values for interpolation, of length nval
# Outputs: 
#   matrix with dimension nval x 2;  2 columns fval, dval are:
#   fval = vector with interpolated f values at xnew
#   dval = vector of interpolated derivatives at xnew
#
pcinterpolate=function(x,fn,deriv,xnew)
{ if(length(x) != length(fn)) stop("x and fn should have same length")
  #if(!is.loaded("pchev")) dyn.load("./dpch.so")
  n=length(x)
  nval=length(xnew)
  out= .Fortran("pchev", as.integer(n), as.double(x), as.double(fn),
   as.double(deriv), as.integer(nval), as.double(xnew), 
   fval=as.double(rep(0,nval)), dval=as.double(rep(0,nval)), 
   ierr=as.integer(0) )
  if(out$ierr!=0) cat("ierr code in pchev is ", out$ierr,"\n")
  cbind(out$fval,out$dval)
}

# test example

#n=21
#n=31
#n=41
#n=51
#n=101
#x=seq(0,pi/2,length=n)
#fn=sin(x)

#der=pcderiv(x,fn)
#print(cbind(der,cos(x)))
#xnew=seq(.05,1.,.05)
#out=pcinterpolate(x,fn,der,xnew)
#fval=sin(xnew)
#print(cbind(out[,1],fval,abs(out[,1]-fval)))
#cat("max err in deriv : ", max(abs(der-cos(x))), "\n")
#cat("max err in interp : ", max(abs(out[,1]-fval)), "\n")
