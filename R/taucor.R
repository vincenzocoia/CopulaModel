# Function for Kendall's tau, computed with O(n log(n)) algorithm;
# also function for Kendall taub with input of 2-way table.

# R interface (variable names as in Algorithm 6.1 of Joe (2014)).
# input to C function
#   x, y  : vectors of data
#   pflag : 0 = silent,
#           1 = print no. ties, exchange count and other summeries
#           2 = print sorted data
#           3 = print sorting iterations
# output to C function
#   x, y  : sorted data
#   n     : sample size (no. of pairs)
#   tau   : Kendall's tau
#   numer : numerator in Kendall's tau,
#             0.5*n*(n-1) - 2*(exchange count) when no ties
#   den  : in Kendall's tau
#   tx   : no. pairs of tied x
#   ty    : no. pairs of tied y
#   txy   : no. pairs of tied (x, y)
#   pflag : print flag

# Knight's method, JASA 1966, v 61, pp 436-439.
# x = n-vector or matrix with n columns
# y = n-vector if x is an n-vector, otherwise set to 0
# pflag = print flag for intermediate calculations
# Output: Kendall tau value or matrix of Kendall tau value as in 
#                              cor(x, method="kendall")
taucor=function(x, y=0, pflag=0)
{ if(is.data.frame(x)) x=as.matrix(x)
  if(is.matrix(x))
  { d=ncol(x); n=nrow(x)
    taumat=matrix(1,d,d)
    pflag=0
    for(j1 in 1:(d-1))
    { for(j2 in (j1+1):d)
      { out= .C("ktau", as.double(x[,j1]), as.double(x[,j2]),
            n=as.integer(n), tau=as.double(0), numer=as.double(0),
            den=as.double(0), tx=as.integer(0), ty=as.integer(0),
            txy=as.integer(0), pflag=as.integer(pflag))
        taumat[j1,j2]=out$tau; taumat[j2,j1]=out$tau
      }
    }
    return(taumat)
  }
  else 
  { # x,y are vectors (no checks)
    m=length(x)
    n=length(y)
    if(m==0 || n==0) stop("both 'x' and 'y' must be non-empty")
    if(m!= n) stop("'x' and 'y' must have the same length")
    if(n>50) pflag = min(pflag, 1)
    out= .C("ktau", x=as.double(x), y=as.double(y),
            n=as.integer(n), tau=as.double(0), numer=as.double(0),
            den=as.double(0), tx=as.integer(0), ty=as.integer(0),
            txy=as.integer(0), pflag=as.integer(pflag))
    return(out$tau)
  }
}

# otab = ordinal 2-way table 
# Output: Kendall taub value
taub=function(otab)
{ nr=nrow(otab)
  nc=ncol(otab)
  n=sum(otab)
  nconc=0
  for(j in 1:(nr-1))
  { for(k in 1:(nc-1))
    { for(j2 in (j+1):nr)
      { for(k2 in (k+1):nc)
        { nconc=nconc+otab[j,k]*otab[j2,k2] }
      }
    }
  }
  ndisc=0
  for(j in 1:(nr-1))
  { for(k in 2:nc)
    { for(j2 in (j+1):nr)
      { for(k2 in 1:(k-1))
        { ndisc=ndisc+otab[j,k]*otab[j2,k2] }
      }
    }
  }
  #cat(nconc,ndisc,nconc-ndisc,"\n")
  nn=n*(n-1)/2
  # ties on variables 1 and 2
  marg1=apply(otab,1,sum)
  marg2=apply(otab,2,sum)
  ties1=sum(marg1*(marg1-1))/2
  ties2=sum(marg2*(marg2-1))/2
  den1=nn-ties1
  den2=nn-ties2
  #cat(nn,ties1,ties2,sqrt(den1*den2),"\n")
  (nconc-ndisc)/sqrt(den1*den2)
}



# R interface (variable names as in Knight's paper)
# input
#   x, y  : vectors of data
#   Pflag : 0 = silent,
#           1 = print no. ties, exchange count and other summeries
#           2 = print sorted data
#           3 = print sorting iterations
# output
#   x, y  : sorted data
#   N     : sample size (no. of pairs)
#   tau   : Kendall's tau
#   S     : numerator in Kendall's tau,
#             0.5*n*(n-1) - 2*(exchange count) when no ties
#   D     : denominator in Kendall's tau
#   T     : no. pairs of tied X
#   U     : no. pairs of tied Y
#   V     : no. pairs of tied (X, Y)
#   Pflag : print flag

# Knight's method, JASA 1966
#taucor= function(x,y,Pflag=0)
#{ if(is.matrix(x))
#  { d=ncol(x); n=nrow(x)
#    taumat=matrix(1,d,d)
#    Pflag=0
#    for(j1 in 1:(d-1))
#    { for(j2 in (j1+1):d)
#      { out= .C("ktau", as.double(x[,j1]), as.double(x[,j2]),
#            N=as.integer(n), tau=as.double(0), S=as.double(0),
#            D=as.double(0), T=as.integer(0), U=as.integer(0),
#            V=as.integer(0), Pflag=as.integer(Pflag))
#        taumat[j1,j2]=out$tau; taumat[j2,j1]=out$tau
#      }
#    }
#    return(taumat)
#  }
#  else 
#  { # x,y are vectors (no checks)
#    m=length(x)
#    n=length(y)
#    if(m==0 || n==0) stop("both 'x' and 'y' must be non-empty")
#    if(m!= n) stop("'x' and 'y' must have the same length")
#    if(n>50) Pflag = min(Pflag, 1)
#    out= .C("ktau", x=as.double(x), y=as.double(y),
#            N=as.integer(n), tau=as.double(0), S=as.double(0),
#            D=as.double(0), T=as.integer(0), U=as.integer(0),
#            V=as.integer(0), Pflag=as.integer(Pflag))
#    return(out$tau)
#  }
#}
