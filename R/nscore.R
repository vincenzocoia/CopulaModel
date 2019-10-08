
# Transform each variable to normal scores 
# data = dataframe or matrix, or vector
# iopt = T to use adjustment 'a' in nscoreOpta()
# Output: matrix or vector of normal scores
nscore=function(data,iopt=F)
{ if(is.vector(data))
  { nr=length(data)
    if(iopt) anormal=nscoreOpta(nr) else anormal=-0.5
    qn=qnorm(((1:nr)+anormal)/(nr+1+2*anormal))
    jj=rank(data)
    out=qn[jj]
  }
  else
  { nc=ncol(data)
    nr=nrow(data)
    if(iopt) anormal=nscoreOpta(nr) else anormal=-0.5
    out=matrix(0,nr,nc)
    qn=qnorm(((1:nr)+anormal)/(nr+1+2*anormal))
    for(j in 1:nc)
    { jj=rank(data[,j])
      tem=qn[jj]
      out[,j]=tem
    }
  }
  out
}

# normal score transform to get variance of 1 and mean of 0
# find optimal value of 'a' between -.75 and -.5
# n = sample size
# mxiter = maximum number of Newton-Raphson iteration
# eps = tolerance for convergence
# iprint = print flag for intermediate results
# Output: the adjustment 'a' for qnorm( (rank+a)/(n+1+2*a) ) scores
nscoreOpta=function(n, mxiter=20,eps=1.e-4,iprint=F)
{ a0=-0.55
  iter=0
  n1=floor((n+3)/2)
  if(iprint) print(c(n1,n))
  ii=(n1:n)
  diff=1
  while(iter<mxiter & abs(diff)>eps) 
  { den=n+1+2*a0
    z=qnorm((ii+a0)/den)
    g=2*sum(z^2)/n-1
    tem=(n+1-2*ii)/den^2
    gp=2*2*sum(tem*z/dnorm(z))/n
    diff=g/gp; a0=a0-diff
    while(a0< -.8 | a0> -.4)
    { diff=diff/2; a0=a0+diff }
    iter=iter+1
    if(iprint) print(c(iter,a0,g,gp))
  }
  if(iter>=mxiter)  cat("did not converge\n")
  a0
}

# Transform each variable to uniform scores 
# data = dataframe or matrix, or vector
# aunif = adjustment 'a' for (rank+a)/(n+1+2*a) as scores. n=sample size 
# Output: matrix or vector of uniform scores
uscore=function(data,aunif=-0.5)
{ if(is.vector(data))
  { nr=length(data)
    us=((1:nr)+aunif)/(nr+1+2*aunif)
    jj=rank(data)
    out=us[jj]
  }
  else
  { nc=ncol(data)
    nr=nrow(data)
    out=matrix(0,nr,nc)
    us=((1:nr)+aunif)/(nr+1+2*aunif)
    for(j in 1:nc)
    { jj=rank(data[,j])
      tem=us[jj]
      out[,j]=tem
    }
  }
  out
}

