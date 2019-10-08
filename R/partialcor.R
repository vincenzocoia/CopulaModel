# functions for partial correlations

# rr = dxd correlation matrix
# Outputs: all partial correlations in 2 forms: a 3-dim array and a matrix,
#    vector min/max partial correlations conditioning on set size
allpcor=function(rr)
{ r=rr
  d=nrow(r)
  ncond=2^d-4  
  # two ways of storing outputs
  aa=array(NA,c(d,d,ncond))
  bb=matrix(0,d*(d-1)/2,ncond)
  bb=as.data.frame(bb)
  no1vec=rep(0,ncond)

  # conditioning on 1 variable
  for(k in 1:d)
  { kk=2^(k-1)
    for(i in 1:(d-1))
    { for(j in (i+1):d)
      { if(i==k | j==k) next
        tem=(r[i,j]-r[i,k]*r[j,k])
        aa[i,j,kk]=tem/sqrt((1-r[i,k]^2)*(1-r[j,k]^2))
        aa[j,i,kk]=aa[i,j,kk]
      }
    }
    diag(aa[,,kk])=1
    aa[k,k,kk]=0
  }

  # conditioning on 2 (or more) variables
  for(k in 1:ncond)
  { kinfo=d2brev(d,k)
    colnames(bb)[k]=kinfo$str
    no1vec[k]=kinfo$no1
    if(kinfo$no1==1) next
    if(kinfo$no1==(d-1)) next
    kk=k
    pp=kinfo$prev
    oo=kinfo$first1
    z=kinfo$z
    ll=length(z)
    #cat("k=", k, " length(z)=", ll,"\n")
    for(i in 1:(ll-1))
    { zi=z[i]
      for(j in (i+1):ll)
      { zj=z[j]
        tem=(aa[zi,zj,pp]-aa[zi,oo,pp]*aa[zj,oo,pp])
        aa[zi,zj,kk]=tem/sqrt((1-aa[zi,oo,pp]^2)*(1-aa[zj,oo,pp]^2))
        aa[zj,zi,kk]=aa[zi,zj,kk]
      }
    }
    for(i in 1:ll) aa[z[i],z[i],kk]=1
  }

  # fill in elements of bb
  cc=0
  for(j in 2:d)
  { for(i in 1:(j-1))
    { cc=cc+1
      bb[cc,]=aa[i,j,]
      rownames(bb)[cc]=paste(i,j,sep=",")
    }
  }

  # get max/min partial corr for each conditioning set size
  mnmx=matrix(0,2,d-2)
  for(sz in 1:(d-2))
  { isz=(no1vec==sz)
    #print(isz)
    tem=c(as.matrix(bb[,isz]))
    #print(tem)
    tem=tem[!is.na(tem)]
    mnmx[1,sz]=min(tem)
    mnmx[2,sz]=max(tem)
  }

  jvec=(no1vec<(d-1))
  bb=bb[,jvec] # eliminate null columns with conditioning set of size d-1
  list(pc3.array=aa,pc2.array=bb,condsize=no1vec[jvec],mnmx=mnmx)
}

#============================================================

# S = covariance or correlation matrix
# given = vector indices for the given or conditioning variables
# j,k = indices for the conditioned variables
# Output: partial correlation of variables j,k given indices in 'given'
partcor=function(S,given,j,k)
{ S11=S[given,given]
  jk=c(j,k)
  S12=S[given,jk]
  S21=S[jk,given]
  S22=S[jk,jk]
  if(length(given)>1) { tem=solve(S11,S12); Om212=S21%*%tem }
  else { tem=S12/S11; Om212=outer(S21,tem) }
  #om11=1-Om212[1,1]; om22=1-Om212[2,2]  # correct only for correlation matrix
  om11=S[j,j]-Om212[1,1]
  om22=S[k,k]-Om212[2,2]
  om12=S[j,k]-Om212[1,2]
  om12/sqrt(om11*om22)
}

#============================================================


# vine array to log determinants for 1-truncation, 2-truncation,... 
# A = dxd vine array with 1:d on diagonal
# rr = dxd correlation matrix
# pc3 = dxdx(2^d-4) with partial correlation matrix in pc3[,,kk]
# Output: vector of log determinants
# This function is not in NAMESPACE, it has been replaced by f90 code.
trvine.logdet= function(A,rr,pc3)
{ d=nrow(rr)
  #pc3=allpcor(rr)
  #pc3=pc3$pc3.array
  pcm=matrix(0,d,d)
  # correlations in row 1
  for(j in 2:d) pcm[1,j]=rr[A[1,j],A[j,j]] 
  # pcor in row ell given A[1:(ell-1),j]
  for(ell in 2:(d-1))
  { ell1=ell-1
    for(j in (ell+1):d)
    { igiven=subset2index(d,A[1:ell1,j])
      #print(c(ell,j,igiven))
      pcm[ell,j]=pc3[A[ell,j],A[j,j],igiven]
    }
  }
  margsumlog=apply(log(1-pcm^2),1,sum)
  logdet=cumsum(margsumlog)
  logdet[1:(d-1)]
}

# inverse of d2brev
# d = dimension
# svec = vector that is subset of 1,...,d
# Output: index of the subset svec in lexicographic ordering
subset2index= function(d,svec)
{ ii=rep(0,d)
  ii[svec]=1
  pw=2^(0:(d-1))
  sum(ii*pw)
}

# This function is used by allpartcor()
# decimal to binary, reverse order (this function not in NAMESPACE)
# vector of positions of 0's, output number of 1's, string of positions of 1,
# first position of 1, decimal code for prev condition with first 1 set to 0
# n = positive integer >=3;
# ii = integer between 0 and 2^n-1 inclusive
#  there is no checking for out of range
# iprint = print flag for intermediate results
# Output: binary vector representation of ii via a list:
#   zeros= positions of zeros in binary rep
#   no1=number of 1s,
#   str=string
#   first1=position of first 1
#   prev=representation for previous decimal where position first1 is 0
d2brev= function(n,ii,iprint=F)
{ d=ii
  jj=rep(0,n)
  z=rep(0,n)
  str=""
  k=0
  no1=0
  for(i in 1:n)
  { jj[i]=d%%2; d=floor(d/2); 
    if(jj[i]==0) 
    { k=k+1; z[k]=i }
    else { str=paste(str,i,sep=","); no1=no1+1 }
  }
  for(i in 1:n) { if(jj[i]==1) { min1=i; break } }
  # decimal if first 1 changed to 0
  # string of 0 positions
  # kk=(1:n)[jj==0],  # not needed ??
  str=substr(str,2,2*no1)
  z=z[1:(n-no1)]
  # decimal code for prev
  prev=ii-2^(min1-1)
  if(iprint)
  { cat("\nindex=", ii,"\n")
    print(jj)
    cat("number of 1s is", no1,"\n")
    cat("position of first 1 is", min1,"\n")
    cat("string is ",str,"\n")
    print(z)
    cat("code for prev using in cond is ", prev,"\n")
  }
  list(zeros=z,no1=no1,str=str,first1=min1,prev=prev)
}


#============================================================
