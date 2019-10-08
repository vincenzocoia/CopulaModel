# vine functions include the following
# varray2M : maximum array computed from vine array A, 
#         used for log-likelihood and simulation
# varray2NO : vine array to natural order
# vbin2array : generate vine array A from given b matrix
# varraycheck :
#   check whether A in natural order is a valid vine array if it has
#   permutation of 1:d on diagonal and permutation of 1;j in A[1:j,j]
#   Natural order also means A[j-1,j]=j-1 and A[j,j]=j

#============================================================

# A = dxd vine array of R-vine in standard order with 1:d on diagonal
# iprint = print flag for intermediate calculations
# str = string to describe vine for printing if iprint=T
# Output: list with two components  
#  M = dxd array with m_{kj}= max a_{k1},..,a_{kj}
#    actually M could be put in lower triangle of A
#  icomp = dxd indicator array on whether back step [k,j] is needed
#    icomp[k-1,m_{kj}=1 if  a_{kj}<m_{kj} for k>=2
varray2M=function(A,iprint=F,str="")
{ d=ncol(A)
  d1=d-1
  M=A
  icomp=matrix(0,d,d)
  for(k in 2:d1)
  { for(j in (k+1):d) M[k,j]=max(M[k-1,j],A[k,j]) }
  if(iprint) { cat("\n",str,"\n"); print(A); print(M) }
  for(k in 2:d1)
  { for(j in (k+1):d) 
    { if(A[k,j]<M[k,j]) icomp[k-1,M[k,j]]=1
    }
  }
  if(iprint) print(icomp)
  list(mxarray=M,icomp=icomp)
}  

#C= matrix(c(1,0,0,0,0, 1,2,0,0,0, 1,2,3,0,0, 1,2,3,4,0, 1,2,3,4,5), 5,5)
#D= matrix(c(1,0,0,0,0, 1,2,0,0,0, 2,1,3,0,0, 3,2,1,4,0, 4,3,2,1,5), 5,5)
#B0=matrix(c(1,0,0,0,0, 1,2,0,0,0, 1,2,3,0,0, 1,2,3,4,0, 1,3,2,4,5), 5,5)
#B1=matrix(c(1,0,0,0,0, 1,2,0,0,0, 1,2,3,0,0, 1,2,3,4,0, 2,1,3,4,5), 5,5)
#B2=matrix(c(1,0,0,0,0, 1,2,0,0,0, 2,1,3,0,0, 1,2,3,4,0, 1,2,3,4,5), 5,5)
#B3=matrix(c(1,0,0,0,0, 1,2,0,0,0, 1,2,3,0,0, 1,3,2,4,0, 2,1,3,4,5), 5,5)
#C6= matrix(c(1,0,0,0,0,0, 1,2,0,0,0,0, 1,2,3,0,0,0, 1,2,3,4,0,0, 1,2,3,4,5,0,
#   1,2,3,4,5,6),6,6)
#D6= matrix(c(1,0,0,0,0,0, 1,2,0,0,0,0, 2,1,3,0,0,0, 3,2,1,4,0,0, 4,3,2,1,5,0,
#   5,4,3,2,1,6),6,6)

#d=5
#A=matrix(0,d,d)
#diag(A)=1:5
#A[1,2:5]=c(1,2,2,4)
#A[2,3:5]=c(1,1,2)
#A[3,4:5]=c(3,1)
#A[4,5]=3
#varray2M(A1,iprint=T,"A")
#varray2M(C ,iprint=T,"C ")
#varray2M(D ,iprint=T,"D ")
#varray2M(B0,iprint=T,"B0")
#varray2M(B1,iprint=T,"B1")
#varray2M(B2,iprint=T,"B2")
#varray2M(B3,iprint=T,"B3")
#varray2M(C6,iprint=T,"C6")
#varray2M(D6,iprint=T,"D6")

#============================================================

# A = dxd vine array; 1:d on diagonal is NOT necessary
# irev = F means A1[d,d]=A[d,d]
# irev = T means A1[d,d]=A[d-1,d] (this option used to check if A is in
#                   equivalence class of size 1 or 2).
# iprint = print flag for intermediate calculations
# Output: vine array that has been converted to natural order;
# natural order means 1:d on diagonal and A[j-1,j]=j-1 for j=2,...,d
# Components of output are
#  $NOa = vine array with A[i,i]=A[i,i+1] without necessarily a sorted diagonal
#  $NO = vine array in natural order with sorted diagonal 1:d
#  $perm = permutation to get to from the NOa to NO
#  $diag = diagonal of NOa
varray2NO=function(A,irev=F,iprint=F)
{ d=nrow(A); d2=d-2; d1=d-1
  A1=matrix(0,d,d)
  T=vpartner(A)
  if(irev) { A1[d,d]=A[d1,d] } else { A1[d,d]=A[d,d] }
  for(k in d:2)
  { x=A1[k,k]
    for(ell in 1:(k-1)) A1[ell,k]=which(T[x,]==ell)
    T[x,]=0; T[,x]=0
    A1[k-1,k-1]=A1[k-1,k]
  }
  # A1 satisfies A[i,i]=A[i,i+1]
  if(iprint) print(A1)
  # now apply permutation
  iorder=order(diag(A1))
  A2=A1
  for(i in 1:d)
  { for(j in i:d) A2[i,j]=iorder[A1[i,j]] }
  if(iprint) print(A2)
  list(NOa=A1,NO=A2,perm=iorder,diag=diag(A1))
}

# A = dxd vine array; 1:d on diagonal is NOT necessary
# This function is not in NAMESPACE (not for direct use).
# Function with T(x,y)=ell if x-y are conditioned partners in tree l
vpartner= function(A)
{ d=nrow(A); 
  T=matrix(0,d,d) 
  for(j in 2:d)
  { x=A[j,j]
    for(k in 1:(j-1))
    { y=A[k,j]; T[x,y]=k; T[y,x]=k }
  }
  T
}

# This function is used for columns 5 or higher in vbin2array
# This function is not in NAMESPACE (not for direct use) 
# b0 = d-vector with length i-3, 
# A =  vine array
# i = column from 4 to ncol(A)
#   A has dimension at least ixi
# Output: ith column of A based on the binary representation b0 for column i
vstepb= function(b0,A,i)
{ itaken=rep(0,i)
  itaken[i]=1
  itaken[i-1]=1
  b=c(1,b0,1,1)
  #ac=i-1  # active column
  ac=i-2  # active column
  A1=A
  A1[i,i]=i
  A1[i-1,i]=i-1
  #for(k in (i-1):1) # older version with ac=i-1
  for(k in (i-2):1)
  { if(b[k]==1)
    { tem=A1[ac,ac]; itaken[tem]=1
      A1[k,i]=tem
      if(k>1)
      { ac=max((1:i)[itaken==0]) }
    }
    else
    { tem=A1[k-1,ac]; A1[k,i]=tem; itaken[tem]=1 }
  }
  #print(A1)
  A1[,i]
}

# Binary representation to vine array
# d = dimension
# b = dxd  matrix
# if(b==0) on input, then in
#    columns 4 to d, binary elements of b are randomly generated 
# iprint = print flag for intermediate calculations
# Calls are made to vstepb to generate a vine array from b[,]
# Output: dxd vine array A
vbin2array=function(d,b=0,iprint=F)
{ irandom=F
  if(class(b)=="numeric" | !is.matrix(b)) irandom=TRUE
  if(is.matrix(b)) { if (nrow(b)!=d | ncol(b)!=d)  irandom=TRUE }
  if(irandom)
  { b=matrix(NA,d,d) 
    b[1,]=1
    diag(b)=1
    for(i in 3:d) b[i-1,i]=1
    for(i in 4:d) b[2:(i-2),i]=rbinom(i-3,1,0.5) 
    if(iprint) print(b)
  }
  A=matrix(0,d,d)
  diag(A)=1:d; A[1,3]=1
  for(i in 2:d) A[i-1,i]=i-1 
  # column 4
  if(b[2,4]==1) # C-vine for first 4 columns 
  { A[1,4]=1; A[2,4]=2 }
  else { A[1,4]=2; A[2,4]=1 }  # D-vine for first 4 columns
  # columns 5 and higher
  for(i in 5:d) 
  { b0=b[2:(i-2),i]  # length i-3
    A[,i]=vstepb(b0,A,i)
  }
  A
}


# Inverse of vstepb function for column i.
# This function is not in NAMESPACE (not for direct use);
#    it is used by varraycheck
# A = vine array, dimension at least ixi
# i = column from 4 to ncol(A) 
# ichk0 > 0 means extra checks for column i
# Output:
#  if A is valid in column i, b is a binary vector with length i; 
#  if A is not valid in column i, -1 is returned.
vinvstepb= function(A,i,ichk0=0)
{ # do these basic checks first
  if(ichk0>0)
  { diagA=diag(A[1:i,1:i])
    if(max(diagA-(1:i))!=0) return(-1)
    for(k in 2:i) { if(A[k-1,k]!=k-1) return(-1) }
    if(A[1][3]!=1) return(-1);
  }

  b=rep(1,i)
  itaken=rep(0,i)
  itaken[i]=1
  itaken[i-1]=1
  ac=i-2  # active column
  for(k in (i-2):1)
  { if(A[k,i]==A[ac,ac])
    { b[k]=1; 
      tem=A[ac,ac]; itaken[tem]=1
      if(k>1) { ac=max((1:i)[itaken==0]) }
    }
    else if(A[k,i]==A[k-1,ac])
    { b[k]=0;
      tem=A[k-1,ac]; itaken[tem]=1 
    }
    else return(-1);   # not valid A in NO(i)
  }
  b
}

# split into varray2bin and varraycheck

# Check whether A in natural order is a valid vine array if it has
#  permutation of 1:d on diagonal and permutation of 1;j in A[1:j,j]
# Natural order also means A[j-1,j]=j-1 and A[j,j]=j
# A = dxd array with 1:d on diagonal
#   Calls to vinvstepb are made for columns 4 to ncol(d).
# Output is binary dxd array b with NA in lower triangle if A is valid
#   otherwise -1 is returned if A is invalid. 
varray2bin=function(A)
{ d=nrow(A)  
  b=matrix(NA,d,d)
  b[1,]=1
  diag(b)=1
  for(i in 3:d) b[i-1,i]=1
  for(i in 4:d) 
  { b0=vinvstepb(A,i)
    #print(b0)
    if(min(b0)==-1) return(-1)
    b[1:i,i]=b0
  }
  b
}

# generate a random vine array in dimension d
# d = dimension
genVineArray=function(d) { vbin2array(d,b=0) }

# Create dxd vine array based on a single number representation of the 
#   binary matrix representation
# The number bnum is converted to the binary array    
# d = integer>=3; if(d==3) bnum=0,
# bnum = integer between 0 and 2^((d-2)*(d-3)/2)-1 ;
#   2^((d-2)*(d-3)/2) is the number of d-dimensional vine arrays
#     in natural order.
#  bnum=0 is the D-vine, bnum=2^((d-2)*(d-3)/2)-1 is the C-vine
# iprint = print flag for intermediate steps
# To get bnum from a binary matrix representation bmat,
#   bvec as a binary vector starts increments from the last flexible position,
#   which is bmat[d-2,d];
# the (d-2)*(d-3)/2 positions in bmat: [2,4], [3,4],[3,5],...,[2,d],...,[d-2,d]
# Output: dxd vine array in natural order form
vnum2array=function(d,bnum=0,iprint=F)
{ dcase=(d-2)*(d-3)/2
  dpow=2^dcase
  A=matrix(0,d,d); diag(A)=1:d; A[1,3]=1
  for(i in 2:d) A[i-1,i]=i-1 
  if(d==3) return(A)

  bnum=floor(bnum)
  if(bnum<0 | bnum>=dpow) bnum=0
  bvec=d2b(dcase,bnum)
  b=matrix(NA,d,d) 
  b[1,]=1
  diag(b)=1
  for(i in 3:d) b[i-1,i]=1
  #for(i in 4:d) b[2:(i-2),i]=rbinom(i-3,1,0.5) 
  ii=0
  for(i in 4:d) { b[2:(i-2),i]=bvec[(ii+1):(ii+i-3)]; ii=ii+i-3 }
  if(iprint) print(b)
  
  # column 4
  if(b[2,4]==1) # C-vine for first 4 columns 
  { A[1,4]=1; A[2,4]=2 }
  else { A[1,4]=2; A[2,4]=1 }  # D-vine for first 4 columns
  # columns 5 and higher
  if(d>=5)
  { for(i in 5:d) 
    { b0=b[2:(i-2),i]  # length i-3
      A[,i]=vstepb(b0,A,i)
    }
  }
  A
}


# checks
#d=6
#d=floor(runif(1,5,10))
#A=vbin2array(d,0,iprint=T)
#print(A)
#b=varray2bin(A)
#print(b)

# various checks for validity of a vine array A with dimension d>=4
# A = dxd array
# Output:
# return 1 for OK
# return -3 for diagonal not 1:d
# return -2 for not permutation of 1:j in column j
# return -1 if cannot find proper binary array from array in natural order
varraycheck=function(A)
{ d=nrow(A)
  if(d!=ncol(A)) return(-1)
  if(sum(abs(sort(diag(A))-(1:d)))!=0) return(-3)
  # convert to 1:d on diagonal
  iorder=order(diag(A))
  A2=A
  for(i in 1:d)
  { for(j in i:d) A2[i,j]=iorder[A[i,j]] }
  #print(A2)
  for(j in 2:d)
  { if(sum(abs(sort(A2[1:(j-1),j])-(1:(j-1))))!=0) return(-2) }
  # next convert to natural order for more checks
  if(d<=3) return(1)
  #ANOobj=varray2NO(A2)
  #print(ANOobj)
  ANOobj=try(varray2NO(A2),silent=T)
  if(class(ANOobj)=="try-error") return(-1)
  b=varray2bin(ANOobj$NO) # if OK, a binary matrix is returned here
  if(is.matrix(b)) return(1) else return(-1)
}

# proper vine array
#A1=matrix(c(1,0,0,0,0,0, 1,2,0,0,0,0, 1,2,3,0,0,0, 1,2,3,4,0,0, 2,1,3,4,5,0,
#   2,1,4,3,5,6),6,6)
#b1=varraycheck(A1)
#print(b1)
## improper vine array
#A2=matrix(c(1,0,0,0,0,0, 1,2,0,0,0,0, 1,2,3,0,0,0, 1,2,3,4,0,0, 2,1,3,4,5,0,
#   2,3,1,4,5,6),6,6)
#b2=varraycheck(A2)
#print(b2)

#============================================================

# Additional functions for vines and vine arrays

# A = dxd vine array
# perm = permutation of 1:d 
# Output: vine array after permutation of indices
varrayperm=function(A,perm) 
{ A2=A
  d=ncol(A)
  for(i in 1:d)
  { for(j in i:d) A2[i,j]=perm[A[i,j]] }
  A2
}

# d = dimension
# iNO = F for standard order
# iNO = T for natural order
# Output: vine array for D-vine in standard order 1-2-..-d 
#   or in natural order
Dvinearray=function(d,iNO=F)
{ D=matrix(0,d,d)
  diag(D)=1:d
  if(iNO) { D=vnum2array(d,bnum=0) }
  else { for(j in 2:d) D[1:(j-1),j]=(j-1):1 }
  D
}

# d = dimension
# Output: vine array for C-vine in standard order 1-2-..-d
Cvinearray=function(d)
{ C=matrix(0,d,d)
  for(j in 1:d) C[j,j:d]=j
  C
}

