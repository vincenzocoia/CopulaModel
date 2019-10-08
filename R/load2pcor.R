# Functions for loading matrix of Gaussian factor models

# Convert from loading matrix amat to partial correlations for second 
#   and subsequent factors
# amat = dxp loading matrix with p>=2
# Output: matrix where columns 2,...,p are partial correlations given previous
#  factors 
load2pcor=function(amat)
{ rhmat=amat
  rhmat[,2]=amat[,2]/sqrt(1-amat[,1]^2)
  p=ncol(amat)
  if(p>2)
  { psiv=1-amat[,1]^2-amat[,2]^2
    for(k in 3:p)
    { rhmat[,k]=amat[,k]/sqrt(psiv)
      psiv=psiv-amat[,k]^2
    }
  }
  rhmat
}

# Partial correlation representation to loadings for p-factor
# rhmat = dxp matrix for pcor2load, correlations with factor 1 in column 1, 
#   partial correlations with factor k given previous factors in column k
# Output: loading matrix
pcor2load=function(rhmat)
{ if(is.vector(rhmat) | ncol(rhmat)==1) return(rhmat)
  # p>=2
  d=nrow(rhmat)
  p=ncol(rhmat)
  amat=matrix(0,d,p)
  amat[,1]=rhmat[,1]
  tem=rep(1,d)
  for(kf in 2:p)
  { temc=sqrt(1-rhmat[,kf-1]^2)
    amat[,kf]=rhmat[,kf]*tem*temc 
    tem=tem*temc
  }
  amat
}

#============================================================

# Apply a sequence of Givens rotations to a loadings matrix to get
# upper right triangle of zeros

# Givens rotations for loading matrix with 2 columns
#  amat = dx2 loading matrix
#  row = index of row to set to 0 in second column
#  iprint = print flag for intermediate steps
# Output: rotated loading matrix
grotate2=function(amat,row,iprint=F)
{ d=nrow(amat)
  a=amat
  a1=amat
  i=row
  # make a[i,2]=0
  den=sqrt(a[i,1]^2+a[i,2]^2)
  g11=a[i,1]/den; g12=a[i,2]/den
  a1[,1]= a[,1]*g11+a[,2]*g12
  a1[,2]= -a[,1]*g12+a[,2]*g11
  a=a1
  # change sign of columns
  if(sum(a[,1]<0)>d/2) a[,1]=-a[,1] 
  if(sum(a[,2]<0)>d/2) a[,2]=-a[,2] 
  #if(iprint) print(a)
  a
}

# Givens rotations for loading matrix with 3 columns
# Givens rotations operator (as right-side operators) on 2 columns at at time
#  amat = dx3 loading matrix
#  row1 = index of row to set to 0 in second and third columns
#  row2 = index of row to set to 0 in third column
#  iprint = print flag for intermediate steps
# Output: rotated loading matrix
grotate3=function(amat,row1=1,row2=2,iprint=F)
{ d=nrow(amat)
  a=amat
  a1=amat
  i=row1
  # make a[i,2]=0
  den=sqrt(a[i,1]^2+a[i,2]^2)
  g11=a[i,1]/den; g12=a[i,2]/den
  a1[,1]= a[,1]*g11+a[,2]*g12
  a1[,2]= -a[,1]*g12+a[,2]*g11
  a=a1
  if(iprint) print(a)
  # make a[i,3]=0
  den=sqrt(a[i,1]^2+a[i,3]^2)
  g11=a[i,1]/den; g13=a[i,3]/den
  a1[,1]= a[,1]*g11+a[,3]*g13
  a1[,3]= -a[,1]*g13+a[,3]*g11
  a=a1
  if(iprint) print(a)
  i2=row2
  # make a[i2,3]=0
  den=sqrt(a[i2,2]^2+a[i2,3]^2)
  g22=a[i2,2]/den; g23=a[i2,3]/den
  a1[,2]= a[,2]*g22+a[,3]*g23
  a1[,3]= -a[,2]*g23+a[,3]*g22
  a=a1
  if(iprint) print(a)
  a
}

