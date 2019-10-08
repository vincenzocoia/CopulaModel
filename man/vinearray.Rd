\name{vinearray}
\Rdversion{1.1}
\alias{vinearray}
\alias{varray2NO}
\alias{vbin2array}
\alias{varray2bin}
\alias{genVineArray}
\alias{vnum2array}
\alias{varrayperm}
\alias{varraycheck}
\alias{Cvinearray}
\alias{Dvinearray}
\title{
Vine array conversion, validation and generation
}
\description{
Vine array conversion, validation and generation
}
\usage{
varray2NO(A,irev=F,iprint=F) # vine array to natural order
vbin2array(d,b=0,iprint=F)   # binary representation to vine aray
varray2bin(A)  # vine array to binary representation
varraycheck(A) # vine array validation check
genVineArray(d)  # generate random vine array
varrayperm(A,perm) # vine array with permutation of indices
vnum2array(d,bnum=0,iprint=F) # decimal (of binary representation) to vine array
Dvinearray(d,iNO=F) # D-vine array
Cvinearray(d)       # C-vine array
}
\arguments{
\item{A}{dxd vine array: only the upper triangle is used.
For varray2bin, A should be in natural order, 1:d on diagonal.
and also A[j-1,j]=j-1 (j=2,...,d).
For varray2NO, varraycheck and varrayperm, A is not required to have 1:d 
on diagonal.}
\item{irev}{irev=F means A1[d,d]=A[d,d], irev=T means A1[d,d]=A[d-1,d]}
\item{iprint}{print flag for some intermediate steps}
\item{d}{dimension of vine array}
\item{b}{binary dxd matrix, the important entries are in columns 4 to d;
rows 2 to j-2 in column j.}
\item{perm}{d-dimension vector, permutation of 1:d}
\item{bnum}{for vnum2array,
integer in the range 0 and 2^((d-2)*(d-3)/2)-1, used to create the 
binary matrix b; vnum2array allows easier enumeration through all vine
arrays in natural order for small d=5,6,7. }
\item{iNO}{Dvinearray in natural order if iNO=T, and in standard 1-2-...-d
order if iNO=F}
}
\value{
\item{Aobj}{varray2NO returns a list with $NOa which is a vine array
with A[i,i]=A[i,i+1] without necessarily a sorted diagonal,
$NO with vine array in natural order with sorted diagonal 1:d,
$perm is the permutation to get to from the NOa to NO,
$diag is the diagonal of NOa.
}
\item{ANO}{if b is an array in
vbin2array, then a vine array in NO=natural order corresponding
to b is returned; if b=0 in vbin2array (or equivalently
genVineArray, then a random vine array is returned.
}
\item{bNO}{if A is a vine array in NO, varray2bin the
binary matrix, otherwise varray2bin returns -1 to indicate that
the array is not in NO.
}
\item{code}{for varraycheck,
code= 1 for valid vine array;
code= -3 for diagonal not 1:d;
code= -2 for not permutation of 1:j in column j;
code= -1 if cannot find proper binary array from array in natural order.
}
\item{A}{for Dvinearray, Cvinearray,
genVineArray and vnum2array, the output is a dxd vine array.}
\item{Aperm}{for varrayperm, the output is a dxd vine array with permuted
indices, perm[] is on the diagonal.}
}
\seealso{
\code{\link{gausstrvine}}
\code{\link{rvinesim}}
}
\references{
Joe H, Cooke RM and Kurowicka D (2011). Regular vines: generation algorithm 
and number of equivalence classes.
In Dependence Modeling: Vine Copula Handbook, pp 219--231.
World Scientific, Singapore.
}
\details{In the VineCopula R package and in
"Dissmann J, Brechmann EC, Czado C and Kurowicka D (2013).
Computational Statistics and Data Analysis, 59, 52-69",
the vine array is abbreviation as RVM and it is A[d:1,d:1]
where d=ncol(A). This conversion is needed if this package and
VineCopula are used.
}
\examples{
C5= matrix(c(1,0,0,0,0, 1,2,0,0,0, 1,2,3,0,0, 1,2,3,4,0, 1,2,3,4,5), 5,5)
D5= matrix(c(1,0,0,0,0, 1,2,0,0,0, 2,1,3,0,0, 3,2,1,4,0, 4,3,2,1,5), 5,5)
D5perm= matrix(c(5,0,0,0,0, 5,4,0,0,0, 4,5,3,0,0, 3,4,5,2,0, 2,3,4,5,1), 5,5)
varrayperm(D5,perm=c(5,4,3,2,1))  # same as D5perm above
print(varray2NO(C5))
print(varray2NO(D5))
DNO=varray2NO(D5perm)
print(DNO)
set.seed(123)
vbin2array(d=5,b=0,iprint=FALSE)  # random vine array in NO=natural order
set.seed(123)
genVineArray(d=5)
b=matrix(0,5,5)
b[2:3,5]=1
vbin2array(d=5,b=b,iprint=F)  # vine array with binary matrix b
for(i in 0:7) { print(vnum2array(5,bnum=i,iprint=FALSE)) }
# more checks
for(i in 1:3)
{ d=floor(runif(1,5,10))
  A=vbin2array(d,0,iprint=TRUE)
  print(A)
  b=varray2bin(A)
  print(b)
}

# interface varraycheck to varray2NO, rename varray2bin
varraycheck(C5)
varraycheck(D5) 
varraycheck(DNO$NO)
# proper vine array
A1=matrix(c(1,0,0,0,0,0, 1,2,0,0,0,0, 1,2,3,0,0,0, 1,2,3,4,0,0, 2,1,3,4,5,0,
   2,1,4,3,5,6),6,6)
b1=varraycheck(A1)
print(b1)
A1perm=varrayperm(A1,c(1,4,5,2,6,3))
varraycheck(A1perm) # also OK
# improper vine array, code=-1
A2=matrix(c(1,0,0,0,0,0, 1,2,0,0,0,0, 1,2,3,0,0,0, 1,2,3,4,0,0, 2,1,3,4,5,0,
   2,3,1,4,5,6),6,6)
b2=varraycheck(A2)
print(b2)
# improper vine array, code=-2
A3=matrix(c(1,0,0,0,0,0, 1,2,0,0,0,0, 1,2,3,0,0,0, 1,2,3,4,0,0, 2,1,3,4,5,0,
   2,1,4,3,3,6),6,6)
b3=varraycheck(A3)

# Dvinearray
Dvinearray(5)
Dvinearray(5,iNO=TRUE)
vnum2array(5,bnum=0) # same as above 
# Cvinearray
Cvinearray(5)
vnum2array(5,bnum=7) # same as above 
}
\keyword{vine}
