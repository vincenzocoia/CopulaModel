
! outer product of two vectors
! input
!   na = length of avec
!   nb = length of bvec
!   avec = vector 1
!   bvec = vector 2
! output
!   abmat = na x nb matrix with outer product of avec and bvec 
subroutine outer(na,nb,avec,bvec,abmat)
  implicit none
  integer na,nb,i
  double precision avec(na),bvec(nb),abmat(na,nb)
  ! na=size(avec); nb=size(bvec)
  do i =1,na  
    abmat(i,:)=avec(i)*bvec
  end do
  return
  end

! decimal to vector
! inputs
!   d = length of vector
!   ncat = number of categories
!   ii = integer between 1 and ncat^d 
! output:
!   jj= integer vector with values in 0,1,...,ncat-1  corresponding to ii 
! for d=2, ncat=3, ii=1,..,9 leads to (0,0), (0,1), ..., (2,2)
! for d=3, ncat=2, ii=1,..,8 leads to (0,0,0), (0,0,1), ..., (1,1,1)
  subroutine d2v(d,ncat,ii,jj)
  implicit none
  integer d,ncat,ii,jj(d)
  integer i,tem
  tem=ii
  do i=d,1,-1
    jj(i)=mod(tem,ncat)
    tem=tem/ncat
  end do
! tem%ncat=mod(tem,ncat) is the remainder from dividing ncat into tem
!    same as tem%%ncat in R
! tem/ncat is integer division rounded down
!    same as floor(tem/ncat) in R
  return
  end

! vector to decimal
! inputs
!   d = length of vector
!   ncat = number of categories
!   jj= integer vector with values in 0,1,...,ncat-1 
! output:
!   integer between 1 and ncat^d, lexicographic order of jj
  integer function v2d(d,ncat,jj)
  implicit none
  integer d,ncat,jj(d)
  integer i,pw,s
  pw=1
  s=0
  do i=d,1,-1
    s=s+pw*jj(i)
    pw=pw*ncat
  end do
  v2d=s
  return
  end

