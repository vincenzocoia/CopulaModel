! conversion of some functions in partialcorr.R  to f90
! gfortran -o pcor pcor.f90

!program mainex
!  implicit none
!  integer d,ii,no1,first1,prev,ncond,i,j,ell1,indx
!  integer, dimension(:), allocatable :: zeros,svec
!  integer, dimension(:,:), allocatable :: A1
!  double precision, dimension(:), allocatable :: logdet
!  double precision, dimension(:,:), allocatable :: r
!  double precision, dimension(:,:,:), allocatable :: pc3

!  d=5;
!  ncond=2**d-4  
!  allocate ( zeros(d),svec(d) )
!  allocate ( pc3(d,d,ncond), r(d,d), A1(d,d), logdet(d-1) )
!  do i=1,d
!    do j=1,d
!       r(i,j)=0.5**(abs(i-j))
!    end do
!  end do
  !do ii=1,5
  !  call d2brev(d, ii, zeros, no1, first1, prev )
  !end do
!  call allpcor(d,ncond,r,pc3)
!  do i=1,d
!    print '(5f10.6)', pc3(i,:,1)
!  end do
!  A1=0
!  A1(1,:)=(/1,1,1,2,3/)
!  A1(2,2:5)=(/2,2,1,1/)
!  A1(3,3:5)=(/3,3,2/)
!  A1(4,4:5)=(/4,4/)
!  A1(5,5)=5
!  do i=1,d
!    print *, A1(i,:)
!  end do
  !do ii=1,10
  !  read *, ell1, (svec(j),j=1,ell1)
  !  call subset2index(d,ell1,svec,indx)
  !  print *, indx
  !end do
!  call trvinelogdet(d,A1,r,pc3,logdet)
!  print *, logdet
!  deallocate ( zeros, pc3, r, A1, svec, logdet )
!  stop
!  end

! decimal to binary, reverse order
! vector of positions of 0's, output number of 1's, string of positions of 1,
! first position of 1, decimal code for prev condition with first 1 set to 0
! inputs
!   n = dimension of vector, 
!   ii = integer between 0 and 2^n-1
! outputs (binary vector representation of ii via components)
!   zeros = positions of zeros in binary representation
!   no1 = number of 1s in the binary representation
!   first1 = position of first 1
!   prev = representation for previous decimal where position first1 is 0 
subroutine d2brev(n, ii, zeros, no1, first1, prev)
  implicit none
  integer n,ii,zeros(n), no1, first1, prev
  integer d,min1,k,i
  integer, dimension(:), allocatable :: jj
  allocate ( jj(n) )
  ! initialize to avoid warning
  min1=0
  d=ii
  jj=0
  zeros=0
  k=0
  no1=0
  do i=1,n
    jj(i)=mod(d,2); d=int(d/2); 
    if(jj(i)==0) then
      k=k+1; zeros(k)=i 
    else 
      no1=no1+1 
    endif
  end do
  A:do i=1,n 
    if(jj(i)==1) then
      min1=i; exit A;  ! breaking from a loop in f90
    endif  
  end do A
 
  ! decimal if first 1 changed to 0
  !zeros=zeros[1:(n-no1)]
  ! decimal code for prev
  prev=ii-2**(min1-1)
  !print *
  !print *,"index=", ii
  !print *,jj
  !print *,"number of 1s is", no1
  !print *,"position of first 1 is", min1
  !print *,zeros
  !print *,"code for prev using in cond is ", prev
  deallocate ( jj ) ! binary vector not returned
  first1=min1
  !list(zeros=z,no1=no1,str=str,first1=min1,prev=prev)
  return 
  end

! input
!   d = dimension
!   ncond = index for the conditioning set
!   r = dxd correlation matrix
! outputs (partial correlations given conditioning set indexed as ncond)
!   aa = aa(:,:,ncond) is completed in the 3-dimensional array
subroutine allpcor(d,ncond,r,aa)
  implicit none
  integer d,ncond
  double precision r(d,d),aa(d,d,ncond)
  integer k,j,kk,i,ll,no1,first1,prev,pp,oo,zi,zj
  integer, dimension(:), allocatable :: no1vec,zeros
  double precision tem
  !ncond=2**d-4  
  allocate ( no1vec(ncond), zeros(d) )
  no1vec=0
  aa=0
  ! conditioning on 1 variable
  do k=1,d
    kk=2**(k-1)
    aa(d,d,kk)=1
    do i=1,d-1
      aa(i,i,kk)=1
      do j= (i+1),d
        if(.not.(i==k .or. j==k)) then
          tem=(r(i,j)-r(i,k)*r(j,k))
          aa(i,j,kk)=tem/sqrt((1-r(i,k)**2)*(1-r(j,k)**2))
          aa(j,i,kk)=aa(i,j,kk)
        end if
      end do
    end do
    aa(k,k,kk)=0
    !print *
    !do i=1,d
    !  print '(5f10.6)', aa(i,:,kk)
    !end do
  end do
  
  ! conditioning on 2 (or more) variables
  !print *, "conditioning on 2 (or more) variables"
  do k = 1,ncond
    call d2brev(d, k, zeros, no1, first1, prev )
    no1vec(k)=no1
    !if(kinfo$no1==1) next
    !if(kinfo$no1==(d-1)) next
    if(no1>1 .and. no1<d-1) then
      kk=k
      pp=prev
      oo=first1
      !ll=length(z) replace by d-kinfo$no1
      ll=d-no1
      do i = 1,(ll-1)
        zi=zeros(i)
        do j = (i+1),ll
          zj=zeros(j)
          tem= aa(zi,zj,pp)-aa(zi,oo,pp)*aa(zj,oo,pp)
          aa(zi,zj,kk)=tem/sqrt((1-aa(zi,oo,pp)**2)*(1-aa(zj,oo,pp)**2))
          aa(zj,zi,kk)=aa(zi,zj,kk)
        end do
      end do
      do i = 1,ll
        aa(zeros(i),zeros(i),kk)=1
      end do
      !print *
      !do i=1,d
      !  print '(5f10.6)', aa(i,:,kk)
      !end do
    end if
  end do
  deallocate ( no1vec, zeros )
  return
  end

! mapping of given index to 3rd dimension of output of allpcor
! inverse of d2brev
! input
!   d = dimension
!   ell1 = length of svec
!   svec = vector that is subset of 1,...,d of length ell1
! output
!   indx = index of the subset svec in lexicographic ordering 
subroutine subset2index(d,ell1,svec,indx)
  implicit none
  integer d,ell1,svec(ell1)
  integer j,pw,indx
  integer, dimension(:), allocatable :: ii
  allocate ( ii(d) )
  ii=0
  do j=1,ell1
   ii(svec(j))=1
  end do
  !  pw=2^(0:(d-1))
  !  sum(ii*pw)
  pw=1; indx=0
  do j=1,d
    indx=indx+pw*ii(j)
    pw=pw*2
  end do
  deallocate ( ii )
  return
  end

!============================================================

! vinearray to logdeterminants for 1-truncation, 2-truncation,... 
! input
!   d = dimension
!   A = dxd vine array with 1;d on diagonal
!   rr = dxd correlation matrix
!   pc3 = d x d x (2^d-4) with partial correlation matrix in pc3[,,kk]
! outputs
!   logdet = vector of log determinants of truncated vines, length d-1
!   pcm = dxd array with partial correlations, position (ell,j)
!          for the jth variable and tree ell (1<=ell<j)
subroutine trvinelogdet(d,A,rr,pc3,logdet,pcm)
  implicit none
  integer d,A(d,d)
  double precision rr(d,d),pc3(d,d,2**d-4),logdet(d-1),pcm(d,d)
  integer ncond,j,ell,ell1,igiven
  !double precision, dimension(:,:), allocatable :: pcm
  double precision margsumlog,cum
  !double precision, dimension(:,:,:), allocatable :: pc3

  ncond=2**d-4  
  !allocate ( pc3(d,d,ncond), pcm(d,d) )
  !allocate ( pcm(d,d) )
  !do i=1,d
  !  print '(5f10.6)', rr(i,:)
  !end do
  !call allpcor(d,ncond,rr,pc3)
  !do i=1,d
  !  print '(5f10.6)', pc3(i,:,1)
  !end do
  pcm=0
  ! correlations in row 1
  do j = 2,d
    pcm(1,j)=rr(A(1,j),A(j,j)) 
  end do
  ! pcor in row ell given A[1:(ell-1),j]
  do ell = 2,(d-1)
    ell1=ell-1
    do j=(ell+1),d
      call subset2index(d,ell1,A(1:ell1,j),igiven)
      !print *, ell,j,igiven
      pcm(ell,j)=pc3(A(ell,j),A(j,j),igiven)
    end do
  end do
  !do ell=1,d
  !  print *, pcm(ell,:)
  !end do
  
  cum=0.d0
  do ell=1,d-1
    margsumlog=sum(log(1.d0-pcm(ell,(ell+1):d)**2))
    cum=cum+margsumlog
    logdet(ell)=cum
  end do
  ! margsumlog=apply(log(1-pcm^2),1,sum)
  !  logdet=cumsum(margsumlog)
  !  logdet[1:(d-1)]
  !deallocate ( pcm )
  return
  end


