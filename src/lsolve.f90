! linear system solver ; Gaussian elimination with partial pivoting
! gfortran -o lsolve lsolve.f90
!program maingepp
!  implicit none
!  integer i,j,n,nrhs
!  double precision, dimension(:,:), allocatable :: a,rhs,soln
!  double precision det,tol
!  tol=1.d-6
!  read *, n,nrhs
!  allocate ( a(n,n),rhs(n,nrhs),soln(n,nrhs) )
!  do i=1,n
!    do j=1,n
!      a(i,j)=0.5**iabs(i-j) 
!    end do
!    do j=1,nrhs
!      rhs(i,j)=i*j
!    end do
!  end do
!  do i=1,n
!    print '(8f10.5)', a(i,:), rhs(i,:)
!  end do
!  call lsolve(n,nrhs,a,rhs,tol,soln,det)
!  print *, "Solution A^{-1} rhs = soln"
!  do i=1,n
!    print '(8f10.5)', soln(i,:)
!  end do
!  print *,det
!  stop
!  end

!a=matrix(c(1,.5,.25,.5,1,.5,.25,.5,1),3,3)
!solve(a,1:3)
! det(a)

! lsolve with gaussian elimination with partial pivoting
! inputs  
!   n = dimension of square matrix a
!   nrhs = #right-hand side vectors
!   a = n x n matrix
!   rhs = n x nrhs matrix
!   tol = tolerance for singularity, set to something like 1.d-5 or smaller
! outputs 
!   soln = n x nrhs matrix with solution to a*soln=rhs
!   det = determinant(a) 
subroutine lsolve(n,nrhs,a,rhs,tol,soln,det)
  implicit none
  integer n,nrhs
  double precision a(n,n),rhs(n,nrhs),soln(n,nrhs),det,tol
  integer parity,i,j,k,l,mrow,nk,n1,ll,k1
  double precision mx,tem,sm
  double precision , dimension(:,:), allocatable :: aa
  nk=n+nrhs
  allocate ( aa(n,nk) )
  n1=n+1
  ! initialize aa = (a | rhs )
  aa(:,1:n)=a
  aa(:,n1:nk)=rhs
  parity=1
  do k=1,n-1
    mx=abs(aa(k,k))
    mrow=k
    k1=k+1
    do l=k1,n
      if(abs(aa(l,k))>mx) then
        mx=abs(aa(l,k))
        mrow=l
      end if
    end do
    if(mx<=tol) then
      det=0.d0
      return
    end if
    if(mrow>k) then
      parity=-parity
      do i=1,nk
        tem=aa(mrow,i)
        aa(mrow,i)=aa(k,i)
        aa(k,i)=tem
      end do
    end if
    do i=k1,n
      aa(i,k)=aa(i,k)/aa(k,k)
      do j=k1,nk
        aa(i,j)=aa(i,j)-aa(i,k)*aa(k,j)
      end do
    end do
  end do
  if(abs(aa(n,n))<=tol) then
    det=0.d0
    return
  end if
  det=parity
  do i=1,n
    det=det*aa(i,i)
  end do
  do ll=n1,nk
    aa(n,ll)=aa(n,ll)/aa(n,n)
    do i=n-1,1,-1
      sm=0.d0
      do j=i+1,n
        sm=sm+aa(i,j)*aa(j,ll)
      end do
      aa(i,ll)=(aa(i,ll)-sm)/aa(i,i)
    end do
  end do
  soln=aa(:,n1:nk)
  deallocate (aa)
  return
  end
