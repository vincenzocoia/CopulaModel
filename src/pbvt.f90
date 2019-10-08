
! sample main program  (this is f90 calling f77)
!  gfortran -c cdfbvt.f
!  gfortran -o pbvt pbvt.f90 cdfbvt.o
!program bvtchk
!  implicit none
!  integer n,i;
!  double precision, dimension(:), allocatable :: z1,z2,r,bcdf
!  integer, dimension(:), allocatable :: nu
!  read *, n
!  allocate ( z1(n), z2(n), r(n), bcdf(n), nu(n) ) 
!  do i=1,n
!    read *,z1(i),z2(i),r(i),nu(i)
!  end do
!  call pbvt(n,z1,z2,r,nu,bcdf)
!  do i=1,n
!    print *,z1(i),z2(i),r(i),nu(i),bcdf(i)
!  end do
!  deallocate (z1, z2, r, bcdf, nu)
!  stop
!  end

! gfortran -fpic -c cdfbvt.f pbvt.f90
! gfortran -shared -o pbvt.so pbvt.o cdfbvt.o
subroutine pbvt(n,z1,z2,r,nu,bcdf)
  implicit none
  integer n,nu(n)
  double precision z1(n),z2(n),r(n),bcdf(n)
  integer i,iinfin(2)
  double precision lower(2),upper(2),mvbvt
  iinfin(1)=0; iinfin(2)=0;   ! lower bounds are -inf
  lower(1)=-10.d0; lower(2)=-10.d0;
  do i=1,n
    upper(1)=z1(i); upper(2)=z2(i);
    if(r(i)<-1 .or. r(i)>1. .or. nu(i)<1) then
       bcdf(i)=-1.
       ! later handle r=1, r=-1? with univariate pt()
    else
       bcdf(i)=mvbvt(nu(i),lower,upper,iinfin,r(i))
    endif
  end do
  return
  end


