! t/Gumbel 2-factor copula models for discrete item response data

! gcc -c ptbeta.c
! gfortran -o irtgum2fact irtgum2fact.f90 irgum12fact.f90  irt12fact.f90 qtnorm.f90 ptbeta.o

! sample main
!program irtgum
!  implicit none
!  integer n,d,nq,i,iq,j,k,npar,ncat
!  integer, dimension(:,:), allocatable :: ydata
!  double precision, dimension(:,:), allocatable :: hess,ucuts
!  double precision, dimension(:), allocatable :: xl,wl,param,grad
!  double precision nllk,nu1,nu2
!  read *,nq
!  allocate ( xl(nq),wl(nq) )
!  do iq =1,nq   ! equidistant and equially weighted for testing
!    xl(iq)=iq/(nq+1.d0)
!    wl(iq)=1.d0/nq
!  end do
!  print *, xl
!  print *, wl
!  nu1=5.d0; nu2=5.d0
!  read *,n,d,ncat
!  print *, "n=", n, " nitems=", d, " ncat=", ncat
!  npar=2*d
!  allocate ( ydata(n,d), param(npar), grad(npar), hess(npar,npar), ucuts(ncat-1,d) )
!  do i=1,n
!    read *, ydata(i,:)
!    print *, ydata(i,:)
!  end do
!  do j=1,d
!    do k=1,ncat-1
!      ucuts(k,j)=(j+k)/(ncat+d+1.d0)
!    end do
!    print *, ucuts(:,j)
!  end do
!  ! 2-factor
!  do j=1,d
!    param(j)=.5d0+j/(2.d0*d+1.d0) ! first factor is t
!    !param(j)=1.d0+.5d0*j ! first factor is Gumbel
!    !param(j+d)=.2d0*j/(d+1.d0) ! second factor is t
!    param(j+d)=1.d0+.1d0*j ! second factor is Gumbel
!  end do
!  print *, param
!  print *, "2-factor t/Gumbel"
!  call irtgum2fact(npar,param,nu1,d,ncat,n,ydata,ucuts,nq,wl,xl,nllk,grad,hess)
!  print *, nllk
!  print "(8f10.6)", grad
!  print "(8f10.6/)", hess
!  deallocate (xl, wl, ydata, param, grad, hess, ucuts)
!  stop
!  end


! 2-factor model with t/Gumbel
! Subroutine irt1derivs is in file irt12fact.f90
! Subroutine irgum2derivs is in file irgum12fact.f90 
! inputs 
!   npar = #parameters,
!   param = parameter vector (dimension 2d=npar)
!   nu1 = Student t shape parameter for latent variable 1
!   d = #items
!   ncat = #categories (ordinal)
!   n = sample size
!   ydata = nxd matrix of integer data, each in {0,1,...,ncat-1}
!   ucuts = (ncat-1)xd matrix of cutpoints on Uniform(0,1) scale
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes
! outputs 
!   nllk = negative log-likelihood, 
!   grad = gradient of nllk, 
!   hess = hessian of nllk
subroutine irtgum2fact(npar,param,nu1,d,ncat,n,ydata,ucuts,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,ncat,nq
  integer ydata(n,d)
  double precision param(npar),nu1,wl(nq),xl(nq),ucuts(ncat-1,d)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,k,iq2
  integer, dimension(:), allocatable :: yvec
  double precision, dimension(:), allocatable :: pmf,tder1,tder2,fval1,integl1,grd
  ! first parameter theta, second parameter gamma, and hence tder1, gder1
  double precision, dimension(:), allocatable :: th,gam,gder1,gder2,tgder
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision, dimension(:,:,:,:), allocatable :: cond,tcond1,tcond2
  double precision, dimension(:,:,:,:), allocatable :: gcond1,gcond2,tgcond2
  double precision, dimension(:), allocatable :: tl1
  double precision, dimension(:,:), allocatable :: tcuts
  double precision fval,integl,ww,w2
  ! t (theta) initial for first factor, g (gamma) initial for second factor 
  double precision tccdf,tcder1,tcder2,gccdf,gcder1,gcder2,gpdf,gpder1,gpder1u
  double precision qt ! function in f90 code qtnorm.f90

  ! npar=2*d
  allocate ( yvec(d), pmf(d), tder1(d), tder2(d), fval1(npar), integl1(npar), grd(npar) )
  allocate ( th(d), gam(d), gder1(d),gder2(d),tgder(d) )
  allocate ( fval2(npar,npar), integl2(npar,npar), hss(npar,npar) )
  allocate ( cond(nq,nq,ncat,d), tcond1(nq,nq,ncat,d), tcond2(nq,nq,ncat,d) )
  allocate ( gcond1(nq,nq,ncat,d), gcond2(nq,nq,ncat,d), tgcond2(nq,nq,ncat,d) )
  ! cond for cond cdf and pmf, cond1 for deriv wrt theta, cond2 for 2nd deriv
  allocate ( tl1(nq), tcuts(ncat-1,d) )
   
  th(1:d)=param(1:d); gam(1:d)=param((d+1):npar)
  cond=0.d0; tcond1=0.d0; tcond2=0.d0
  gcond1=0.d0; gcond2=0.d0; tgcond2=0.d0
  cond(:,:,ncat,:)=1.d0  ! condcdf is 1 for last category ncat
  do iq=1,nq ! convert to t(nu1) scale
    tl1(iq)=qt(xl(iq),nu1)
  end do
  !print "(8f10.6/)", tl1
  do j=1,d
    do k=1,ncat-1
      tcuts(k,j)=qt(ucuts(k,j),nu1)
    end do
  end do

  ! get \p copcond(ucuts[k]|v,th[j]) / \p th[j] , k=1,...,ncat ; v=latent
  ! below might not be most efficient
  do iq=1,nq ! loop over quadrature points
  do iq2=1,nq ! loop over quadrature points
    do j=1,d ! loop over items
      do k=1,ncat-1 ! loop of categories
        call irt1derivs(tl1(iq),tcuts(k,j),th(j),nu1,tccdf,tcder1,tcder2)
        !if(iq==1 .and. iq2==1 .and. k==1) print "(4f11.6)", th(j),tccdf,tcder1,tcder2
        call irgum2derivs(xl(iq2),tccdf,gam(j),gccdf,gcder1,gcder2,gpdf,gpder1,gpder1u)
        !if(iq==1 .and. iq2==1 .and. k==1) print "(7f11.6)", gam(j),gccdf,gcder1,gcder2,gpdf,gpder1,gpder1u
        cond(iq,iq2,k,j)=gccdf
        gcond1(iq,iq2,k,j)=gcder1
        gcond2(iq,iq2,k,j)=gcder2
        tcond1(iq,iq2,k,j)=gpdf*tcder1
        tcond2(iq,iq2,k,j)=gpdf*tcder2+gpder1u*tcder1*tcder1
        tgcond2(iq,iq2,k,j)=gpder1*tcder1
      end do
    end do
  end do
  end do
  ! differencing
  do iq=1,nq 
  do iq2=1,nq 
    do j=1,d
      do k=ncat,2,-1
        cond(iq,iq2,k,j)=cond(iq,iq2,k,j)-cond(iq,iq2,k-1,j)
        tcond1(iq,iq2,k,j)=tcond1(iq,iq2,k,j)-tcond1(iq,iq2,k-1,j)
        tcond2(iq,iq2,k,j)=tcond2(iq,iq2,k,j)-tcond2(iq,iq2,k-1,j)
        gcond1(iq,iq2,k,j)=gcond1(iq,iq2,k,j)-gcond1(iq,iq2,k-1,j)
        gcond2(iq,iq2,k,j)=gcond2(iq,iq2,k,j)-gcond2(iq,iq2,k-1,j)
        tgcond2(iq,iq2,k,j)=tgcond2(iq,iq2,k,j)-tgcond2(iq,iq2,k-1,j)
      end do
    end do
  end do
  end do
  !print "(8f10.6)", cond(1,1,:,d)
  !print "(8f10.6)", tcond1(1,1,:,d)
  !print "(8f10.6)", gcond1(1,1,:,d)
  ! looping through data
  nllk=0.d0; grad=0.d0; hess=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    yvec=ydata(i,:)+1
    !if(i==1) print "(8i3)", yvec
    integl=0.d0; integl1=0.d0; integl2=0.d0; ! initialize for integrals
    do iq=1,nq 
    do iq2=1,nq 
      ! update integrand, and derivs wrt copula parameters 
      do j=1,d
        pmf(j)=cond(iq,iq2,yvec(j),j)
        tder1(j)=tcond1(iq,iq2,yvec(j),j)  ! deriv wrt th
        tder2(j)=tcond2(iq,iq2,yvec(j),j)
        gder1(j)=gcond1(iq,iq2,yvec(j),j)  ! deriv wrt gam
        gder2(j)=gcond2(iq,iq2,yvec(j),j)
        tgder(j)=tgcond2(iq,iq2,yvec(j),j)
      end do
      fval=product(pmf)
      ! fval1: vector of partial deriv wrt theta[j], j=1,...,d
      ! fval2: matrix of 2nd partial deriv wrt theta[j], theta[j2], 
      do j=1,d
        fval1(j)=fval*tder1(j)/pmf(j)
        fval1(j+d)=fval*gder1(j)/pmf(j)
        fval2(j,j)=fval*tder2(j)/pmf(j)
        fval2(j+d,j+d)=fval*gder2(j)/pmf(j)
        ! cross th , gam terms for fval2
        fval2(j,j+d)=fval* tgder(j)/pmf(j) 
        fval2(j+d,j) = fval2(j,j+d)
      end do
      do j=1,d
        do j2=1,d
          if(j2 /= j) then
            fval2(j,j2)=fval1(j)*tder1(j2)/pmf(j2)  
            fval2(j+d,j2+d)=fval1(j+d)*gder1(j2)/pmf(j2)
            fval2(j+d,j2)=fval1(j+d)*tder1(j2)/pmf(j2) 
            fval2(j,j2+d)=fval1(j)*gder1(j2)/pmf(j2) 
          endif
        end do
      end do
      ! update quadrature loops
      ww=wl(iq); w2=wl(iq2)
      integl=integl+fval*ww*w2
      integl1=integl1+fval1*ww*w2
      integl2=integl2+fval2*ww*w2
    end do
    end do
    ! update contribution to negative log-likelihood
    nllk=nllk-log(integl)
    grd=integl1/integl; grad=grad-grd;
    do j=1,npar
      do j2=1,npar
        hss(j,j2)=integl2(j,j2)/integl - grd(j)*grd(j2); 
      end do
    end do
    hess=hess-hss
  end do
  deallocate (yvec, pmf, tder1, tder2, fval1, fval2, integl1, integl2, grd, hss)
  deallocate (cond, tcond1, tcond2, gcond1, gcond2, tgcond2)
  deallocate (th, gam, gder1, gder2, tgder)
  deallocate ( tl1, tcuts )
  return
  end

