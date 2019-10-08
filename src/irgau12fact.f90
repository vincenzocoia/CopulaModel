! Gaussian 1-factor and 2-factor copula models for discrete item response data

! gfortran -o irgau12fact irgau12fact.f90 pnorm.f90 qtnorm.f90

! sample main
!program irgau
!  implicit none
!  integer n,d,nq,i,iq,j,k,npar,ncat
!  integer, dimension(:,:), allocatable :: ydata
!  double precision, dimension(:,:), allocatable :: hess,ucuts
!  double precision, dimension(:), allocatable :: xl,wl,param,grad
!  double precision nllk
!  read *,nq
!  allocate ( xl(nq),wl(nq) )
!  do iq =1,nq   ! equidistant and equially weighted for testing
!    xl(iq)=iq/(nq+1.d0)
!    wl(iq)=1.d0/nq
!  end do
!  read *,n,d,ncat
!  npar=d
!  allocate ( ydata(n,d), param(npar), grad(npar), hess(npar,npar), ucuts(ncat-1,d) )
!  do i=1,n
!    read *, ydata(i,:)
!  end do
!  do j=1,d
!    do k=1,ncat-1
!      ucuts(k,j)=(j+k)/(ncat+d+1.d0)
!    end do
!  end do
!  !print *, ucuts
!  do j=1,npar
!    param(j)=.5d0+j/(2.d0*d+1.d0)
!  end do
!  print *, "1-factor"
!  call irgau1fact(npar,param,d,ncat,n,ydata,ucuts,nq,wl,xl,nllk,grad,hess)
!  print *, nllk
!  print "(8f10.6)", grad
!  print "(8f10.6/)", hess
!  deallocate (param,grad,hess)
!  npar=2*d
!  allocate ( param(npar), grad(npar), hess(npar,npar) )
!  ! 2-factor
!  do j=1,d
!    param(j)=.5d0+j/(2.d0*d+1.d0)
!    param(j+d)=.2d0*j/(d+1.d0) ! second factor
!  end do
!  print *, param
!  print *, "2-factor"
!  call irgau2fact(npar,param,d,ncat,n,ydata,ucuts,nq,wl,xl,nllk,grad,hess)
!  print *, nllk
!  print "(8f10.6)", grad
!  print "(8f10.6/)", hess
!  deallocate (xl, wl, ydata, param, grad, hess, ucuts)
!  stop
!  end

! 1-factor model with Gaussian, 
! inputs 
!   npar = #parameters,
!   th = rho = parameter vector (dimension d=npar),  -1<th[j]<1
!   d = #items
!   ncat = #categories (ordinal)
!   n = sample size
!   ydata = nxd matrix of integer data, each in {0,1,...,ncat-1}
!   ucuts = (ncat-1)xd matrix of cutpoints on N(0,1) scale [not U(0,1)]
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes gl$nodes
! outputs 
!   nllk = negative log-likelihood, 
!   grad = gradient of nllk, 
!   hess = hessian of nllk
subroutine irgau1fact(npar,th,d,ncat,n,ydata,ucuts,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,ncat,nq
  integer ydata(n,d)
  double precision th(npar),wl(nq),xl(nq),ucuts(ncat-1,d)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,k
  integer, dimension(:), allocatable :: yvec
  double precision, dimension(:), allocatable :: pmf,der1,der2,fval1,integl1,grd
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision, dimension(:,:,:), allocatable :: cond,cond1,cond2
  double precision, dimension(:), allocatable :: zl
  double precision, dimension(:,:), allocatable :: zcuts
  double precision fval,integl,ww
  double precision qnorms  ! in qtnorm.f90

  ! npar=d
  allocate ( yvec(d), pmf(d), der1(d), der2(d), fval1(d), integl1(d), grd(d) )
  allocate ( fval2(d,d), integl2(d,d), hss(d,d) )
  allocate ( cond(nq,ncat,d), cond1(nq,ncat,d), cond2(nq,ncat,d) )
  allocate ( zl(nq), zcuts(ncat-1,d) )
  ! cond for cond cdf and pmf, cond1 for deriv wrt theta, cond2 for 2nd deriv
   
  !print *, "n=", n, " nitems=", d, " ncat=", ncat
  !print *, ucuts
  !print *, zl
  cond=0.d0; cond1=0.d0; cond2=0.d0
  cond(:,ncat,:)=1.d0  ! condcdf is 1 for last category ncat
  do iq=1,nq ! convert to N(0,1) scale
    zl(iq)=qnorms(xl(iq))
  end do
  do j=1,d
    do k=1,ncat-1
      zcuts(k,j)=qnorms(ucuts(k,j))
    end do
  end do

  ! get \p copcond(ucuts[k]|v,th[j]) / \p th[j] , k=1,...,ncat ; v=latent
  do iq=1,nq ! loop over quadrature points
    do j=1,d ! loop over items
      do k=1,ncat-1 ! loop of categories
        call irgau1derivs(zl(iq),zcuts(k,j),th(j),cond(iq,k,j),cond1(iq,k,j),cond2(iq,k,j))
      end do
    end do
  end do
  ! differencing
  do iq=1,nq 
    do j=1,d
      do k=ncat,2,-1
        cond(iq,k,j)=cond(iq,k,j)-cond(iq,k-1,j)
        cond1(iq,k,j)=cond1(iq,k,j)-cond1(iq,k-1,j)
        cond2(iq,k,j)=cond2(iq,k,j)-cond2(iq,k-1,j)
      end do
    end do
  end do
  !print "(8f10.6)", cond(1,:,d)
  !print "(8f10.6)", cond1(1,:,d)
  !print "(8f10.6)", cond2(1,:,d)
  ! looping through data
  nllk=0.d0; grad=0.d0; hess=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    yvec=ydata(i,:)+1
    !if(i==1) print "(8i3)", yvec
    integl=0.d0; integl1=0.d0; integl2=0.d0; ! initialize for integrals
    do iq=1,nq 
      ! update integrand, and derivs wrt copula parameters 
      do j=1,d
        pmf(j)=cond(iq,yvec(j),j)
        der1(j)=cond1(iq,yvec(j),j)
        der2(j)=cond2(iq,yvec(j),j)
      end do
      fval=product(pmf)
      ! fval1: vector of partial deriv wrt theta[j], j=1,...,d
      ! fval2: matrix of 2nd partial deriv wrt theta[j], theta[j2], 
      do j=1,d
        fval1(j)=fval*der1(j)/pmf(j)
        fval2(j,j)=fval*der2(j)/pmf(j)
        do j2=1,d
          if(j2 /= j) fval2(j,j2)=fval1(j)*der1(j2)/pmf(j2)
        end do
      end do
      ! update quadrature loops
      ww=wl(iq)
      integl=integl+fval*ww
      integl1=integl1+fval1*ww
      integl2=integl2+fval2*ww
    end do
    ! update contribution to negative log-likelihood
    nllk=nllk-log(integl)
    grd=integl1/integl; grad=grad-grd;
    do j=1,d
      do j2=1,d
        hss(j,j2)=integl2(j,j2)/integl - grd(j)*grd(j2);
      end do
    end do
    hess=hess-hss
  end do
  deallocate (yvec, pmf, der1, der2, fval1, fval2, integl1, integl2, grd, hss)
  deallocate (cond, cond1, cond2)
  deallocate ( zl, zcuts )
  return
  end

! 2-factor model with Gaussian
! inputs 
!   npar = #parameters,
!   param = parameter vector (dimension 2d=npar)
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
subroutine irgau2fact(npar,param,d,ncat,n,ydata,ucuts,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,ncat,nq
  integer ydata(n,d)
  double precision param(npar),wl(nq),xl(nq),ucuts(ncat-1,d)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,k,iq2
  integer, dimension(:), allocatable :: yvec
  double precision, dimension(:), allocatable :: pmf,tder1,tder2,fval1,integl1,grd
  ! first parameter theta, second parameter gamma, and hence tder1, gder1
  double precision, dimension(:), allocatable :: th,gam,gder1,gder2,tgder
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision, dimension(:,:,:,:), allocatable :: cond,tcond1,tcond2
  double precision, dimension(:,:,:,:), allocatable :: gcond1,gcond2,tgcond2
  double precision, dimension(:), allocatable :: zl
  double precision, dimension(:,:), allocatable :: zcuts
  double precision fval,integl,ww,w2
  ! t (theta) initial for first factor, g (gamma) initial for second factor 
  double precision tccdf,tcder1,tcder2,gccdf,gcder1,gcder2,gpdf,gpder1,gpder1u
  double precision qnorms ! function in f90 code qtnorm.f90

  ! npar=2*d
  allocate ( yvec(d), pmf(d), tder1(d), tder2(d), fval1(npar), integl1(npar), grd(npar) )
  allocate ( th(d), gam(d), gder1(d),gder2(d),tgder(d) )
  allocate ( fval2(npar,npar), integl2(npar,npar), hss(npar,npar) )
  allocate ( cond(nq,nq,ncat,d), tcond1(nq,nq,ncat,d), tcond2(nq,nq,ncat,d) )
  allocate ( gcond1(nq,nq,ncat,d), gcond2(nq,nq,ncat,d), tgcond2(nq,nq,ncat,d) )
  ! cond for cond cdf and pmf, cond1 for deriv wrt theta, cond2 for 2nd deriv
  allocate ( zl(nq), zcuts(ncat-1,d) )
   
  th(1:d)=param(1:d); gam(1:d)=param((d+1):npar)
  cond=0.d0; tcond1=0.d0; tcond2=0.d0
  gcond1=0.d0; gcond2=0.d0; tgcond2=0.d0
  cond(:,:,ncat,:)=1.d0  ! condcdf is 1 for last category ncat
  do iq=1,nq ! convert to N(0,1) scale
    zl(iq)=qnorms(xl(iq))
  end do
  do j=1,d
    do k=1,ncat-1
      zcuts(k,j)=qnorms(ucuts(k,j))
    end do
  end do

  ! get \p copcond(ucuts[k]|v,th[j]) / \p th[j] , k=1,...,ncat ; v=latent
  ! below might not be most efficient
  do iq=1,nq ! loop over quadrature points
  do iq2=1,nq ! loop over quadrature points
    do j=1,d ! loop over items
      do k=1,ncat-1 ! loop of categories
        call irgau1derivs(zl(iq),zcuts(k,j),th(j),tccdf,tcder1,tcder2)
        !if(iq==1 .and. iq2==1 .and. k==1) print "(4f11.6)", th(j),tccdf,tcder1,tcder2
        tccdf=qnorms(tccdf)
        call irgau2derivs(zl(iq2),tccdf,gam(j),gccdf,gcder1,gcder2,gpdf,gpder1,gpder1u)
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
  deallocate ( zl, zcuts )
  return
  end


! Function with Gaussian condcdf and its deriv wrt rho, order 1 and 2
!  conditioning is z2|z1
! inputs
!   z1,z2 = values in N(0,1) scale in (-oo,oo)
!   rho = scalar in (-1,1)
! outputs
!   ccdf = conditional cdf, 
!   cder1 = \p ccdf/\p rho
!   cder2 = \p^2 ccdf/\p rho^2
subroutine irgau1derivs(z1,z2,rho,ccdf,cder1,cder2)
  implicit none
  double precision z1,z2,rho,ccdf,cder1,cder2
  double precision pnorms,dnorms
  double precision r1,r1sq,zcond,cpdf,tem
  r1=1.d0-rho*rho; r1sq=sqrt(r1)
  zcond=(z2-rho*z1)/r1sq
  ccdf=pnorms(zcond); cpdf=dnorms(zcond)
  cder1=cpdf*(rho*z2-z1)/r1/r1sq
  tem=((rho*z2-z1)**2) *(rho*z1-z2) + (z2+2.d0*z2*rho*rho-3.d0*rho*z1)*r1
  cder2=cpdf*tem/(r1**3)/r1sq
  return
  end

! Function with Gaussian condcdf and its deriv wrt rho, order 1 and 2
!   also pdf, pdf1rh= \p pdf/\p rho, pdf1v= \p pdf/ \p u2 (second arg)
! inputs
!   z1,z2 = values in N(0,1) scale in (-oo,oo)
!   rho = scalar in (-1,1)
! outputs
!   ccdf = conditional cdf, 
!   cder1 = \p ccdf/\p rho
!   cder2 = \p^2 ccdf/\p rho^2
!   pdf = Gaussian pdf
!   pdf1rh= \p pdf/\p rho, 
!   pdf1v= \p pdf/ \p u2
subroutine irgau2derivs(z1,z2,rho,ccdf,cder1,cder2,pdf,pdf1rh,pdf1v)
  implicit none
  double precision z1,z2,rho,ccdf,cder1,cder2,cpdf,pdf,pdf1rh,pdf1v
  double precision pnorms,dnorms
  double precision r1,r1sq,zcond,tem
  double precision qf,lder1rh,lder1v
  r1=1.d0-rho*rho; r1sq=sqrt(r1)
  zcond=(z2-rho*z1)/r1sq
  ccdf=pnorms(zcond); cpdf=dnorms(zcond)
  cder1=cpdf*(rho*z2-z1)/r1/r1sq
  tem=((rho*z2-z1)**2) *(rho*z1-z2) + (z2+2.d0*z2*rho*rho-3.d0*rho*z1)*r1
  cder2=cpdf*tem/(r1**3)/r1sq
  ! pdf, lpdf1rh= \p lpdf/\p rho, lpdf1v= \p lpdf/ \p u2
  qf=z1*z1+z2*z2-2.d0*rho*z1*z2;
  lder1rh=rho/r1 - rho*qf/r1**2 + z1*z2/r1
  lder1v= (-z2+rho*z1)/r1+z2; lder1v=lder1v/dnorms(z2)
  ! next line is copula pdf
  pdf=exp(-.5d0*(qf/r1-z1*z1-z2*z2))/sqrt(r1) 
  pdf1rh=pdf*lder1rh; pdf1v=pdf*lder1v
  return
  end

