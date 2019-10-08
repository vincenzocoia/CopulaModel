! Student t 1-factor and 2-factor copula models for discrete item response data

! gcc -c ptbeta.c
! gfortran -o irt12fact irt12fact.f90 qtnorm.f90 ptbeta.o

! sample main
!program irt
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
!  npar=d
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
!  ! 1-factor
!  do j=1,d
!    param(j)=.5d0+j/(2.*d+1.d0)
!  end do
!  print *, param
!  print *, "1-factor"
!  call irt1fact(npar,param,nu1,d,ncat,n,ydata,ucuts,nq,wl,xl,nllk,grad,hess)
!  print *, nllk
!  print "(8f10.6)", grad
!  print "(8f10.6/)", hess
!  deallocate (param,grad,hess)
!  npar=2*d
!  allocate ( param(npar), grad(npar), hess(npar,npar) )
!  ! 2-factor
!  do j=1,d
!    param(j)=.5+j/(2.*d+1.)
!    param(j+d)=.2*j/(d+1.) ! second factor
!  end do
!  print *, param
!  print *, "2-factor"
!  call irt2fact(npar,param,nu1,nu2,d,ncat,n,ydata,ucuts,nq,wl,xl,nllk,grad,hess)
!  print *, nllk
!  print "(8f10.6)", grad
!  print "(8f10.6/)", hess
!  deallocate (xl, wl, ydata, param, grad, hess, ucuts)
!  stop
!  end


! 1-factor model with t, 
! inputs 
!   npar = #parameters,
!   nu1 = constant shape parameter for Student 
!   th = parameter vector (dimension d=npar), -1<th[j]<1
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
subroutine irt1fact(npar,th,nu1,d,ncat,n,ydata,ucuts,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,ncat,nq
  integer ydata(n,d)
  double precision th(npar),nu1,wl(nq),xl(nq),ucuts(ncat-1,d)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,k
  integer, dimension(:), allocatable :: yvec
  double precision, dimension(:), allocatable :: pmf,der1,der2,fval1,integl1,grd
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision, dimension(:,:,:), allocatable :: cond,cond1,cond2
  double precision, dimension(:), allocatable :: tl1
  double precision, dimension(:,:), allocatable :: tcuts
  double precision fval,integl,ww
  double precision qt ! function in f90 code qtnorm.f90

  ! npar=d
  allocate ( yvec(d), pmf(d), der1(d), der2(d), fval1(d), integl1(d), grd(d) )
  allocate ( fval2(d,d), integl2(d,d), hss(d,d) )
  allocate ( cond(nq,ncat,d), cond1(nq,ncat,d), cond2(nq,ncat,d) )
  allocate ( tl1(nq), tcuts(ncat-1,d) )
  ! cond for cond cdf and pmf, cond1 for deriv wrt theta, cond2 for 2nd deriv
   
  !print *, "n=", n, " nitems=", d, " ncat=", ncat
  !print *, "nq=", nq
  !print "(8f10.6/)", wl
  cond=0.d0; cond1=0.d0; cond2=0.d0
  cond(:,ncat,:)=1.d0  ! condcdf is 1 for last category ncat
  do iq=1,nq ! convert to t(nu1) scale
    tl1(iq)=qt(xl(iq),nu1)
  end do
  do j=1,d
    do k=1,ncat-1
      tcuts(k,j)=qt(ucuts(k,j),nu1)
    end do
  end do
  ! get \p copcond(ucuts[k]|v,th[j]) / \p th[j] , k=1,...,ncat ; v=latent
  do iq=1,nq ! loop over quadrature points
    do j=1,d ! loop over items
      do k=1,ncat-1 ! loop of categories
        call irt1derivs(tl1(iq),tcuts(k,j),th(j),nu1,cond(iq,k,j),cond1(iq,k,j),cond2(iq,k,j))
      end do
    end do
  end do
  !print "(8f10.6)", cond(1,:,d)
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
      !if(iq==1 .and. i==1) print "(4f10.6)", pmf(1),der1(1),der2(1)
      !if(iq==1 .and. i==1) print "(8f10.6)", pmf,fval
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
    !if(i==1) print "(8f10.6/)",  integl, integl1
    ! update contribution to negative log-likelihood
    nllk=nllk-log(integl)
    !print *, pmf
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
  deallocate ( tl1, tcuts )
  return
  end

! 2-factor model with t
! inputs 
!   npar = #parameters,
!   param = parameter vector (dimension 2d=npar)
!   nu1 = Student t shape parameter for latent variable 1
!   nu2 = Student t shape parameter for latent variable 2
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
subroutine irt2fact(npar,param,nu1,nu2,d,ncat,n,ydata,ucuts,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,ncat,nq
  integer ydata(n,d)
  double precision param(npar),nu1,nu2,wl(nq),xl(nq),ucuts(ncat-1,d)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,k,iq2
  integer, dimension(:), allocatable :: yvec
  double precision, dimension(:), allocatable :: pmf,tder1,tder2,fval1,integl1,grd
  ! first parameter theta, second parameter gamma, and hence tder1, gder1
  double precision, dimension(:), allocatable :: th,gam,gder1,gder2,tgder
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision, dimension(:,:,:,:), allocatable :: cond,tcond1,tcond2
  double precision, dimension(:,:,:,:), allocatable :: gcond1,gcond2,tgcond2
  double precision, dimension(:), allocatable :: tl1,tl2
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
  allocate ( tl1(nq), tl2(nq), tcuts(ncat-1,d) )
   
  !print *, "n=", n, " nitems=", d, " ncat=", ncat
  !print *, "nq=", nq
  !print "(8f10.6/)", wl
  th(1:d)=param(1:d); gam(1:d)=param((d+1):npar)
  cond=0.d0; tcond1=0.d0; tcond2=0.d0
  gcond1=0.d0; gcond2=0.d0; tgcond2=0.d0
  cond(:,:,ncat,:)=1.d0  ! condcdf is 1 for last category ncat
  do iq=1,nq ! convert to t(nu1), t(nu2) scale
    tl1(iq)=qt(xl(iq),nu1)
    tl2(iq)=qt(xl(iq),nu2)
  end do
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
        tccdf=qt(tccdf,nu2)
        call irt2derivs(tl2(iq2),tccdf,gam(j),nu2,gccdf,gcder1,gcder2,gpdf,gpder1,gpder1u)
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
  deallocate ( tl1, tl2, tcuts )
  return
  end


! function with t condcdf and its deriv wrt rho, order 1 and 2
!  conditioning is u2|u1
! inputs
!   t1,t2 = values with qt(; nu1) transforms in calling routine
!   rho = scalar in (-1,1)
!   nu = degree of freedom parameter
! outputs
!   ccdf = conditional cdf, 
!   cder1 = \p ccdf/\p rho
!   cder2 = \p^2 ccdf/\p rho^2
subroutine irt1derivs(t1,t2,rho,nu,ccdf,cder1,cder2)
  implicit none
  double precision t1,t2,rho,rho2,nu,ccdf,cder1,cder2
  double precision logpi,r2,cr,xt,tem1,tem2,const,gtem1,gtem2
  double precision pt   ! function in C code ptbeta.c
  double precision lgamma   ! function in C code lgamma.c, in gnu fortran?
  logpi=1.1447298858494d0
  rho2 = rho*rho;
  r2 = sqrt((nu+t1*t1)/(nu+1.d0));
  cr = sqrt(1.d0-rho2)
  !xt = (t1-rho*t2)/cr/r2; ! reverse roles of t1,t2 cimpared with ctderivs
  xt = (t2-rho*t1)/cr/r2;
  ccdf = pt(xt,nu+1.d0); 
  tem1 = (t2*rho-t1)/cr**3/r2;
  tem2 = xt*xt/(nu+1.d0);
  const = exp(lgamma(0.5d0*nu+1.d0)-lgamma(0.5d0*nu+0.5)-0.5d0*logpi-0.5d0*log(nu+1.d0));
  cder1 = const*(1.d0+tem2)**(-0.5d0*nu-1.d0)*tem1;
  gtem1 = (t2+2.d0*t2*rho2-3.d0*t1*rho)/cr**5/r2;
  gtem2 = 2.d0*(rho*(t1*t1+t2*t2)-t1*t2*(rho2+1))/(nu+1.d0)/cr**4/r2**2;
  cder2 = -const*(0.5d0*nu+1.d0)*(1.d0+tem2)**(-0.5d0*nu-2.d0)*gtem2*tem1;
  cder2 = cder2+const*(1.d0+tem2)**(-0.5d0*nu-1.d0)*gtem1;
  return
  end

! Function with t condcdf and its deriv wrt rho, order 1 and 2
!  also pdf, pdf1rh= \p pdf/\p rho, pdf1v= \p pdf/ \p u2 (second arg)
! inputs
!   t1,t2 = values with qt(; nu1) transforms in calling routine
!   rho = scalar in (-1,1)
!   nu = degree of freedom parameter
!     t2 is used for the ccdf with latent variable 1
! outputs
!   ccdf = conditional cdf, 
!   cder1 = \p ccdf/\p rho
!   cder2 = \p^2 ccdf/\p rho^2
!   pdf = t pdf
!   pdf1rh= \p pdf/\p rho, 
!   pdf1v= \p pdf/ \p u2
subroutine irt2derivs(t1,t2,rho,nu,ccdf,cder1,cder2,pdf,pdf1rh,pdf1v)
  implicit none
  double precision t1,t2,rho,rho2,nu,ccdf,cder1,cder2,pdf,pdf1rh,pdf1v
  double precision logpi,r2,cr,xt,tem1,tem2,const,gtem1,gtem2
  double precision pt   ! function in C code ptbeta.c
  double precision lpdf,lder1,ltder1v,ltder1t1
  double precision pi,con,lgdif,coef,dt1,dt2,den,den1
  double precision t2sq,dd2,reg,coef3
  double precision lgamma   ! function in C code lgamma.c, in gnu fortran?
  ! from ctderivs
  logpi=1.1447298858494d0
  rho2 = rho*rho;
  r2 = sqrt((nu+t1*t1)/(nu+1.d0));
  cr = sqrt(1.d0-rho2)
  !xt = (t1-rho*t2)/cr/r2;
  xt = (t2-rho*t1)/cr/r2;
  ccdf = pt(xt,nu+1.d0); 
  tem1 = (t2*rho-t1)/cr**3/r2;
  tem2 = xt*xt/(nu+1.d0);
  const = exp(lgamma(0.5d0*nu+1.d0)-lgamma(0.5d0*nu+0.5d0)-0.5*logpi-0.5d0*log(nu+1.d0));
  cder1 = const*(1.d0+tem2)**(-0.5d0*nu-1.d0)*tem1;
  gtem1 = (t2+2.d0*t2*rho2-3.d0*t1*rho)/cr**5/r2;
  gtem2 = 2.d0*(rho*(t1*t1+t2*t2)-t1*t2*(rho2+1.d0))/(nu+1.d0)/cr**4/r2**2;
  cder2 = -const*(0.5d0*nu+1.d0)*(1.d0+tem2)**(-0.5d0*nu-2.d0)*gtem2*tem1;
  cder2 = cder2+const*(1.d0+tem2)**(-0.5d0*nu-1.d0)*gtem1;
  ! pdf and deriv wrt theta, u1
  ! from lt1derivs
  pi=3.14159265358979323846d0; ! pi
  con=0.5d0*log(pi*nu);
  lgdif=lgamma(0.5d0*(nu+1.d0))-lgamma(0.5d0*nu);
  coef = 1.d0-rho2; ! coef2 = coef*coef;
  den = 1.d0+(t1*t1-2.d0*rho*t1*t2+t2*t2)/(nu*coef);
  den1 = 2.d0*rho*(den-1.d0)/coef-2.d0*t1*t2/(nu*coef);
  !den2 = 2*rho*den1/coef+2*(den-1)*(2-coef)/coef2-4*rho*t1*t2/(nu*coef2);
  dt1 = lgdif-con-0.5d0*(nu+1.d0)*log(1.d0+t1*t1/nu);
  dt2 = lgdif-con-0.5d0*(nu+1.d0)*log(1.d0+t2*t2/nu);
  lpdf = -log(2.d0*pi)-0.5d0*log(coef)-0.5d0*(nu+2.d0)*log(den)-dt1-dt2;
  lder1 = rho/coef-0.5d0*(nu+2.d0)*den1/den; ! this is symmetric in t1,t2
  ! from lt2derivs 
  !t1sq = t1*t1;
  !reg = t1-rho*t2;
  !coef3 = den*coef*nu;
  !dd1 = exp(dt1);
  !ltder1t1 = -(nu+2)*reg/coef3 + (nu+1)*t1/(nu+t1sq);
  !ltder1u =ltder1t1/dd1;
  !pdf=exp(lpdf); pdf1rh=pdf*lder1; pdf1u=pdf*ltder1u
  ! switch t1,t2
  t2sq = t2*t2;
  reg = t2-rho*t1;
  coef3 = den*coef*nu;
  dd2 = exp(dt2);
  ltder1t1 = -(nu+2.d0)*reg/coef3 + (nu+1.d0)*t2/(nu+t2sq);
  ltder1v =ltder1t1/dd2;
  pdf=exp(lpdf); pdf1rh=pdf*lder1; pdf1v=pdf*ltder1v
  return
  end


