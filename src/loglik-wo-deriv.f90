! negative log-likelihoods without derivatives

! sample main program to check cond option
! gfortran -o llik loglik-wo-deriv.f90 ../src/qtnorm.f90 ptbeta.o lgamma.o
!program llkwoderiv
!  implicit none
!  integer n,d,nq,i,iq,j
!  double precision, dimension(:), allocatable :: xl,wl,th
!  double precision, dimension(:,:), allocatable :: udata,param
!  double precision nu,lpdf,ccdf
!  read *,nq
!  allocate ( xl(nq),wl(nq) )
!  do iq =1,nq   ! equidistant and equally weighted for testing
!    xl(iq)=iq/(nq+1.d0)
!    wl(iq)=1.d0/nq
!  end do
!  read *,n,d
!  allocate ( udata(n,d), th(d), param(2,d) )
!  do i=1,n
!    read *, udata(i,:)
!  end do
!  do j=1,d
!    th(j)=1.d0+.5d0*j
!    param(1,j)=.1d0*j
!    param(2,j)=1.d0+.1d0*j
!  end do
!  nu=3.d0
!  do j=1,d
!    call lfrk(udata(1,j),xl(1),th(j),1,lpdf,ccdf)
!    print *,"frk",j,lpdf,ccdf
!    call lfrk(udata(1,j),xl(1),th(j),0,lpdf,ccdf)
!    print *, lpdf,ccdf
!    call lgum(udata(1,j),xl(2),th(j),1,lpdf,ccdf)
!    print *, "gum",j,lpdf,ccdf
!    call lgum(udata(1,j),xl(2),th(j),0,lpdf,ccdf)
!    print *, lpdf,ccdf
!    call lbb1(udata(1,j),xl(2),param(:,j),1,lpdf,ccdf)
!    print *, "bb1",j,lpdf,ccdf
!    call lbb1(udata(1,j),xl(2),param(:,j),0,lpdf,ccdf)
!    print *, lpdf,ccdf
!    call lt(udata(1,j),xl(1),param(1,j),nu,1,lpdf,ccdf)
!    print *, "t-3",j,lpdf,ccdf
!    call lt(udata(1,j),xl(1),param(1,j),nu,0,lpdf,ccdf)
!    print *, lpdf,ccdf
!  end do
!  deallocate (xl, wl, udata, param)
!  stop
!  end


! Frank 1-factor copula
! inputs 
!   npar = #parameters,
!   th = parameter vector (dimension d=npar), -oo<th[j]<oo
!   d = #variables
!   n = sample size
!   udata = nxd matrix of uniform scores
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes
! output 
!   nllk = negative log-likelihood
subroutine lfrk1lik(npar,th,d,n,udata,nq,wl,xl,nllk)
  implicit none
  integer npar,d,n,nq 
  double precision nllk, th(npar),udata(n,d),wl(nq),xl(nq)
  integer i,iq,j
  double precision, dimension(:), allocatable :: uvec,lpdf,ccdf
  double precision fval,integl,ww

  ! npar=d
  allocate ( uvec(d), lpdf(d), ccdf(d) )
  nllk=0.d0;  ! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    integl=0.d0; ! initialize for integrals
    do iq=1,nq ! loop over quadrature points
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        call lfrk(uvec(j),xl(iq),th(j),0,lpdf(j),ccdf(j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=exp(sum(lpdf))
      ! update quadrature loops
      ww=wl(iq)
      integl=integl+fval*ww
    end do
    ! update contribution to negative log-likelihood
    nllk=nllk-log(integl)
  end do
  deallocate (uvec, lpdf, ccdf)
  return
  end

! Gumbel 1-factor copula
! inputs 
!   npar = #parameters,
!   th = parameter vector (dimension d=npar), th[j]>1
!   d = #variables
!   n = sample size
!   udata = nxd matrix of uniform scores
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes
! output 
!   nllk = negative log-likelihood
subroutine lgum1lik(npar,th,d,n,udata,nq,wl,xl,nllk)
  implicit none
  integer npar,d,n,nq 
  double precision nllk, th(npar),udata(n,d),wl(nq),xl(nq)
  integer i,iq,j
  double precision, dimension(:), allocatable :: uvec,lpdf,ccdf
  double precision fval,integl,ww

  ! npar=d
  allocate ( uvec(d), lpdf(d), ccdf(d) )
  nllk=0.d0;  ! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    integl=0.d0; ! initialize for integrals
    do iq=1,nq ! loop over quadrature points
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        call lgum(uvec(j),xl(iq),th(j),0,lpdf(j),ccdf(j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=exp(sum(lpdf))
      ! update quadrature loops
      ww=wl(iq)
      integl=integl+fval*ww
    end do
    ! update contribution to negative log-likelihood
    nllk=nllk-log(integl)
  end do
  deallocate (uvec, lpdf, ccdf)
  return
  end

! Student t 1-factor copula 
! inputs 
!   npar = #parameters,
!   th = parameter vector (dimension d=npar), -1<th[j]<1
!   nu = degree of freedom >0
!   d = #variables
!   n = sample size
!   tdata = nxd matrix of transformed uniform scores, qt(udata,nu)
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   tl = vector of transformed quadrature nodes, qt(xl,nu)
! output 
!   nllk = negative log-likelihood
subroutine lt1lik(npar,th,nu,d,n,tdata,nq,wl,tl,nllk)
  implicit none
  integer npar,d,n,nq 
  double precision nllk, th(npar),nu,tdata(n,d),wl(nq),tl(nq)
  integer i,iq,j
  double precision, dimension(:), allocatable :: tvec,lpdf,ccdf
  double precision fval,integl,ww

  ! npar=d
  allocate ( tvec(d), lpdf(d), ccdf(d) )
  nllk=0.d0;  ! initialize
  do i=1,n ! loop over rows of data set
    tvec=tdata(i,:)
    integl=0.d0; ! initialize for integrals
    do iq=1,nq ! loop over quadrature points
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        call lt(tvec(j),tl(iq),th(j),nu,0,lpdf(j),ccdf(j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=exp(sum(lpdf))
      ! update quadrature loops
      ww=wl(iq)
      integl=integl+fval*ww
    end do
    ! update contribution to negative log-likelihood
    nllk=nllk-log(integl)
  end do
  deallocate (tvec, lpdf, ccdf)
  return
  end
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! BB1 1-factor copula 
! inputs 
!   npar = #parameters,
!   param0 = parameter vector (dimension 2d=npar), theta>0, delta>1
!   d = #variables
!   n = sample size
!   udata = nxd matrix of uniform scores
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes
! output 
!   nllk = negative log-likelihood
subroutine lbb11lik(npar,param0,d,n,udata,nq,wl,xl,nllk)
  implicit none
  integer npar,d,n,nq 
  double precision nllk, param0(npar),udata(n,d),wl(nq),xl(nq)
  integer i,iq,j
  double precision, dimension(:), allocatable :: uvec,lpdf,ccdf
  double precision fval,integl,ww
  double precision, dimension(:,:), allocatable :: param

  ! npar=d*2
  allocate ( uvec(d), lpdf(d), ccdf(d) )
  allocate( param(2,d) )
  param=reshape(param0,(/2,d/))
  nllk=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    integl=0.d0; ! initialize for integrals
    do iq=1,nq ! loop over quadrature points
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        call lbb1(uvec(j),xl(iq),param(:,j),0,lpdf(j),ccdf(j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=exp(sum(lpdf))
      ! update quadrature loops
      ww=wl(iq)
      integl=integl+fval*ww
    end do
    ! update contribution to negative log-likelihood
    nllk=nllk-log(integl)
  end do
  deallocate (uvec, lpdf, ccdf)
  return
  end

! Frank 2-factor copula
! inputs 
!   npar = #parameters,
!   param = parameter vector (dimension 2d=npar), -oo<param[j]<oo
!   d = #variables
!   n = sample size
!   udata = nxd matrix of uniform scores
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes
! output 
!   nllk = negative log-likelihood
subroutine lfrk2lik(npar,param,d,n,udata,nq,wl,xl,nllk)
  implicit none
  integer npar,d,n,nq 
  double precision nllk, param(npar),udata(n,d),wl(nq),xl(nq)
  integer i,iq,j,iq2
  double precision, dimension(:), allocatable :: uvec,tlpdf
  double precision fval,integl,ww,w2,epsvl
  double precision, dimension(:), allocatable :: th,gam,glpdf
  double precision, dimension(:), allocatable :: tccdf,gccdf

  ! npar=2*d
  allocate ( uvec(d), tlpdf(d) )
  allocate ( th(d), gam(d), glpdf(d) ) 
  allocate ( tccdf(d), gccdf(d) )
  th(1:d)=param(1:d); gam(1:d)=param((d+1):npar)
  nllk=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    integl=0.d0; ! initialize for integrals
    !print *, integl1
    do iq=1,nq ! loop over quadrature points
    do iq2=1,nq ! loop over quadrature points
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        ! factor 1 lpdf and ccdf
        call lfrk(uvec(j),xl(iq),th(j),1,tlpdf(j),tccdf(j))
        ! factor 2 lpdf
        !epsvl = 0.00000001  
        epsvl = 1.d-5 !to prevent the first argument of copula to be close to zero or one
        if(tccdf(j)<epsvl)   tccdf(j)=epsvl ! helps to avoid "singular hessian" error 
        if(tccdf(j)>1-epsvl) tccdf(j)=1-epsvl   !
        call lfrk(tccdf(j),xl(iq2),gam(j),0,glpdf(j),gccdf(j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=exp(sum(tlpdf)+sum(glpdf))
      ! update quadrature loops
      ww=wl(iq); w2=wl(iq2)
      integl=integl+fval*ww*w2
    end do
    end do
    ! update contribution to negative log-likelihood
    nllk=nllk-log(integl)
    end do
  deallocate (uvec, tlpdf)
  deallocate (th, gam, glpdf)
  deallocate (tccdf, gccdf)
  return
  end


! Gumbel 2-factor copula
! inputs 
!   npar = #parameters,
!   param = parameter vector (dimension 2d=npar), param[j]>1
!   d = #variables
!   n = sample size
!   udata = nxd matrix of uniform scores
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes
! output 
!   nllk = negative log-likelihood
subroutine lgum2lik(npar,param,d,n,udata,nq,wl,xl,nllk)
  implicit none
  integer npar,d,n,nq 
  double precision nllk, param(npar),udata(n,d),wl(nq),xl(nq)
  integer i,iq,j,iq2
  double precision, dimension(:), allocatable :: uvec,tlpdf
  double precision fval,integl,ww,w2,epsvl
  double precision, dimension(:), allocatable :: th,gam,glpdf
  double precision, dimension(:), allocatable :: tccdf,gccdf

  ! npar=2*d
  allocate ( uvec(d), tlpdf(d) )
  allocate ( th(d), gam(d), glpdf(d) ) 
  allocate ( tccdf(d), gccdf(d) )
  th(1:d)=param(1:d); gam(1:d)=param((d+1):npar)
  nllk=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    integl=0.d0;  ! initialize for integrals
    !print *, integl1
    do iq=1,nq ! loop over quadrature points
    do iq2=1,nq ! loop over quadrature points
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        ! factor 1 lpdf and ccdf
        call lgum(uvec(j),xl(iq),th(j),1,tlpdf(j),tccdf(j))
        ! factor 2 lpdf
        epsvl = 1.d-5                            !to prevent the first argument of copula to be close to zero or one
        if(tccdf(j)<epsvl)   tccdf(j)=epsvl     !helps to avoid "singular hessian" error for Gumbel copula
        if(tccdf(j)>1-epsvl) tccdf(j)=1-epsvl   !
        call lgum(tccdf(j),xl(iq2),gam(j),0,glpdf(j),gccdf(j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=exp(sum(tlpdf)+sum(glpdf))
      ! update quadrature loops
      ww=wl(iq); w2=wl(iq2)
      integl=integl+fval*ww*w2
    end do
    end do
    ! update contribution to negative log-likelihood
    nllk=nllk-log(integl)
    end do
  deallocate (uvec, tlpdf)
  deallocate (th, gam, glpdf)
  deallocate (tccdf, gccdf)
  return
  end

! Student 2-factor copula
! inputs 
!   npar = #parameters,
!   param = parameter vector (dimension 2d=npar), -1<param[j]<1
!   nu1 = degree of freedom >0
!   nu2 = degree of freedom >0
!   d = #variables
!   n = sample size
!   tdata = nxd matrix of transformed uniform scores : t(nu1) scale
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   tl1 = vector of transformed quadrature nodes wrt nu1
!   tl2 = vector of transformed quadrature nodes wrt nu2
! output 
!   nllk = negative log-likelihood
subroutine lt2lik(npar,param,nu1,nu2,d,n,tdata,nq,wl,tl1,tl2,nllk)
  implicit none
  integer npar,d,n,nq
  double precision nllk, param(npar),nu1,nu2,tdata(n,d),wl(nq),tl1(nq),tl2(nq)
  integer i,iq,j,iq2 
  double precision, dimension(:), allocatable :: tvec,tlpdf
  double precision fval,integl,ww,w2,qtccdf
  double precision, dimension(:), allocatable :: th,gam,glpdf
  double precision, dimension(:), allocatable :: tccdf,gccdf
  double precision qt   ! function in f90 code qtnorm.f90

  ! npar=2*d
  allocate ( tvec(d), tlpdf(d) )
  allocate ( th(d), gam(d), glpdf(d) ) 
  allocate ( tccdf(d), gccdf(d) )
  th(1:d)=param(1:d); gam(1:d)=param((d+1):npar)
  nllk=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    tvec=tdata(i,:)
    integl=0.d0; ! initialize for integrals
    do iq=1,nq ! loop over quadrature points
    do iq2=1,nq ! loop over quadrature points
      ! get \p copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        ! factor 1 lpdf and ccdf
        call lt(tvec(j),tl1(iq),th(j),nu1,1,tlpdf(j),tccdf(j))
        ! print "(3f10.6)", pdfv(j),der1(j),der2(j)
        ! factor 2 lpdf 
        qtccdf = qt(tccdf(j),nu2); ! **Student quantile is needed here**  
        call lt(qtccdf,tl2(iq2),gam(j),nu2,0,glpdf(j),gccdf(j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=exp(sum(tlpdf)+sum(glpdf))
      ! update quadrature loops
      ww=wl(iq); w2=wl(iq2)
      integl=integl+fval*ww*w2
    end do
    end do
    ! update contribution to negative log-likelihood
    nllk=nllk-log(integl)
    end do
  deallocate (tvec, tlpdf)
  deallocate (th, gam, glpdf)
  deallocate (tccdf, gccdf)
  return
  end

! BB1/Frank 2-factor copula
! inputs 
!   npar = #parameters,
!   param0 = parameter vector (dimension 3d=npar), theta>0, delta>1, -oo<gam<oo
!   d = #variables
!   n = sample size
!   udata = nxd matrix of uniform scores
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes
! output 
!   nllk = negative log-likelihood
subroutine lbb1frklik(npar,param0,d,n,udata,nq,wl,xl,nllk)
  implicit none
  integer npar,d,n,nq
  double precision param0(npar),udata(n,d),wl(nq),xl(nq)
  double precision nllk
  integer i,iq,j,d2,iq2 
  double precision, dimension(:), allocatable :: uvec,tlpdf,tccdf
  double precision, dimension(:), allocatable :: gam,glpdf,gccdf
  double precision, dimension(:,:), allocatable :: param
  double precision fval,integl,ww,w2
  
  !npar=d*3
  allocate ( uvec(d), tlpdf(d) )
  allocate ( gam(d), glpdf(d) ) 
  allocate ( tccdf(d), gccdf(d) )
  allocate ( param(2,d) )
  d2=2*d
  param=reshape(param0(1:d2),(/2,d/))
  gam(1:d)=param0((d2+1):npar)
  nllk=0.d0; ! initialize
  do i=1,n   ! loop over rows of data set
    uvec=udata(i,:)
    integl=0.d0; ! initialize for integrals
    do iq=1,nq   ! loop over quadrature points
    do iq2=1,nq  ! loop over quadrature points
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        ! factor 1 lpdf and ccdf
        call lbb1(uvec(j),xl(iq),param(:,j),1,tlpdf(j),tccdf(j))
        ! factor 2 lpdf 
        if (tccdf(j) < 1.d-5) tccdf(j) = 1.d-5         !!! to avoid overflow
        if (tccdf(j) > 1.d0-1.d-5) tccdf(j) = 1.d0-1.d-5  !!! to avoid overflow
        call lbb1(tccdf(j),xl(iq2),gam(j),0,glpdf(j),gccdf(j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=exp(sum(tlpdf)+sum(glpdf))
      ! update quadrature loops
      ww=wl(iq); w2=wl(iq2)
      integl=integl+fval*ww*w2
    end do
    end do
    ! update contribution to negative log-likelihood
    nllk=nllk-log(integl)
    end do
  deallocate (uvec, tlpdf)
  deallocate (param, gam, glpdf)
  deallocate (tccdf, gccdf)
  return
  end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! function with Frank lpdf and ccdf for the factor copulas 
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (-oo,oo) but could be underflow problems for cpar>35
!   icond = 1 for computing conditional cdf, =0 for not
! outputs 
!   lpdf = log pdf 
!   ccdf = conditional pdf if icond=1
subroutine lfrk(u1,u2,cpar,icond,lpdf,ccdf)
  implicit none
  double precision u1,u2,cpar,ccdf,lpdf
  double precision den,t0,t1,t2,lccdf
  integer icond
  t0 = exp(-cpar);  t1 = exp(-cpar*u1); t2 = exp(-cpar*u2);
  den = t1+t2-t0-t1*t2;
  lpdf = log(cpar*(1.d0-t0))-cpar*(u1+u2)-2.d0*log(abs(den)); 
  ccdf = 0.d0;  
  if(icond>=1) then 
    !lccdf = -cpar*u2+log(1.d0-t1)-log(den);  
    ! below for positive or negative cpar
    lccdf = -cpar*u2+log((1.d0-t1)/den);  
    ccdf = exp(lccdf);
  end if
  return
  end

! function with Gumbel lpdf and ccdf for the factor copulas 
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar >1
!   icond = 1 for computing conditional cdf, =0 for not
! outputs 
!   lpdf = log pdf 
!   ccdf = conditional pdf if icond=1
subroutine lgum(u1,u2,cpar,icond,lpdf,ccdf)
  implicit none
  double precision u1,u2,cpar,ccdf,lpdf,lccdf
  double precision x,y,tx,ty,xd,yd,s,m,den,logm
  integer icond
  x = -log(u1); y = -log(u2);  tx = log(x); ty = log(y); 
  xd = x**cpar; yd = y**cpar;  s = xd+yd; m = s**(1.d0/cpar)
  den = m+cpar-1.d0; logm=log(m); 
  lpdf = -m+log(den)+(1.d0-2.d0*cpar)*logm+(cpar-1.d0)*(tx+ty)+x+y;
  ccdf = 0.d0;  
  if(icond >= 1) then
    lccdf = y-m+(1.d0-cpar)*(logm-ty);
    ccdf = exp(lccdf);
  end if  
  return
  end

! function with Student t lpdf and ccdf for the factor copulas 
! inputs
!   t1,t2 = real values
!   rho = scalar in (-1,1)
!   nu = degree of freedom >0
!   icond = 1 for computing conditional cdf, =0 for not
! outputs 
!   lpdf = log pdf 
!   ccdf = conditional pdf if icond=1
subroutine lt(t1,t2,rho,nu,icond,lpdf,ccdf)
  implicit none
  double precision t1,t2,t1sq,t2sq,rho,rho2,nu,ccdf,lpdf
  double precision dt1,dt2,pi,r2,cr,xt,coef,con,den,lgdif
  double precision pt   ! function in C code ptbeta.c
  double precision lgamma   ! function in C code lgamma.c, in gnu fortran?
  integer icond
  pi=3.14159265358979323846d0; ! pi
  lgdif=lgamma(0.5d0*(nu+1.d0))-lgamma(0.5d0*nu);  
  coef = 1.d0-rho*rho;   t1sq = t1*t1; t2sq = t2*t2;
  con=0.5d0*log(pi*nu);
  den = 1.d0+(t1sq-2.d0*rho*t1*t2+t2sq)/(nu*coef);  
  dt1 = lgdif-con-0.5d0*(nu+1.d0)*log(1.d0+t1sq/nu);
  dt2 = lgdif-con-0.5d0*(nu+1.d0)*log(1.d0+t2sq/nu);
  lpdf = -log(2.d0*pi)-0.5d0*log(coef)-0.5d0*(nu+2.d0)*log(den)-dt1-dt2;
  ccdf = 0.d0;  
  if(icond>=1) then
    rho2 = rho*rho;
    r2 = sqrt((nu+t2*t2)/(nu+1));
    cr = sqrt(1.d0-rho2);
    xt = (t1-rho*t2)/cr/r2;
    ccdf = pt(xt,nu+1.d0);
  end if  
  return
  end

! function with BB1 lpdf and ccdf for the factor copulas 
! inputs
!   u1,u2 = values in (0,1)
!   param = (theta,delta), theta>0, delta>1
!   icond = 1 for computing conditional cdf, =0 for not
! outputs 
!   lpdf = log pdf 
!   ccdf = conditional pdf if icond=1
subroutine lbb1(u1,u2,param,icond,lpdf,ccdf)
  implicit none
  double precision u1,u2,param(2),theta,delta,m,lpdf,ccdf
  double precision t1,t2,tu1,tu2,td01,td02,lccdf
  double precision s,t1a,t2a,cf1,dl1,dl1n,mp1,lmp1
  double precision thtem,dltem,coef,smlog
  integer icond
  theta=param(1); delta=param(2);
  t1 = u1**(-theta); t2 = u2**(-theta)
  t1a=t1-1.d0; t2a=t2-1.d0; 
  tu1 = -log(u1); tu2 = -log(u2); 
  td01 = (t1a)**delta; td02 = (t2a)**delta
  dl1n=1.d0-delta; dl1=delta-1.d0;
  thtem=2.d0+1.d0/theta; dltem=1.d0-2.d0*delta;
  coef = theta*delta+1.d0; smlog=log(t1a)+log(t2a);
  s = td01+td02;  m = s**(1.d0/delta); 
  mp1=1.d0+m;   lmp1=log(mp1);
  lpdf = -(thtem)*log(mp1)+(dltem)*log(m)+log(theta*(dl1)+coef*m);
  lpdf = lpdf+(dl1)*(smlog)+(theta+1.d0)*(tu1+tu2);
  ccdf = 0.d0;
  if(icond>=1) then
    cf1 = 1.d0+1.d0/theta;
    lccdf = -cf1*lmp1+(dl1n)*(log(m)-log(t2a))+(theta+1.d0)*tu2;
    ccdf = exp(lccdf);
  end if
  return
  end
