! Gumbel 1-factor and 2-factor copula models
! This file can serve as a template for 1-factor and 2-factor copula
!  models with other bivariate linking copulas to the latent variables.
! Partial derivatives are taken with respect to log bivariate copula density.

! 1-factor model with Gumbel, 
! inputs 
!   npar = #parameters,
!   th = parameter vector (dimension d=npar), th[j]>1
!   d = #variables
!   n = sample size
!   udata = nxd matrix of uniform scores
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes
! outputs 
!   nllk = negative log-likelihood, 
!   grad = gradient of nllk, 
!   hess = hessian of nllk
subroutine gum1fact(npar,th,d,n,udata,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,nq
  double precision th(npar),udata(n,d),wl(nq),xl(nq)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2
  double precision, dimension(:), allocatable :: uvec,lpdf,der1,der2,fval1,integl1,grd
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision fval,integl,ww

  ! npar=d
  allocate ( uvec(d), lpdf(d), der1(d), der2(d), fval1(d), integl1(d), grd(d) )
  allocate ( fval2(d,d), integl2(d,d), hss(d,d) )
  nllk=0.d0; grad=0.d0; hess=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    integl=0.d0; integl1=0.d0; integl2=0.d0; ! initialize for integrals
    do iq=1,nq ! loop over quadrature points
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        call lgum1derivs(xl(iq),uvec(j),th(j),lpdf(j),der1(j),der2(j))
        ! if(iq==1) print "(3f10.6)", lpdf(j),der1(j),der2(j)
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=exp(sum(lpdf))
      ! fval1: vector of partial deriv wrt theta[j], j=1,...,d
      ! fval2: matrix of 2nd partial deriv wrt theta[j], theta[j2], 
      do j=1,d
        fval1(j)=fval*der1(j)
        fval2(j,j)=fval*(der2(j)+der1(j)*der1(j))
        do j2=1,d
          if(j2 /= j) fval2(j,j2)=fval1(j)*der1(j2)
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
        hss(j,j2)=integl2(j,j2)/integl - grd(j)*grd(j2); ! any outer in f90? 
      end do
    end do
    hess=hess-hss
  end do
  deallocate (uvec, lpdf, der1, der2, fval1, fval2, integl1, integl2, grd, hss)
  return
  end

! 2-factor model with Gumbel
! inputs 
!   npar = #parameters,
!   param = parameter vector (dimension 2d=npar)
!   d = #variables
!   n = sample size
!   udata = nxd matrix of uniform scores
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes
! outputs 
!   nllk = negative log-likelihood, 
!   grad = gradient of nllk, 
!   hess = hessian of nllk
subroutine gum2fact(npar,param,d,n,udata,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,nq
  double precision param(npar),udata(n,d),wl(nq),xl(nq)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,iq2
  double precision, dimension(:), allocatable :: uvec,tlpdf,tder1,tder2,fval1,integl1,grd
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision fval,integl,ww,w2
  double precision, dimension(:), allocatable :: th,gam,glpdf,gder1,gder2,gder1u,gder2u,gdermix
  double precision, dimension(:), allocatable :: ccdf,cder1,cder2,cdertem

  ! npar=2*d
  allocate ( uvec(d), tlpdf(d), tder1(d), tder2(d), fval1(npar), integl1(npar), grd(npar) )
  allocate ( th(d), gam(d), glpdf(d), gder1(d), gder2(d), gder1u(d), gder2u(d), gdermix(d) ) 
  allocate ( ccdf(d), cder1(d), cder2(d), cdertem(d) )
  allocate ( fval2(npar,npar), integl2(npar,npar), hss(npar,npar) )
  th(1:d)=param(1:d); gam(1:d)=param((d+1):npar)
  nllk=0.d0; grad=0.d0; hess=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    integl=0.d0; integl1=0.d0; integl2=0.d0; ! initialize for integrals
    !print *, integl1
    do iq=1,nq ! loop over quadrature points
    do iq2=1,nq ! loop over quadrature points
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        ! factor 1 derivatives
        call lgum1derivs(uvec(j),xl(iq),th(j),tlpdf(j),tder1(j),tder2(j))
        ! next lines for testing only
        !call lgum2derivs(uvec(j),xl(iq),th(j),glpdf(j),gder1(j),gder2(j),gder1u(j),gder2u(j),gdermix(j))
        !if(iq==1 .and. iq2==1) print "(7f11.6)", th(j),glpdf(j),gder1(j),gder2(j),gder1u(j),gder2u(j),gdermix(j)
        ! conditionals for factor 1
        call cgumderivs(uvec(j),xl(iq),th(j),ccdf(j),cder1(j),cder2(j))
        !if(iq==1.and.iq2==1) print "(4f11.6)", th(j),ccdf(j),cder1(j),cder2(j)
        ! factor 2 derivatives
        if (ccdf(j) < 1.d-5) ccdf(j) = 1.d-5         !!! to avoid overflow
        if (ccdf(j) > 1.d0-1.d-5) ccdf(j) = 1.d0-1.d-5     !!! to avoid overflow
        call lgum2derivs(ccdf(j),xl(iq2),gam(j),glpdf(j),gder1(j),gder2(j),gder1u(j),gder2u(j),gdermix(j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=exp(sum(tlpdf)+sum(glpdf))
      ! fval1: vector of partial deriv wrt param[j], j=1,...,npar=2*d
      ! fval2: matrix of 2nd partial deriv wrt param[j], param[j2], 
      fval1=0.d0; fval2=0.d0
      do j=1,d
        cdertem(j)= tder1(j) + cder1(j)*gder1u(j) 
        fval1(j+d)=fval*gder1(j)
        fval1(j)=fval*cdertem(j)
        fval2(j+d,j+d)=fval*(gder2(j)+gder1(j)*gder1(j))
        fval2(j,j)=fval*( cdertem(j)*cdertem(j) + tder2(j)   &
           +cder1(j)*cder1(j)*gder2u(j) +cder2(j)*gder1u(j) )
        ! cross th , gam terms for fval2
        fval2(j,j+d)=fval* (cdertem(j)*gder1(j)+ cder1(j)*gdermix(j) ) 
        fval2(j+d,j) = fval2(j,j+d)
      end do
      do j=1,d
        do j2=1,d
          if(j2 /= j) then
            fval2(j+d,j2+d)=fval1(j+d)*gder1(j2)
            fval2(j+d,j2)=fval1(j2)*gder1(j)
            fval2(j,j2)=fval1(j2)*cdertem(j)
            fval2(j,j2+d)=fval1(j)*gder1(j2)
          endif
        end do
      end do
      !if(iq==1 .and. iq2==2) print *, fval
      !if(iq==1 .and. iq2==2) print "(8f10.6,/)", fval1 
      !if(iq==1 .and. iq2==2) print "(8f10.6)", fval2
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
  deallocate (uvec, tlpdf, tder1, tder2, fval1, fval2, integl1, integl2, grd, hss)
  deallocate (th, gam, glpdf, gder1, gder2, gder1u, gder2u, gdermix)
  deallocate (ccdf, cder1, cder2,cdertem)
  return
  end

!--------------------------------------------------------------------
! Remaining code below is translated from R code written by Pavel Krupskii

! Function with Gumbel lpdf, lder1, lder2 (partial wrt cpar and 2nd order)
!  for factor 1.
! For derivatives, lpdf is function of cpar and m=(x^cpar+y^cpar)^(1/cpar)
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (1,oo)
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p cpar
!   lder2 = \p^2 lpdf/\p cpar^2
subroutine lgum1derivs(u1,u2,cpar,lpdf,lder1,lder2)
  implicit none
  double precision u1,u2,cpar,lpdf,lder1,lder2
  double precision x,y,tx,ty,xd,yd,s,m,logs,logm,msq,dlsq,dlcu
  double precision sder1,sder2,mder1,mder2,den,den2
  x = -log(u1); y = -log(u2);
  tx = log(x); ty = log(y); 
  xd = x**cpar; yd = y**cpar; s = xd+yd; m = s**(1.d0/cpar);
  logs=log(s);
  dlsq=cpar*cpar; dlcu=dlsq*cpar
  sder1 = xd*tx+yd*ty;
  sder2 = xd*tx*tx+yd*ty*ty;
  mder1 = m*sder1/(s*cpar)-m*logs/dlsq;
  mder2 = -mder1*logs/dlsq-2.d0*m*sder1/(s*dlsq)+2.d0*m*logs/dlcu;
  mder2 = mder2+sder2*m/(s*cpar)+(mder1/s-m*sder1/(s*s))*sder1/cpar;
  den = m+cpar-1.d0; den2=den*den
  logm=log(m); msq=m*m
  lpdf = -m+log(den)+(1.d0-2.d0*cpar)*logm+(cpar-1.d0)*(tx+ty)+x+y;
  lder1 = -mder1+(mder1+1)/den-2*logm+(1-2*cpar)*mder1/m+tx+ty;
  lder2 = -mder2+mder2/den-(mder1+1.d0)**2/den2-4.d0*mder1/m+(1.d0-2.d0*cpar)*(mder2/m-mder1*mder1/msq);
  return
  end

! Function with Gumbel log pdf c(u1,u2,cpar) for factor 2, and 
!   lder1, lder2 (partial wrt cpar, 1st and 2nd order),
!   also lder1u, lder2u (partial wrt u1, 1st and 2nd order).
!   and ldermix  (partial wrt cpar and u1)
! For derivatives, lpdf is function of cpar and m=(x^cpar+y^cpar)^(1/cpar).
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (1,oo)
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p cpar
!   lder2 = \p^2 lpdf/\p cpar^2
!   lder1u = \p lpdf/\p u1
!   lder2u = \p^2 lpdf/\p u1^2
!   ldermix = \p^2 lpdf/\p u1 \p cpar
subroutine lgum2derivs(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermix)
  implicit none
  double precision u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermix
  double precision x,y,tx,ty,xd,yd,s,m,logs,logm,msq,dlsq,dlcu
  double precision sder1,sder2,mder1,mder2,den,den2
  double precision mu,m2u,u1sq,muder1
  x = -log(u1); y = -log(u2);
  tx = log(x); ty = log(y); 
  xd = x**cpar; yd = y**cpar; s = xd+yd; m = s**(1.d0/cpar);
  logs=log(s);
  dlsq=cpar*cpar; dlcu=dlsq*cpar
  ! for 1-factor and 2-factor models
  sder1 = xd*tx+yd*ty;
  sder2 = xd*tx*tx+yd*ty*ty;
  mder1 = m*sder1/(s*cpar)-m*logs/dlsq;
  mder2 = -mder1*logs/dlsq-2.d0*m*sder1/(s*dlsq)+2.d0*m*logs/dlcu;
  mder2 = mder2+sder2*m/(s*cpar)+(mder1/s-m*sder1/(s*s))*sder1/cpar;
  den = m+cpar-1.d0; den2=den*den
  logm=log(m); msq=m*m
  lpdf = -m+log(den)+(1.d0-2.d0*cpar)*logm+(cpar-1.d0)*(tx+ty)+x+y;
  lder1 = -mder1+(mder1+1.d0)/den-2.d0*logm+(1.d0-2.d0*cpar)*mder1/m+tx+ty;
  lder2 = -mder2+mder2/den-(mder1+1.d0)**2/den2-4.d0*mder1/m+(1.d0-2.d0*cpar)*(mder2/m-mder1*mder1/msq);
  ! for 2-factor model
  u1sq=u1*u1;
  mu = -m*xd/(u1*s*x); 
  m2u = (1.d0-cpar)*m*xd*xd/(u1*s*x)**2+(cpar-1.d0)*m*xd/(s*x*x*u1sq)+m*xd/(s*x*u1sq);
  muder1 = -(mder1/s-m*sder1/(s*s))*xd/(x*u1)-m*xd*tx/(s*u1*x);
  lder1u = -mu+mu/den+(1.d0-2.d0*cpar)*mu/m-(cpar-1.d0)/(u1*x)-1.d0/u1;
  lder2u = -m2u+m2u/den-mu*mu/den2+(1.d0-2.d0*cpar)*(m2u/m-mu*mu/msq)+(cpar-1.d0)/(u1sq*x); 
  lder2u = lder2u-(cpar-1.d0)/(x*x*u1sq)+1.d0/u1sq;
  ldermix = -muder1+muder1/den-mu*(mder1+1.d0)/den2-2.d0*mu/m+(1.d0-2.d0*cpar)*(muder1/m-mu*mder1/msq)-1.d0/(x*u1);
  return
  end

! Function with Gumbel condcdf C_{1|2}(u1|u2;cpar) for factor 1, and 
! ccdf cder1, cder2 (partial wrt cpar, 1st and 2nd order)
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (1,oo)
! outputs
!   ccdf = conditional cdf, 
!   cder1 = \p ccdf/\p cpar
!   cder2 = \p^2 ccdf/\p cpar^2
subroutine cgumderivs(u1,u2,cpar,ccdf,cder1,cder2)
  implicit none
  double precision u1,u2,cpar,ccdf,cder1,cder2
  double precision x,y,tx,ty,xd,yd,s,m,logs,logm,msq,dlsq,dlcu
  double precision sder1,sder2,mder1,mder2
  double precision lccdf,lcder1,lcder2
  x = -log(u1); y = -log(u2);
  tx = log(x); ty = log(y); 
  xd = x**cpar; yd = y**cpar; s = xd+yd; m = s**(1.d0/cpar);
  logs=log(s);
  dlsq=cpar*cpar; dlcu=dlsq*cpar
  sder1 = xd*tx+yd*ty;
  sder2 = xd*tx*tx+yd*ty*ty;
  mder1 = m*sder1/(s*cpar)-m*logs/dlsq;
  mder2 = -mder1*logs/dlsq-2.d0*m*sder1/(s*dlsq)+2.d0*m*logs/dlcu;
  mder2 = mder2+sder2*m/(s*cpar)+(mder1/s-m*sder1/(s*s))*sder1/cpar;
  logm=log(m); msq=m*m
  lccdf = y-m+(1.d0-cpar)*(logm-ty)
  !lcder1 = -mder1+(1.d0-cpar)*mder1/m-logm+tx;
  lcder1 = -mder1+(1.d0-cpar)*mder1/m-logm+ty;  ! given y
  lcder2 = -mder2-2.d0*mder1/m+(1.d0-cpar)*(mder2/m-mder1*mder1/msq);
  ccdf = exp(lccdf)
  cder1 = ccdf*lcder1
  cder2 = ccdf*(lcder1*lcder1+lcder2)
  return
  end

