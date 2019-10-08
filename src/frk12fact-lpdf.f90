! Frank 1-factor and 2-factor copula models
! This code follows the style of gum12fact.f90
! Code written by Pavel Krupskii

! 1-factor copula with Frank
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
subroutine lfrk1fact(npar,th,d,n,udata,nq,wl,xl,nllk,grad,hess)
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
        call lfrk1derivs(xl(iq),uvec(j),th(j),lpdf(j),der1(j),der2(j))
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

! 2-factor copula with Frank
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
subroutine lfrk2fact(npar,param,d,n,udata,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,nq
  double precision param(npar),udata(n,d),wl(nq),xl(nq)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,iq2
  double precision, dimension(:), allocatable :: uvec,tlpdf,tder1,tder2,fval1,integl1,grd
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision fval,integl,ww,w2,epsvl
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
        call lfrk1derivs(uvec(j),xl(iq),th(j),tlpdf(j),tder1(j),tder2(j))
        ! next lines for testing only
        !call lfrk2derivs(uvec(j),xl(iq),th(j),glpdf(j),gder1(j),gder2(j),gder1u(j),gder2u(j),gdermix(j))
        !if(iq==1 .and. iq2==1) print "(7f11.6)", th(j),glpdf(j),gder1(j),gder2(j),gder1u(j),gder2u(j),gdermix(j)
        ! conditionals for factor 1
        call clfrkderivs(uvec(j),xl(iq),th(j),ccdf(j),cder1(j),cder2(j))
        !if(iq==1.and.iq2==1) print "(4f11.6)", th(j),ccdf(j),cder1(j),cder2(j)
        ! factor 2 derivatives
        epsvl = 0.00000001d0                    !to prevent the first argument of copula to be close to zero or one
        if(ccdf(j)<epsvl)   ccdf(j)=epsvl     !helps to avoid "singular hessian" error for Gumbel copula
        if(ccdf(j)>1.d0-epsvl) ccdf(j)=1.d0-epsvl   !

        call lfrk2derivs(ccdf(j),xl(iq2),gam(j),glpdf(j),gder1(j),gder2(j),gder1u(j),gder2u(j),gdermix(j))
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

! Function with Frank lpdf, lder1, lder2 (partial wrt cpar and 2nd order)
!  for factor 1.
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (-oo,oo) but could be underflow problems for cpar>35
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p cpar
!   lder2 = \p^2 lpdf/\p cpar^2
subroutine lfrk1derivs(u1,u2,cpar,lpdf,lder1,lder2)
  implicit none
  double precision u1,u2,cpar,lpdf,lder1,lder2
  double precision den,den1,den2,t0,t1,t2
  t0 = exp(-cpar);
  t1 = exp(-cpar*u1);
  t2 = exp(-cpar*u2);
  den = t1+t2-t0-t1*t2;
  den1 = -u1*t1-u2*t2+t0+(u1+u2)*t1*t2;
  den2 = u1*u1*t1+u2*u2*t2-t0-(u1+u2)*(u1+u2)*t1*t2;
  lpdf = log(abs(cpar))+log(abs(1.d0-t0))-cpar*(u1+u2)-2.d0*log(abs(den));
  lder1 = 1.d0/cpar+t0/(1.d0-t0)-(u1+u2)-2.d0*den1/den;
  lder2 = -1.d0/(cpar*cpar)-t0/((1.d0-t0)*(1.d0-t0))-2.d0*den2/den+2.d0*den1*den1/(den*den); 
  return
  end

! Function with Frank log pdf c(u1,u2,cpar) for factor 2, and 
!   lder1, lder2 (partial wrt cpar, 1st and 2nd order)
!   also lder1u lder2u (partial wrt u1, 1st and 2nd order)
!   and ldermix  (partial wrt cpar and u1) 
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (-oo,oo)
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p cpar
!   lder2 = \p^2 lpdf/\p cpar^2
!   lder1u = \p lpdf/\p u1
!   lder2u = \p^2 lpdf/\p u1^2
!   ldermix = \p^2 lpdf/\p u1 \p cpar
subroutine lfrk2derivs(u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermix)
  implicit none
  double precision u1,u2,cpar,lpdf,lder1,lder2,lder1u,lder2u,ldermix
  double precision den,den1,den2,den1u,den2u,denmix,t0,t1,t2

  t0 = exp(-cpar);
  t1 = exp(-cpar*u1);
  t2 = exp(-cpar*u2);
  den = t1+t2-t0-t1*t2;
  den1 = -u1*t1-u2*t2+t0+(u1+u2)*t1*t2;
  den2 = u1*u1*t1+u2*u2*t2-t0-(u1+u2)*(u1+u2)*t1*t2;
  den1u = -cpar*t1*(1.d0-t2);  
  den2u = cpar*cpar*t1*(1.d0-t2);
  denmix = t1*(-1.d0+cpar*u1+t2-(u1+u2)*cpar*t2); 
   
  !lpdf = log(abs(cpar))+log(abs(1-t0))-cpar*(u1+u2)-2*log(abs(den));
  ! 121023 maybe later add the limits as cpar->0
  ! pdf = cpar*(1-t0)/den^2 where
  !    1-t0 has same sign as cpar,  den has same sign as cpar
  lpdf = log(cpar*(1.d0-t0))-cpar*(u1+u2)-2.d0*log(abs(den));
  lder1 = 1.d0/cpar+t0/(1.d0-t0)-(u1+u2)-2.d0*den1/den;
  lder2 = -1.d0/(cpar*cpar)-t0/((1.d0-t0)*(1.d0-t0))-2.d0*den2/den+2.d0*den1*den1/(den*den); 
  lder1u = -cpar-2.d0*den1u/den;
  lder2u = -2.d0*den2u/den+2.d0*den1u*den1u/(den*den);
  ldermix = -1.d0-2.d0*denmix/den+2.d0*den1u*den1/(den*den);
  return
  end

! Function with Frank condcdf C_{1|2}(u1|u2;cpar) for factor 1, and
! ccdf cder1, cder2 (partial wrt cpar, 1st and 2nd order)
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (-oo,oo)
! outputs
!   ccdf = conditional cdf, 
!   cder1 = \p ccdf/\p cpar
!   cder2 = \p^2 ccdf/\p cpar^2
subroutine clfrkderivs(u1,u2,cpar,ccdf,cder1,cder2)
  implicit none
  double precision u1,u2,cpar,ccdf,cder1,cder2
  double precision t0,t1,t2,den,den1,den2,lccdf,lcder1,lcder2
  
  t0 = exp(-cpar);
  t1 = exp(-cpar*u1);
  t2 = exp(-cpar*u2);
  den = t1+t2-t0-t1*t2;
  den1 = -u1*t1-u2*t2+t0+(u1+u2)*t1*t2;
  den2 = u1*u1*t1+u2*u2*t2-t0-(u1+u2)*(u1+u2)*t1*t2;
  !lccdf = -cpar*u2+log(1-t1)-log(den);
  lccdf = -cpar*u2+log((1.d0-t1)/den);
  lcder1 = -u2+u1*t1/(1.d0-t1)-den1/den;
  lcder2 = -u1*u1*t1/((1.d0-t1)*(1.d0-t1))-den2/den+den1*den1/(den*den);
  ccdf = exp(lccdf)
  cder1 = ccdf*lcder1
  cder2 = ccdf*(lcder1*lcder1+lcder2)
  return
  end

