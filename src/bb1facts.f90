! BB1 1-factor and BB1/Frank 2-factor copula models
! This file can serve as a template for 1-factor and 2-factor copula
!  when the linking copula for the first factor has 2 parameters and
!  the linking copula for the second factor has 1 parameter.

! sample main program
!program bb1
!  implicit none
!  integer n,d,nq,i,iq,j,npar
!  double precision, dimension(:,:), allocatable :: udata,hess
!  double precision, dimension(:), allocatable :: xl,wl,grad,param
!  double precision nllk
!  read *,nq
!  allocate ( xl(nq),wl(nq) )
!  do iq =1,nq   ! equidistant and equially weighted for testing
!    xl(iq)=iq/(nq+1.d0)
!    wl(iq)=1.d0/nq
!  end do
!  read *,n,d
!  npar=d*3
!  allocate ( udata(n,d), param(npar), grad(npar), hess(npar,npar) )
!  do i=1,n
!    read *, udata(i,:)
!  end do
!  do j=1,d
!    param(j)=.1d0*j
!    param(j+d)=1.d0+.1d0*j
!    param(j+d*2)=0.5d0+.1d0*j
!  end do
!  call bb1frk2fact(npar,param,d,n,udata,nq,wl,xl,nllk,grad,hess)
!  print *, nllk
!  print "(8f10.5,/)", grad
!  print "(8f10.5)", hess
!  deallocate (xl, wl, param, udata, grad, hess)
!  stop
!  end


! 1-factor model with BB1 
! inputs 
!   npar = #parameters,
!   param0 = parameter vector (dimension d=npar), 
!       param0 has th1,dl1,th2,dl2, ... as vector theta>0, delta>1
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
subroutine bb11fact(npar,param0,d,n,udata,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,nq
  double precision param0(npar),udata(n,d),wl(nq),xl(nq)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,jj(2),jj2(2)
  double precision, dimension(:), allocatable :: uvec,lpdf,fval1,integl1,grd
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss,der1,der2
  double precision, dimension(:,:), allocatable :: param
  double precision fval,integl,ww,omat(2,2)

  ! npar=d*2
  allocate ( uvec(d), lpdf(d), der1(2,d), der2(3,d), fval1(npar), integl1(npar), grd(npar) )
  allocate ( fval2(npar,npar), integl2(npar,npar), hss(npar,npar), param(2,d) )
  param=reshape(param0,(/2,d/))
  nllk=0.d0; grad=0.d0; hess=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    integl=0.d0; integl1=0.d0; integl2=0.d0; ! initialize for integrals
    do iq=1,nq ! loop over quadrature points
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        call lbb1derivs(xl(iq),uvec(j),param(:,j),lpdf(j),der1(:,j),der2(:,j))
        !if(iq==1) print "(6f11.6)", lpdf(j),der1(:,j),der2(:,j)
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=exp(sum(lpdf))
      ! fval1: vector of partial deriv wrt param[j], j=1,...,np
      ! fval2: matrix of 2nd partial deriv wrt param[j], param[j2], 
      do j=1,d
        jj=(/2*j-1,2*j/)
        fval1(jj)=fval*der1(:,j)
        call outer(2,2,der1(:,j),der1(:,j),omat)
        omat(1,1)=omat(1,1)+der2(1,j)
        omat(1,2)=omat(1,2)+der2(2,j)
        omat(2,1)=omat(2,1)+der2(2,j)
        omat(2,2)=omat(2,2)+der2(3,j)
        fval2(jj,jj)=fval*omat
        do j2=1,d
          if(j2 /= j) then
            call outer(2,2,der1(:,j),der1(:,j2),omat)
            jj2=(/2*j2-1,2*j2/)
            fval2(jj,jj2)=fval*omat
          endif
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
    do j=1,npar
      do j2=1,npar
        hss(j,j2)=integl2(j,j2)/integl - grd(j)*grd(j2); 
      end do
    end do
    hess=hess-hss
  end do
  deallocate (uvec, lpdf, der1, der2, fval1, fval2, integl1, integl2, grd, hss)
  deallocate ( param )
  return
  end

! 2-factor model with BB1/Frank (linking copula family for latent1/latent2) 
! The subroutine lfrk2derivs is in file frk12fact-lpdf.f90
! inputs 
!   npar = #parameters,
!   param0 = parameter vector (dimension d=npar), 
!         BB1 parameters th1,dl1,th2,dl2, ...  theta>0, delta>1
!         Frank parameters gam1,gam2, ...  -oo<gamma<oo
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
subroutine bb1frk2fact(npar,param0,d,n,udata,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,nq
  double precision param0(npar),udata(n,d),wl(nq),xl(nq)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,d2,iq2,jj(2),jj2(2)
  double precision, dimension(:), allocatable :: uvec,tlpdf,fval1,integl1,grd
  double precision, dimension(:), allocatable :: gam,glpdf,gder1,gder2,gder1u,gder2u,gdermix
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss,tder1,tder2
  double precision, dimension(:,:), allocatable :: param,cder1,cder2,cdertem
  double precision fval,integl,ww,w2,omat(2,2),omatc(2,2)
  double precision, dimension(:), allocatable :: ccdf

  !npar=d*3
  allocate ( uvec(d), tlpdf(d), tder1(2,d), tder2(3,d), fval1(npar), integl1(npar), grd(npar) )
  allocate ( gam(d), glpdf(d), gder1(d), gder2(d), gder1u(d), gder2u(d), gdermix(d) ) 
  allocate ( ccdf(d), cder1(2,d), cder2(3,d), cdertem(2,d) )
  allocate ( fval2(npar,npar), integl2(npar,npar), hss(npar,npar), param(2,d) )
  d2=2*d
  param=reshape(param0(1:d2),(/2,d/))
  gam(1:d)=param0((d2+1):npar)
  nllk=0.d0; grad=0.d0; hess=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    integl=0.d0; integl1=0.d0; integl2=0.d0; ! initialize for integrals
    do iq=1,nq ! loop over quadrature points
    do iq2=1,nq ! loop over quadrature points
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        call lbb1derivs(xl(iq),uvec(j),param(:,j),tlpdf(j),tder1(:,j),tder2(:,j))
        !if(iq==1.and.iq2==1) print "(6f11.6)", tlpdf(j),tder1(:,j),tder2(:,j)
        call cbb1derivs(uvec(j),xl(iq),param(:,j),ccdf(j),cder1(:,j),cder2(:,j))
        !if(iq==1.and.iq2==1) print "(8f11.6)", param(:,j),ccdf(j),cder1(:,j),cder2(:,j)
        ! factor 2 derivatives
        if (ccdf(j) < 1.d-5) ccdf(j) = 1.d-5         !!! to avoid overflow
        if (ccdf(j) > 1.d0-1.d-5) ccdf(j) = 1.d0-1.d-5     !!! to avoid overflow
        call lfrk2derivs(ccdf(j),xl(iq2),gam(j),glpdf(j),gder1(j),gder2(j),gder1u(j),gder2u(j),gdermix(j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=exp(sum(tlpdf)+sum(glpdf))
      ! fval1: vector of partial deriv wrt param[j], j=1,...,npar=2*d
      ! fval2: matrix of 2nd partial deriv wrt param[j], param[j2], 
      fval1=0.d0; fval2=0.d0
      do j=1,d
        cdertem(:,j)= tder1(:,j) + cder1(:,j)*gder1u(j) 
        fval1(j+d2)=fval*gder1(j)
        jj=(/2*j-1,2*j/)
        !fval1(j)=fval*cdertem(j)
        fval1(jj)=fval*cdertem(:,j)
        fval2(j+d2,j+d2)=fval*(gder2(j)+gder1(j)*gder1(j))
        !fval2(j,j+d)=fval* (cdertem(j)*gder1(j)+ cder1(j)*gdermix(j) ) 
        fval2(jj,j+d2)=fval* (cdertem(:,j)*gder1(j)+ cder1(:,j)*gdermix(j) ) 
        !fval2(j+d2,j) = fval2(j,j+d)
        fval2(j+d2,2*j-1) = fval2(2*j-1,j+d2)
        fval2(j+d2,2*j) = fval2(2*j,j+d2)
        !! cross th , gam terms for fval2
        !fval2(j,j)=fval*( cdertem(j)*cdertem(j) + tder2(j)   &
        !   +cder1(j)*cder1(j)*gder2u(j) +cder2(j)*gder1u(j) )
        !call outer(2,2,tder1(:,j),tder1(:,j),omat)
        call outer(2,2,cdertem(:,j),cdertem(:,j),omat)
        call outer(2,2,cder1(:,j),cder1(:,j),omatc)
        omat(1,1)=omat(1,1)+tder2(1,j) +cder2(1,j)*gder1u(j) +omatc(1,1)*gder2u(j)
        omat(1,2)=omat(1,2)+tder2(2,j) +cder2(2,j)*gder1u(j) +omatc(1,2)*gder2u(j)
        omat(2,1)=omat(2,1)+tder2(2,j) +cder2(2,j)*gder1u(j) +omatc(2,1)*gder2u(j)
        omat(2,2)=omat(2,2)+tder2(3,j) +cder2(3,j)*gder1u(j) +omatc(2,2)*gder2u(j)
      end do
      do j=1,d
        jj=(/2*j-1,2*j/)
        do j2=1,d
          if(j2 /= j) then
            jj2=(/2*j2-1,2*j2/)
            !fval2(j+d,j2+d)=fval1(j+d)*gder1(j2)
            !fval2(j+d2,j2+d2)=fval1(j+d)*gder1(j2)  ! error
            fval2(j+d2,j2+d2)=fval1(j+d2)*gder1(j2)
            !fval2(j+d,j2)=fval1(j2)*gder1(j)
            fval2(j+d2,jj2)=fval1(jj2)*gder1(j)
            !fval2(j,j2+d)=fval1(j)*gder1(j2)
            fval2(jj,j2+d2)=fval1(jj)*gder1(j2)
            !fval2(j,j2)=fval1(j2)*cdertem(j)
            call outer(2,2,cdertem(:,j),cdertem(:,j2),omat)
            fval2(jj,jj2)=fval*omat
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
  deallocate (param, gam, glpdf, gder1, gder2, gder1u, gder2u, gdermix)
  deallocate (ccdf, cder1, cder2,cdertem)
  return
  end

!--------------------------------------------------------------------
! Remaining code below is translated from R code written by Pavel Krupskii

! Function with BB1 lpdf, lder1, lder2 (partial wrt param and 2nd order)
!  for factor 1;
! lder1 is 2-vector, lder2 is 2x2 matrix
! inputs
!   u1,u2 = values in (0,1)
!   param = (theta, delta), theta>0, delta>1
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p param
!   lder2 = \p^2 lpdf/\p param \p param^T
subroutine lbb1derivs(u1,u2,param,lpdf,lder1,lder2)
  implicit none
  double precision u1,u2,param(2),lpdf,lder1(2),lder2(3)
  double precision theta,delta,der1th,der1dl,der2th,der2dl,derthd
  double precision t10,t1der1th,t1der1dl,t1der2th,t1der2dl,t1derthd
  double precision t20,t2der1th,t2der1dl,t2der2th,t2der2dl,t2derthd
  double precision t30,t3der1th,t3der1dl,t3der2th,t3der2dl,t3derthd
  double precision t40,t4der1th,t4der1dl,t4der2th,t4der2dl,t4derthd
  double precision m,mder1th,mder1dl,mder2th,mder2dl,mderthd
  double precision den,coef,t1,t2,tu1,tu2
  double precision mp1,msq,thsq,den2,thtem,dltem,dl1,t1a,t2a,smlog

  theta=param(1); delta=param(2);
  !!if(theta < 1e-5) theta = 1e-5;     !!to avoid overflow
  !!if(delta < 1+1e-5) delta = 1+1e-5; !!to avoid overflow
  !!didn't work: too flat likelihood 

  call mderivs(u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd)
  mp1=1.d0+m; msq=m*m; thsq=theta*theta
  thtem=2.d0+1.d0/theta
  t10 = -(thtem)*log(mp1)
  t1der1th = log(mp1)/thsq - (thtem)*mder1th/mp1
  t1der1dl = -(thtem)*mder1dl/mp1
  t1der2th = -2.d0*log(mp1)/theta**3+(2.d0/thsq)*mder1th/mp1
  t1der2th = t1der2th-(thtem)*(mder2th/mp1-mder1th*mder1th/(mp1*mp1))
  t1der2dl = -(thtem)*(mder2dl/mp1-mder1dl*mder1dl/(mp1*mp1))
  t1derthd = mder1dl/(thsq*(mp1))-(thtem)*(mderthd/(mp1)-mder1th*mder1dl/(mp1*mp1))

  dltem=1.d0-2.d0*delta; dl1=delta-1.d0
  t20 = (dltem)*log(m)
  t2der1th = (dltem)*mder1th/m
  t2der1dl = -2.d0*log(m)+(dltem)*mder1dl/m
  t2der2th = (dltem)*(mder2th/m-mder1th*mder1th/msq)
  t2der2dl = -4.d0*mder1dl/m+(dltem)*(mder2dl/m-mder1dl*mder1dl/msq)
  t2derthd = -2.d0*mder1th/m+(dltem)*(mderthd/m-mder1th*mder1dl/msq)

  coef = theta*delta+1.d0
  den = theta*(dl1)+(coef)*m
  den2=den*den
  t30 = log(theta*(dl1)+coef*m)
  t3der1th = (dl1+delta*m+coef*mder1th)/den
  t3der1dl = (theta+theta*m+coef*mder1dl)/den
  t3der2th = (2.d0*delta*mder1th+coef*mder2th)/den
  t3der2th = t3der2th-(dl1+delta*m+coef*mder1th)**2/den2
  t3der2dl = (2.d0*theta*mder1dl+coef*mder2dl)/den
  t3der2dl = t3der2dl-(theta+theta*m+coef*mder1dl)**2/den2;
  t3derthd = (1.d0+m+theta*mder1th+delta*mder1dl+coef*mderthd)/den;
  t3derthd = t3derthd - (theta+theta*m+coef*mder1dl)*(dl1+delta*m+coef*mder1th)/den2

  t1 = u1**(-theta); t2 = u2**(-theta)
  t1a=t1-1.d0; t2a=t2-1.d0; smlog=log(t1a)+log(t2a)
  tu1 = -log(u1); tu2 = -log(u2)
  t40 = (dl1)*(smlog)+(theta+1.d0)*(tu1+tu2)
  t4der1th = (dl1)*(t1*tu1/t1a+t2*tu2/t2a)+tu1+tu2  
  t4der1dl = smlog
  t4der2th = -(dl1)*(t1*tu1*tu1/t1a**2+t2*tu2*tu2/t2a**2)
  t4der2dl = 0.d0
  t4derthd = t1*tu1/t1a+t2*tu2/t2a

  lpdf = t10+t20+t30+t40
  der1th = t1der1th+t2der1th+t3der1th+t4der1th
  der1dl = t1der1dl+t2der1dl+t3der1dl+t4der1dl
  der2th = t1der2th+t2der2th+t3der2th+t4der2th
  der2dl = t1der2dl+t2der2dl+t3der2dl+t4der2dl
  derthd = t1derthd+t2derthd+t3derthd+t4derthd
  lder1(1)=der1th; lder1(2)=der1dl  ! gradient of log pdf
  lder2(1)=der2th; lder2(2)=derthd; lder2(3)=der2dl; ! Hessian terms
  return
  end

! Function with BB1 condcdf C_{1|2}(u1|u2;theta,delta) for factor 1, and
! ccdf cder1(2), cder2(3) (partial wrt param, 1st and 2nd order),
! second deriv has theta^2, mix, delta^2
! inputs
!   u1,u2 = values in (0,1)
!   param = (theta, delta), theta>0, delta>1
! outputs
!   ccdf = conditional cdf, 
!   cder1 = \p ccdf/\p param
!   cder2 = \p^2 ccdf/\p param \p param^T
subroutine cbb1derivs(u1,u2,param,ccdf,cder1,cder2)
  implicit none
  double precision u1,u2,param(2),ccdf,cder1(2),cder2(3)
  double precision theta,delta
  double precision m,mder1th,mder1dl,mder2th,mder2dl,mderthd
  double precision t2,tu2,cf1,logm
  double precision mp1,msq,thsq,dl1n,t2a,lt2a,lmp1
  double precision lcdf,lder1th,lder1dl,lder2th,lder2dl,lderthd
  double precision cder1th,cder1dl,cder2th,cder2dl,cderthd

  theta=param(1); delta=param(2)
  !!if(theta < 1e-5) theta = 1e-5;     !!to avoid overflow
  !!if(delta < 1+1e-5) delta = 1+1e-5; !!to avoid overflow
  !!didn't work: too flat likelihood

  call mderivs(u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd)
  mp1=1.d0+m; msq=m*m; thsq=theta*theta
  lmp1=log(mp1); logm=log(m);

  t2 = u2**(-theta); tu2 = -log(u2);
  t2a=t2-1.d0; lt2a=log(t2a)
  cf1 = 1.d0+1.d0/theta;
  dl1n=1.d0-delta

  lcdf = -cf1*lmp1+(dl1n)*(logm-lt2a)+(theta+1.d0)*tu2;
  lder1th = lmp1/thsq-cf1*mder1th/mp1;
  lder1th = lder1th+(dl1n)*(mder1th/m-t2*tu2/t2a)+tu2;
  lder1dl = -cf1*mder1dl/mp1+(dl1n)*mder1dl/m-logm+lt2a;
  lder2th = -2.d0*lmp1/(thsq*theta)+(2.d0/thsq)*mder1th/mp1-cf1*(mder2th/mp1-mder1th**2/(mp1*mp1));
  lder2th = lder2th+(dl1n)*(mder2th/m-mder1th**2/msq+tu2*tu2*t2/(t2a*t2a));
  lder2dl = -cf1*(mder2dl/mp1-mder1dl**2/(mp1*mp1))-2.d0*mder1dl/m;
  lder2dl = lder2dl+(dl1n)*(mder2dl/m-mder1dl**2/msq);
  lderthd = (1.d0/thsq)*mder1dl/mp1-cf1*(mderthd/mp1-mder1th*mder1dl/(mp1*mp1))-mder1th/m;
  lderthd = lderthd+(dl1n)*(mderthd/m-mder1th*mder1dl/msq)+t2*tu2/t2a;

  ccdf = exp(lcdf);
  cder1th=ccdf*lder1th;
  cder1dl=ccdf*lder1dl;
  cder2dl=ccdf*(lder2dl+lder1dl**2);
  cder2th=ccdf*(lder2th+lder1th**2);
  cderthd=ccdf*(lderthd+lder1dl*lder1th);
  cder1(1)=cder1th; cder1(2)=cder1dl  ! gradient of log pdf
  cder2(1)=cder2th; cder2(2)=cderthd; cder2(3)=cder2dl; ! Hessian terms
  return
  end

! 1st and 2nd order derivatives of m = s^(1/delta) wrt theta and delta
! inputs
!   u1,u2 = values in (0,1)
!   theta >0
!   delta >1
! outputs
!   m = s^(1/delta); s =  [u1^(-theta)-1]^delta + [u2^(-theta)-1]^delta
!   mder1th = \p m/\p theta
!   mder1dl = \p m/\p delta
!   mder2th = \p^2 m/\p theta^2
!   mder2dl = \p^2 m/\p delta^2
!   mderthd = \p^2 m/\p theta \p delta
subroutine mderivs(u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd)
  implicit none
  double precision u1,u2,theta,delta,m,mder1th,mder1dl,mder2th,mder2dl,mderthd
  double precision t1,t2,tu1,tu2,ttu1,ttu2,td01,td02,td11,td12,td21,td22
  double precision s,sder1th,sder1dl,sder2th,sder2dl,sderthd,m1,ts,m1der1dl 
  double precision t1a,t2a,dlsq,dlcu

  t1 = u1**(-theta); t2 = u2**(-theta)
  t1a=t1-1.d0; t2a=t2-1.d0
  tu1 = -log(u1); tu2 = -log(u2)
  ttu1 = log(t1a); ttu2 = log(t2a)
  td01 = (t1a)**delta; td02 = (t2a)**delta
  !td11 = (t1-1)**(delta-1); td12 = (t2-1)**(delta-1)
  !td21 = (t1-1)**(delta-2); td22 = (t2-1)**(delta-2)
  td11=td01/t1a; td12=td02/t2a;
  td21=td11/t1a; td22=td12/t2a;
 
  s = td01+td02
  sder1th = delta*(td11*t1*tu1+td12*t2*tu2)
  sder1dl = td01*ttu1+td02*ttu2
  sder2th = delta*(delta-1.d0)*(td21*t1*t1*tu1*tu1+td22*t2*t2*tu2*tu2)
  sder2th = sder2th+delta*(td11*t1*tu1*tu1+td12*t2*tu2*tu2)
  sder2dl = td01*ttu1*ttu1+td02*ttu2*ttu2
  sderthd = sder1th/delta+delta*(td11*ttu1*tu1*t1+td12*ttu2*tu2*t2)

  m = s**(1.d0/delta); m1 = m/s
  ts = log(s)
  dlsq=delta*delta; dlcu=delta*dlsq
  mder1th = m1*sder1th/delta
  mder1dl = m1*sder1dl/delta - m*ts/dlsq
  m1der1dl = mder1dl/s - m*sder1dl/s**2
  mder2th = (1.d0-delta)*m1*sder1th**2/(dlsq*s)+m1*sder2th/delta
  mder2dl = 2.d0*m*ts/dlcu-mder1dl*ts/dlsq-2.d0*m1*sder1dl/dlsq
  mder2dl = mder2dl+sder2dl*m1/delta+sder1dl*m1der1dl/delta
  mderthd = -m1*sder1th/dlsq+sder1th*m1der1dl/delta+m1*sderthd/delta
  return
  end

