! Frank 1-factor and 2-factor copula models
! This file is a template for 1-factor and 2-factor copula when
! partial derivatives are taken with respect to the bivariate copula density
! It is different from the other *12fact*.f90 codes.

! Code written by Pavel Krupskii

! gfortran -fpic -c frk12fact.f90
!  out= .Fortran("frk2fact",
!  out= .Fortran("frk1fact",
!      as.integer(npar), as.double(th), as.integer(d), as.integer(n),
!      as.double(udata), 
!      as.integer(nq), as.double(wl), as.double(xl), 
!      nllk=as.double(0.),grad=as.double(rep(0,npar)),
!      hess=as.double(rep(0,npar*npar))  )

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
subroutine frk1fact(npar,th,d,n,udata,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,nq
  double precision th(npar),udata(n,d),wl(nq),xl(nq)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2
  double precision, dimension(:), allocatable :: uvec,pdfv,der1,der2,fval1,integl1,grd
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision fval,integl,ww

  ! npar=d
  allocate ( uvec(d), pdfv(d), der1(d), der2(d), fval1(d), integl1(d), grd(d) )
  allocate ( fval2(d,d), integl2(d,d), hss(d,d) )
  nllk=0.d0; grad=0.d0; hess=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    !print "(i4,10f7.3)", i, uvec
    integl=0.d0; integl1=0.d0; integl2=0.d0; ! initialize for integrals
    do iq=1,nq ! loop over quadrature points
      ! get \p copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        call frk1derivs(xl(iq),uvec(j),th(j),pdfv(j),der1(j),der2(j))
        ! print "(3f10.6)", pdfv(j),der1(j),der2(j)
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=product(pdfv)
      ! fval1: vector of partial deriv wrt theta[j], j=1,...,d
      ! fval2: matrix of 2nd partial deriv wrt theta[j], theta[j2], 
      do j=1,d
        fval1(j)=fval*der1(j)/pdfv(j)
        fval2(j,j)=fval*der2(j)/pdfv(j)
        do j2=1,d
          !if(j2.ne.j) fval2(j,j2)=fval1(j)*der1(j2)/pdfv(j2)
          if(j2 /= j) fval2(j,j2)=fval1(j)*der1(j2)/pdfv(j2)
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
  deallocate (uvec, pdfv, der1, der2, fval1, fval2, integl1, integl2, grd, hss)
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
subroutine frk2fact(npar,param,d,n,udata,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,nq
  double precision param(npar),udata(n,d),wl(nq),xl(nq)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,iq2
  double precision, dimension(:), allocatable :: uvec,tpdf,tder1,tder2,fval1,integl1,grd
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision fval,integl,ww,w2
  double precision, dimension(:), allocatable :: th,gam,gpdf,gder1,gder2,gder1u,gder2u,gdermix
  double precision, dimension(:), allocatable :: ccdf,cder1,cder2

  ! npar=2*d
  allocate ( uvec(d), tpdf(d), tder1(d), tder2(d), fval1(npar), integl1(npar), grd(npar) )
  allocate ( th(d), gam(d), gpdf(d), gder1(d), gder2(d), gder1u(d), gder2u(d), gdermix(d) ) 
  allocate ( ccdf(d), cder1(d), cder2(d) )
  allocate ( fval2(npar,npar), integl2(npar,npar), hss(npar,npar) )
  th(1:d)=param(1:d); gam(1:d)=param((d+1):npar)
  nllk=0.d0; grad=0.d0; hess=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    integl=0.d0; integl1=0.d0; integl2=0.d0; ! initialize for integrals
    do iq=1,nq ! loop over quadrature points
    do iq2=1,nq ! loop over quadrature points
      ! get \p copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        ! factor 1 derivatives
        call frk1derivs(uvec(j),xl(iq),th(j),tpdf(j),tder1(j),tder2(j))
        ! print "(3f10.6)", pdfv(j),der1(j),der2(j)
        ! conditionals for factor 1
        call cfrkderivs(uvec(j),xl(iq),th(j),ccdf(j),cder1(j),cder2(j))
        ! factor 2 derivatives
        call frk2derivs(ccdf(j),xl(iq2),gam(j),gpdf(j),gder1(j),gder2(j),gder1u(j),gder2u(j),gdermix(j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=product(tpdf)*product(gpdf)
      ! fval1: vector of partial deriv wrt param[j], j=1,...,npar=2*d
      ! fval2: matrix of 2nd partial deriv wrt param[j], param[j2], 
      fval1=0.d0; fval2=0.d0
      do j=1,d
        fval1(j+d)=fval*gder1(j)/gpdf(j)
        fval1(j)=fval*( tder1(j)/tpdf(j) + cder1(j)*gder1u(j)/gpdf(j) )
        fval2(j+d,j+d)=fval*gder2(j)/gpdf(j)
        fval2(j,j)=fval*( tder2(j)/tpdf(j) +  &
            2.d0*tder1(j)*cder1(j)*gder1u(j)/tpdf(j)/gpdf(j) + &
            cder1(j)*cder1(j)*gder2u(j)/gpdf(j) + &
            cder2(j)*gder1u(j)/gpdf(j)  )
        ! cross th , g terms for fval2
        fval2(j,j+d)=fval* ( tder1(j)*gder1(j)/tpdf(j) + &
               cder1(j)*gdermix(j) ) /  gpdf(j)    
        fval2(j+d,j) = fval2(j,j+d)
      end do
      do j=1,d
        do j2=1,d
          if(j2 /= j) then
            fval2(j+d,j2+d)=fval1(j+d)*gder1(j2)/gpdf(j2)
            fval2(j+d,j2)=fval1(j2)*gder1(j)/gpdf(j)
            fval2(j,j2)=fval1(j2)*(tder1(j)/tpdf(j)+cder1(j)*gder1u(j)/gpdf(j))
            fval2(j,j2+d)=fval1(j)*gder1(j2)/gpdf(j2)
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
        hss(j,j2)=integl2(j,j2)/integl - grd(j)*grd(j2); ! any outer in f90? 
      end do
    end do
    hess=hess-hss
  end do
  deallocate (uvec, tpdf, tder1, tder2, fval1, fval2, integl1, integl2, grd, hss)
  deallocate (th, gam, gpdf, gder1, gder2, gder1u, gder2u, gdermix)
  deallocate (ccdf, cder1, cder2)
  return
  end


! Function with Frank pdf, der1, der2 (partial wrt cpar, 1st and 2nd order)
!  for factor 1.
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (-oo,oo)
! outputs 
!   pdf = copula density
!   der1 = \p pdf/\p cpar
!   der2 = \p^2 pdf/\p cpar^2
subroutine frk1derivs(u1,u2,cpar,pdf,der1,der2)
  implicit none
  double precision u1,u2,cpar,pdf,der1,der2
  double precision t00,t0,t1,t2,tem,num,den,den2,den3,dden,dnum,d2den,d2num
  t0 = exp(-cpar); t00=1.d0-t0;
  t1 = exp(-u1*cpar); t2 = exp(-u2*cpar);
  pdf=cpar*t1*t2*t00;
  tem=t00-(1.d0-t1)*(1.d0-t2);
  pdf=pdf/(tem*tem);
  num = cpar*(1.d0-t0)*t1*t2;
  den = -t0+t1+t2-t1*t2; den2=den*den; den3=den2*den
  dden = t0-u1*t1-u2*t2+(u1+u2)*t1*t2;
  dnum = t1*t2*(1.d0-t0-cpar*(u1+u2)+cpar*(u1+u2+1)*t0);
  der1 = dnum/den2 - 2.d0*num*dden/den3;   
  d2den = -t0+u1*u1*t1+u2*u2*t2-(u1+u2)**2*t1*t2; ! power is ** not ^
  d2num = t1*t2*(2.d0*t0*(u1+u2+1.d0)+cpar*(u1+u2)**2-2.d0*(u1+u2)-cpar*t0*(u1+u2+1.d0)**2);
  der2 = d2num/den2 - (2.d0*num*d2den + 4.d0*dnum*dden)/den3 + 6.d0*dden*dden*num/(den2*den2);
  return
  end

! Function with Frank pdf c(u1,u2,gam) for factor 2, and
!   der1, der2 (partial wrt gam, 1st and 2nd order)
!   also der1u der2u (partial wrt u1, 1st and 2nd order)
!   and dermix  (partial wrt gam and u1)
! inputs
!   u1,u2 = values in (0,1)
!   gam = scalar in (-oo,oo)
! outputs 
!   pdf = copula density
!   der1 = \p pdf/\p gam
!   der2 = \p^2 pdf/\p gam^2
!   der1u = \p pdf/\p u1
!   der2u = \p^2 pdf/\p u1^2
!   dermix = \p^2 pdf/\p u1 \p gam
subroutine frk2derivs(u1,u2,gam,pdf,der1,der2,der1u,der2u,dermix)
  implicit none
  double precision u1,u2,gam,pdf,der1,der2,der1u,der2u,dermix
  double precision t00,t0,t1,t2,tem,num,den,den2,den3,dden,dnum,d2den,d2num
  double precision duden,d2uden,gam2,dudden,tem1,tem2,den4
  t0 = exp(-gam); t00=1.d0-t0;
  t1 = exp(-u1*gam); t2 = exp(-u2*gam);
  pdf=gam*t1*t2*t00;
  tem=t00-(1.d0-t1)*(1.d0-t2);
  pdf=pdf/(tem*tem);
  num = gam*(1.d0-t0)*t1*t2;
  den = -t0+t1+t2-t1*t2; den2=den*den; den3=den2*den; den4=den2*den2;
  dden = t0-u1*t1-u2*t2+(u1+u2)*t1*t2;
  dnum = t1*t2*(1.d0-t0-gam*(u1+u2)+gam*(u1+u2+1.d0)*t0);
  der1 = dnum/den2 - 2.d0*num*dden/den3;   
  d2den = -t0+u1*u1*t1+u2*u2*t2-(u1+u2)**2*t1*t2; ! power is ** not ^
  d2num = t1*t2*(2.d0*t0*(u1+u2+1.d0)+gam*(u1+u2)**2-2.d0*(u1+u2)-gam*t0*(u1+u2+1.d0)**2);
  der2 = d2num/den2 - (2.d0*num*d2den + 4.d0*dnum*dden)/den3 + 6.d0*dden*dden*num/den4;
  ! derivatives wrt u1
  duden = -gam*t1*(1.d0-t2);
  der1u = -gam*(1.d0-t0)*t1*t2*(gam/den2+2.d0*duden/den3);
  gam2=gam*gam;
  d2uden = gam2*t1*(1.d0-t2);
  der2u = gam2/den2+(4.d0*gam*duden-2.d0*d2uden)/den3+6.d0*duden*duden/den4;
  der2u = der2u*gam*(1.d0-t0)*t1*t2;
  ! derivative wrt u1 and gam
  dudden = t1*(t2-1.d0+gam*u1-t2*gam*(u1+u2));
  tem1 = -(1.d0-t0+gam*t0-gam*(u1+u2)*(1.d0-t0))*(gam/den2+2.d0*duden/den3);
  tem2 = gam*(1.d0-t0)*(-1.d0/den2+2.d0*gam*dden/den3-2.d0*dudden/den3+6.d0*duden*dden/den4);
  dermix= t1*t2*(tem1+tem2);
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
subroutine cfrkderivs(u1,u2,cpar,ccdf,cder1,cder2)
  implicit none
  double precision u1,u2,cpar,ccdf,cder1,cder2
  double precision t00,t0,t1,t2,tem1,tem2,den,den2,den3,dden,d2den
  t0 = exp(-cpar); t00=1.d0-t0;
  t1 = exp(-u1*cpar); t2 = exp(-u2*cpar);
  ccdf=t2/(t00/(1.d0-t1)-(1.d0-t2));
  den = -t0+t1+t2-t1*t2;  den2=den*den; den3=den2*den
  dden = t0-u1*t1-u2*t2+(u1+u2)*t1*t2;
  cder1 = t2*(((u1+u2)*t1-u2)/den -(1.d0-t1)*dden/den2);
  d2den = -t0+u1*u1*t1+u2*u2*t2-(u1+u2)**2*t1*t2; 
  tem1 = (2.d0*dden*(u2-(u1+u2)*t1)-(1.d0-t1)*d2den)/den2;
  tem2 = (u2*u2-(u1+u2)**2*t1)/den+2.d0*(1.d0-t1)*dden*dden/den3;
  cder2 = t2*(tem1+tem2);
  return
  end

