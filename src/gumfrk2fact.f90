! Gumbel/Frank 2-factor copula model, Gumbel factor 1, Frank factor 2
! Subroutine lgum1derivs cgumderivs are in file gum12fact.f90
! Subroutine frk2derivs is in file frk12fact-pdf.f90

! inputs 
!   npar = #parameters,
!   param = parameter vector (dimension 2d=npar) theta>1, -oo<gam<oo
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
subroutine gumfrk2fact(npar,param,d,n,udata,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,nq
  double precision param(npar),udata(n,d),wl(nq),xl(nq)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,iq2
  double precision, dimension(:), allocatable :: uvec,tlpdf,tder1,tder2,fval1,integl1,grd
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision fval,integl,ww,w2
  double precision, dimension(:), allocatable :: th,gam,gpdf,gder1,gder2,gder1u,gder2u,gdermix
  double precision, dimension(:), allocatable :: ccdf,cder1,cder2,cdertem

  ! npar=2*d
  allocate ( uvec(d), tlpdf(d), tder1(d), tder2(d), fval1(npar), integl1(npar), grd(npar) )
  allocate ( th(d), gam(d), gpdf(d), gder1(d), gder2(d), gder1u(d), gder2u(d), gdermix(d) ) 
  allocate ( ccdf(d), cder1(d), cder2(d), cdertem(d) )
  allocate ( fval2(npar,npar), integl2(npar,npar), hss(npar,npar) )
  th(1:d)=param(1:d); gam(1:d)=param((d+1):npar)
  nllk=0.; grad=0.; hess=0.; ! initialize
  do i=1,n ! loop over rows of data set
    uvec=udata(i,:)
    integl=0.; integl1=0.; integl2=0.; ! initialize for integrals
    !print *, integl1
    do iq=1,nq ! loop over quadrature points
    do iq2=1,nq ! loop over quadrature points
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        ! factor 1 derivatives
        call lgum1derivs(uvec(j),xl(iq),th(j),tlpdf(j),tder1(j),tder2(j))
        ! conditionals for factor 1
        call cgumderivs(uvec(j),xl(iq),th(j),ccdf(j),cder1(j),cder2(j))
        !if(iq==1.and.iq2==1) print "(4f11.6)", th(j),ccdf(j),cder1(j),cder2(j)
        ! factor 2 derivatives
        call frk2derivs(ccdf(j),xl(iq2),gam(j),gpdf(j),gder1(j),gder2(j),gder1u(j),gder2u(j),gdermix(j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      fval=exp(sum(tlpdf))*product(gpdf)
      ! fval1: vector of partial deriv wrt param[j], j=1,...,npar=2*d
      ! fval2: matrix of 2nd partial deriv wrt param[j], param[j2], 
      fval1=0.; fval2=0.
      ! convert from derivative of log pdf for Gumbel to deriv of pdf ?
      do j=1,d
        cdertem(j)= tder1(j) + cder1(j)*gder1u(j)/gpdf(j) 
        fval1(j+d)=fval*gder1(j)/gpdf(j)
        fval1(j)=fval*cdertem(j)
        fval2(j+d,j+d)=fval*gder2(j)/gpdf(j)
        fval2(j,j)=fval*( tder1(j)*tder1(j) + tder2(j)   &
           +2.d0*tder1(j)*cder1(j)*gder1u(j)/gpdf(j) &
           +cder1(j)*cder1(j)*gder2u(j)/gpdf(j) &
           +cder2(j)*gder1u(j)/gpdf(j) )
        ! cross th , gam terms for fval2
        fval2(j,j+d)=fval* (tder1(j)*gder1(j)+ cder1(j)*gdermix(j) )/gpdf(j) 
        fval2(j+d,j) = fval2(j,j+d)
      end do
      do j=1,d
        do j2=1,d
          if(j2 /= j) then
            fval2(j+d,j2+d)=fval1(j+d)*gder1(j2)/gpdf(j2)
            fval2(j+d,j2)=fval1(j2)*gder1(j)/gpdf(j)
            fval2(j,j2)=fval1(j2)*cdertem(j) 
            fval2(j,j2+d)=fval1(j)*gder1(j2)/gpdf(j2)
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
  deallocate (th, gam, gpdf, gder1, gder2, gder1u, gder2u, gdermix)
  deallocate (ccdf, cder1, cder2,cdertem)
  return
  end

! function with lpdf, lder1, lder2 (partial wrt delta and 2nd order)
! subroutine lgum1derivs(u1,u2,delta,lpdf,lder1,lder2)

! function with log pdf c(u1,u2,delta) for factor 2 
! der1, der2 (partial wrt delta, 1st and 2nd order)
! also der1u der2u (partial wrt u1, 1st and 2nd order)
!   and dermix  (partial wrt delta and u1 
! subroutine lgum2derivs(u1,u2,delta,lpdf,lder1,lder2,lder1u,lder2u,ldermix)

! function with condcdf C_{1|2}(u1|u2;delta) for factor 1 
! ccdf cder1, cder2 (partial wrt delta, 1st and 2nd order)
! subroutine cgumderivs(u1,u2,delta,ccdf,cder1,cder2)

! function with pdf, der1, der2 (partial wrt theta, 1st and 2nd order)
! subroutine frk1derivs(u1,u2,theta,pdf,der1,der2)

! function with pdf c(u1,u2,gam) for factor 2 
! der1, der2 (partial wrt gam, 1st and 2nd order)
! also der1u der2u (partial wrt u1, 1st and 2nd order)
!   and dermix  (partial wrt gam and u1 
! subroutine frk2derivs(u1,u2,gam,pdf,der1,der2,der1u,der2u,dermix)

