! Student 1-factor and 2-factor copula models; df's are fixed 
! for simplicity, assume df is nu1 for all variables with factor1
!                          and nu2 for all variables with factor2
! nu1,nu2 cannot be less than 1, 
!   also version of qt here not accurate for 1<nu<2

! 2-factor copula with t and monotone interpolation for t cdf
! The subroutines lt1derivs lt2derivs are in file t12fact.f90
! The subroutine  ctderivsipol is in file appr-ctcdf.f90

! tl1() have been transformed to t(nu1) scale = qt(gl$nodes,nu1)
! tl2() have been transformed to t(nu2) scale = qt(gl$nodes,nu2)
! tdata have been transformed to t(nu) scale = qt(udata,nu)
! inputs 
!   npar = #parameters,
!   param = parameter vector of rhos (dimension 2d=npar), -1<rho[j]<1
!   nu1 = degree of freedom >0
!   nu2 = degree of freedom >0
!   d = #variables
!   n = sample size
!   tdata = nxd matrix of transformed U(0,1) data : t(nu1) scale
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   tl1 = vector of transformed quadrature nodes wrt nu1
!   tl2 = vector of transformed quadrature nodes wrt nu2
!   nipol = # points used for the interpolation
!   qq = quantiles of t(nu(1)+1)
!   pp = probabilities of t(nu(1)+1)
!   pder = derivs for the interpolation
! outputs 
!   nllk = negative log-likelihood, 
!   grad = gradient of nllk, 
!   hess = hessian of nllk
subroutine t2ipol(npar,param,nu1,nu2,d,n,tdata,nq,wl,tl1,tl2, &
  nipol,qq,pp,pder, nllk,grad,hess)
  implicit none
  integer npar,d,n,nq, nipol ! added for interpol
  double precision param(npar),nu1,nu2,tdata(n,d),wl(nq),tl1(nq),tl2(nq)
  double precision qq(nipol),pp(nipol),pder(nipol)  ! added for interpol
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,iq2
  double precision, dimension(:), allocatable :: tvec,tpdf,tder1,tder2,fval1,integl1,grd
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision fval,integl,ww,w2,qtccdf
  double precision, dimension(:), allocatable :: th,gam,gpdf,gder1,gder2,gder1u,gder2u,gdermix
  double precision, dimension(:), allocatable :: ccdf,cder1,cder2,cdertem
  double precision qt   ! function in f90 code qtnorm.f90

  ! npar=2*d
  allocate ( tvec(d), tpdf(d), tder1(d), tder2(d), fval1(npar), integl1(npar), grd(npar) )
  allocate ( th(d), gam(d), gpdf(d), gder1(d), gder2(d), gder1u(d), gder2u(d), gdermix(d) ) 
  allocate ( ccdf(d), cder1(d), cder2(d), cdertem(d) )
  allocate ( fval2(npar,npar), integl2(npar,npar), hss(npar,npar) )
  th(1:d)=param(1:d); gam(1:d)=param((d+1):npar)
  nllk=0.; grad=0.; hess=0.; ! initialize
  do i=1,n ! loop over rows of data set
    tvec=tdata(i,:)
    integl=0.; integl1=0.; integl2=0.; ! initialize for integrals
    do iq=1,nq ! loop over quadrature points
    do iq2=1,nq ! loop over quadrature points
      ! get \p copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        ! factor 1 derivatives
        call lt1derivs(tvec(j),tl1(iq),th(j),nu1,tpdf(j),tder1(j),tder2(j))
        ! print "(3f10.6)", pdfv(j),der1(j),der2(j)
        ! conditionals for factor 1
        ! call ctderivs(tvec(j),tl1(iq),th(j),nu1,ccdf(j),cder1(j),cder2(j))
        ! modified function for interpol (in appr-ctcdf.f90)
        call ctderivsipol(tvec(j),tl1(iq),th(j),nu1, &
            nipol,qq,pp,pder, ccdf(j),cder1(j),cder2(j))  !C_{ij|V_0}
        if(ccdf(j) < 0.00001) ccdf(j) = 0.00001
        if(ccdf(j) > 0.99999) ccdf(j) = 0.99999
        ! factor 2 derivatives
        qtccdf = qt(ccdf(j),nu2); ! **Student quantile is needed here**  
        call lt2derivs(qtccdf,tl2(iq2),gam(j),nu2,gpdf(j),gder1(j),gder2(j),gder1u(j),gder2u(j),gdermix(j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      !fval=product(tpdf)*product(gpdf)
      fval=exp(sum(tpdf)+sum(gpdf))
      ! fval1: vector of partial deriv wrt param[j], j=1,...,npar=2*d
      ! fval2: matrix of 2nd partial deriv wrt param[j], param[j2], 
      fval1=0.; fval2=0.
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
  deallocate (tvec, tpdf, tder1, tder2, fval1, fval2, integl1, integl2, grd, hss)
  deallocate (th, gam, gpdf, gder1, gder2, gder1u, gder2u, gdermix)
  deallocate (ccdf, cder1, cder2, cdertem)
  return
  end

