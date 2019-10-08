! Student 1-factor and 2-factor copula models; df's are fixed 
! for simplicity, assume df is nu1 for all variables with factor1
!                          and nu2 for all variables with factor2
! nu1,nu2 cannot be less than 1, 
!   also version of qt here not accurate for 1<nu<2

! gfortran  t12fact.f90 ptbeta.o qtnorm.o
!program t2factchk
!  implicit none
!  integer n,d,nq,i,iq,j,npar
!  double precision, dimension(:,:), allocatable :: tdata,hess
!  double precision, dimension(:), allocatable :: tl,wl,param,grad
!  double precision nllk,nu1,nu2
!  read *,nq
!  allocate ( tl(nq),wl(nq) )
!  do iq =1,nq   ! equidistant and equally weighted for testing
!    tl(iq)=iq/(nq+1.d0)
!    wl(iq)=1.d0/nq
!  end do
!  read *,n,d
!  npar=d*2
!  allocate ( tdata(n,d), param(npar), grad(npar), hess(npar,npar) )
!!  do i=1,n
!    read *, tdata(i,:)
!  end do
!  do j=1,d
!    param(j)=j*0.12d0
!    param(j+d)=j*0.1d0
!    if(param(j)>=1.) param(j)=0.5d0
!  end do
!  nu1=3.d0; nu2=4.d0
!  print "(10f8.2)", param, nu1,nu2
!  !call t1fact(npar,param,nu1,d,n,tdata,nq,wl,tl,nllk,grad,hess)
!  call t2fact(npar,param,nu1,nu2,d,n,tdata,nq,wl,tl,tl,nllk,grad,hess)
!  print *, nllk
!  print "(8f10.6)", grad
!  print "(8f10.6)", hess
!  deallocate (tl, wl, tdata, grad, hess)
!  stop
!  end

! 1-factor model with t(nu)
! tl() have been transformed to t(nu) scale = qt(gl$nodes,nu)
! tdata have been transformed to t(nu) scale = qt(udata,nu)
! inputs 
!   npar = #parameters,
!   rho = parameter vector (dimension d=npar), -1<rho[j]<1
!   nu = degree of freedom >0
!   d = #variables
!   n = sample size
!   tdata = nxd matrix of transformed U(0,1) data
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   tl = vector of transformed quadrature nodes
! outputs 
!   nllk = negative log-likelihood, 
!   grad = gradient of nllk, 
!   hess = hessian of nllk
subroutine t1fact(npar,rho,nu,d,n,tdata,nq,wl,tl,nllk,grad,hess)
  implicit none
  integer npar,d,n,nq
  double precision rho(d),nu,tdata(n,d),wl(nq),tl(nq) 
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2
  double precision, dimension(:), allocatable :: tvec,lpdf,lder1,lder2,fval1,integl1,grd
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision fval,integl,ww,txl

  ! npar=d
  allocate ( tvec(d), lpdf(d), lder1(d), lder2(d), fval1(d), integl1(d), grd(d) )
  allocate ( fval2(d,d), integl2(d,d), hss(d,d) )
  nllk=0.d0; grad=0.d0; hess=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    tvec=tdata(i,:)
    !print "(i2,10f7.3)", i, tvec
    integl=0.d0; integl1=0.d0; integl2=0.d0; ! initialize for integrals
    do iq=1,nq ! loop over quadrature points
      ! get \p log copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      txl = tl(iq);
      do j=1,d
        call lt1derivs(txl,tvec(j),rho(j),nu,lpdf(j),lder1(j),lder2(j))
        !if(iq==1) print "(3f10.6)", lpdf(j),lder1(j),lder2(j)
        end do
      ! update integrand, and derivs wrt copula parameters 
      fval=exp(sum(lpdf))
      ! fval1: vector of partial deriv wrt theta[j], j=1,...,d
      ! fval2: matrix of 2nd partial deriv wrt theta[j], theta[j2], 
      do j=1,d
        fval1(j)=fval*lder1(j)
        fval2(j,j)=fval*(lder2(j)+lder1(j)*lder1(j))
        do j2=1,d
          if(j2 /= j) fval2(j,j2)=fval1(j)*lder1(j2)
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
  deallocate (tvec, lpdf, lder1, lder2, fval1, fval2, integl1, integl2, grd, hss)
  return
  end

! 2-factor copula with t
! inputs 
!   npar = #parameters,
!   param = parameter vector of rhos (dimension 2d=npar)
!   nu1 = degree of freedom >0
!   nu2 = degree of freedom >0
!   d = #variables
!   n = sample size
!   tdata = nxd matrix of transformed U(0,1) data : t(nu1) scale
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   tl1 = vector of transformed quadrature nodes wrt nu1
!   tl2 = vector of transformed quadrature nodes wrt nu2
! outputs 
!   nllk = negative log-likelihood, 
!   grad = gradient of nllk, 
!   hess = hessian of nllk
subroutine t2fact(npar,param,nu1,nu2,d,n,tdata,nq,wl,tl1,tl2,nllk,grad,hess)
  implicit none
  integer npar,d,n,nq
  double precision param(npar),nu1,nu2,tdata(n,d),wl(nq),tl1(nq),tl2(nq)
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
  nllk=0.d0; grad=0.d0; hess=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    tvec=tdata(i,:)
    integl=0.d0; integl1=0.d0; integl2=0.d0; ! initialize for integrals
    do iq=1,nq ! loop over quadrature points
    do iq2=1,nq ! loop over quadrature points
      ! get \p copden(u[j],v,th[j]) / \p th[j] , j=1,...,d
      do j=1,d
        ! factor 1 derivatives
        call lt1derivs(tvec(j),tl1(iq),th(j),nu1,tpdf(j),tder1(j),tder2(j))
        ! print "(3f10.6)", pdfv(j),der1(j),der2(j)
        ! conditionals for factor 1
        call ctderivs(tvec(j),tl1(iq),th(j),nu1,ccdf(j),cder1(j),cder2(j))
        ! factor 2 derivatives
        qtccdf = qt(ccdf(j),nu2); ! **Student quantile is needed here**  
        call lt2derivs(qtccdf,tl2(iq2),gam(j),nu2,gpdf(j),gder1(j),gder2(j),gder1u(j),gder2u(j),gdermix(j))
      end do
      ! update integrand, and derivs wrt copula parameters 
      !fval=product(tpdf)*product(gpdf)
      fval=exp(sum(tpdf)+sum(gpdf))
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

! Function with t lpdf, lder1, lder2 (partial wrt rho and 2nd order)
! inputs
!   t1 = qt(u1,nu) transformed to t scale
!   t2 = qt(u2,nu) transformed to t scale
!   rho = scalar in (-1,1)
!   nu = degree of freedom >0
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p rho
!   lder2 = \p^2 lpdf/\p rho^2
subroutine lt1derivs(t1,t2,rho,nu,lpdf,lder1,lder2)
  implicit none
  double precision t1,t2,rho,nu,lpdf,lder1,lder2
  double precision pi,con,lgdif,coef,coef2,dt1,dt2,den,den1,den2
  double precision lgamma   ! function in C code lgamma.c, in gnu fortran?
  pi=3.14159265358979323846d0; ! pi
  con=0.5d0*log(pi*nu);
  !lgdif=log(gamma(0.5*(nu+1)))-log(gamma(0.5*nu));
  lgdif=lgamma(0.5d0*(nu+1.d0))-lgamma(0.5d0*nu);
  coef = 1.d0-rho*rho; coef2 = coef*coef;
  den = 1.d0+(t1*t1-2.d0*rho*t1*t2+t2*t2)/(nu*coef);
  den1 = 2.d0*rho*(den-1.d0)/coef-2.d0*t1*t2/(nu*coef);
  den2 = 2.d0*rho*den1/coef+2.d0*(den-1.d0)*(2.d0-coef)/coef2-4.d0*rho*t1*t2/(nu*coef2);
  dt1 = lgdif-con-0.5d0*(nu+1.d0)*log(1.d0+t1*t1/nu);
  dt2 = lgdif-con-0.5d0*(nu+1.d0)*log(1.d0+t2*t2/nu);
  lpdf = -log(2.d0*pi)-0.5d0*log(coef)-0.5d0*(nu+2.d0)*log(den)-dt1-dt2;
  lder1 = rho/coef-0.5d0*(nu+2.d0)*den1/den;
  lder2 = (2-coef)/coef2-0.5d0*(nu+2.d0)*(den2/den-den1*den1/(den*den));
  return
  end

! Function with t log pdf = log f(t1,t2,rho,nu) for factor 2, and
!   lder1, lder2 (partial wrt rho, 1st and 2nd order)
!   also lder1u lder2u (partial wrt t1, 1st and 2nd order)
!   and ldermix  (partial wrt rho and t1 )
! inputs
!   t1 = qt(u1,nu) transformed to t scale
!   t2 = qt(u2,nu) transformed to t scale
!   rho = scalar in (-1,1)
!   nu = degree of freedom >0
! outputs 
!   lpdf = log pdf, 
!   lder1 = \p lpdf/\p rho
!   lder2 = \p^2 lpdf/\p rho^2
!   ltder1u = \p lpdf/\p t1
!   ltder2u = \p^2 lpdf/\p t1^2
!   ltdermix = \p^2 lpdf/\p u1 \p rho
subroutine lt2derivs(t1,t2,rho,nu,lpdf,lder1,lder2,ltder1u,ltder2u,ltdermix)
  implicit none
  double precision t1,t1sq,t2,t2sq,rho,nu,lpdf,lder1,lder2,ltder1u,ltder2u,ltdermix
  double precision pi,con,lgdif,coef,coef2,coef3,dt1,dt2,den,den1,den2
  double precision dd1,dd12,ltder1t1,ltder2t1,ltdert1mx, reg
  double precision lgamma   ! function in C code lgamma.c, in gnu fortran

  pi=3.14159265358979323846d0; ! pi
  con=0.5d0*log(pi*nu);
  !lgdif=log(gamma(0.5*(nu+1)))-log(gamma(0.5*nu));
  lgdif=lgamma(0.5d0*(nu+1.d0))-lgamma(0.5d0*nu);
  coef = 1.d0-rho*rho; coef2 = coef*coef; 
  t1sq = t1*t1; t2sq = t2*t2;
  den = 1.d0+(t1sq-2.d0*rho*t1*t2+t2sq)/(nu*coef);
  den1 = 2.d0*rho*(den-1.d0)/coef-2.d0*t1*t2/(nu*coef);
  den2 = 2.d0*rho*den1/coef+2.d0*(den-1.d0)*(2.d0-coef)/coef2-4.d0*rho*t1*t2/(nu*coef2);
  dt1 = lgdif-con-0.5d0*(nu+1.d0)*log(1.d0+t1sq/nu);
  dt2 = lgdif-con-0.5d0*(nu+1.d0)*log(1.d0+t2sq/nu);
  dd1 = exp(dt1); dd12 = dd1*dd1; 
  reg = t1-rho*t2;  ! added temporary variable
  coef3 = den*coef*nu; ! moved
  lpdf = -log(2.d0*pi)-0.5d0*log(coef)-0.5d0*(nu+2.d0)*log(den)-dt1-dt2;
  lder1 = rho/coef-0.5d0*(nu+2.d0)*den1/den;
  lder2 = (2.d0-coef)/coef2-0.5d0*(nu+2.d0)*(den2/den-den1*den1/(den*den));
  !ltder1t1 = -(nu+2)*(t1-rho*t2)/coef3 + (nu+1)*t1/(nu+t1sq);
  !ltder2t1 = -(nu+2)/coef3+2*(nu+2)*(t1-rho*t2)**2/coef3**2+(nu+1)*(nu-t1sq)/(nu+t1sq)**2;
  !ltdert1mx = (nu+2)*(t2/coef3-2*rho*(t1-rho*t2)/coef3/coef+(t1-rho*t2)*den1/den/coef3);
  ltder1t1 = -(nu+2.d0)*reg/coef3 + (nu+1.d0)*t1/(nu+t1sq);
  ltder2t1 = -(nu+2.d0)/coef3+2.d0*(nu+2.d0)*(reg*reg)/coef3**2+(nu+1.d0)*(nu-t1sq)/(nu+t1sq)**2;
  ltdert1mx = (nu+2.d0)*(t2/coef3-2.d0*rho*reg/coef3/coef+reg*den1/den/coef3);
  ltder1u =ltder1t1/dd1;
  ltder2u =ltder2t1/dd12+(nu+1.d0)*t1*ltder1t1/dd12/(nu+t1sq);
  ltdermix=ltdert1mx/dd1;
  return
  end


! function with t condcdf T_{1|2}(t1|t2;rho,nu+1) for factor 1, and 
! ccdf cder1, cder2 (partial wrt rho, 1st and 2nd order)
! inputs
!   t1 = qt(u1,nu) transformed to t scale
!   t2 = qt(u2,nu) transformed to t scale
!   rho = scalar in (-1,1)
!   nu = degree of freedom >0
! outputs 
!   ccdf = conditional cdf, 
!   cder1 = \p ccdf/\p rho
!   cder2 = \p^2 ccdf/\p rho^2
subroutine ctderivs(t1,t2,rho,nu,ccdf,cder1,cder2)
  implicit none
  double precision t1,t2,rho,rho2,nu,ccdf,cder1,cder2
  double precision logpi,r2,cr,xt,tem1,tem2,const,gtem1,gtem2
  double precision pt   ! function in C code ptbeta.c
  double precision lgamma   ! function in C code lgamma.c, in gnu fortran
  logpi=1.1447298858494d0
  rho2 = rho*rho;
  r2 = sqrt((nu+t2*t2)/(nu+1.d0));
  cr = sqrt(1.d0-rho2)
  xt = (t1-rho*t2)/cr/r2;
  ccdf = pt(xt,nu+1.d0); 
  tem1 = (t1*rho-t2)/cr**3/r2;
  tem2 = xt*xt/(nu+1.d0);
  const = exp(lgamma(0.5d0*nu+1.d0)-lgamma(0.5d0*nu+0.5d0)-0.5d0*logpi-0.5d0*log(nu+1.d0));
  cder1 = const*(1.d0+tem2)**(-0.5d0*nu-1.d0)*tem1;
  gtem1 = (t1+2.d0*t1*rho2-3.d0*t2*rho)/cr**5/r2;
  gtem2 = 2.d0*(rho*(t1*t1+t2*t2)-t1*t2*(rho2+1.d0))/(nu+1.d0)/cr**4/r2**2;
  cder2 = -const*(0.5d0*nu+1.d0)*(1.d0+tem2)**(-0.5d0*nu-2.d0)*gtem2*tem1;
  cder2 = cder2+const*(1.d0+tem2)**(-0.5d0*nu-1.d0)*gtem1;
  return
  end

