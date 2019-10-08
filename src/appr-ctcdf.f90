

! conditional cdf of t copula with pt computed via monotone interpolation
! inputs
!   t1 = qt(u1,nu) transformed to t scale
!   t2 = qt(u2,nu) transformed to t scale
!   rho = scalar in (-1,1)
!   nu = degree of freedom >0
!   nipol = # points used for the interpolation
!   qq = quantiles of t(nu(1)+1)
!   pp = probabilities of t(nu(1)+1)
!   pder = derivs for the interpolation
! output
!   ccdf = conditional cdf
subroutine ctcdfipol(t1,t2,rho,nu, nipol,qq,pp,pder, ccdf)
  implicit none
  integer nipol,ierr ! added for interpol
  double precision qq(nipol),pp(nipol),pder(nipol),pdernew  ! added for interpol
  double precision t1,t2,rho,rho2,nu,ccdf
  double precision r2,cr,xt
  !double precision pt   ! function in C code ptbeta.c
  rho2 = rho*rho;
  r2 = sqrt((nu+t2*t2)/(nu+1));
  cr = sqrt(1-rho2)
  xt = (t1-rho*t2)/cr/r2;
  if(xt<qq(1)) xt=qq(1)
  if(xt>qq(nipol)) xt=qq(nipol)
  !ccdf = pt(xt,nu+1); 
  !print *,"1",ccdf
  ! replacement for interpol
  call pchev(nipol,qq,pp,pder, 1,xt,ccdf,pdernew,ierr)
  !print *,"2",ccdf
  return
  end

! function with condcdf T_{1|2}(t1|t2;rho,nu+1) for factor 1, and
! ccdf cder1, cder2 (partial wrt rho, 1st and 2nd order)
! where pt is computed via monotone interpolation
! inputs
!   t1 = qt(u1,nu) transformed to t scale
!   t2 = qt(u2,nu) transformed to t scale
!   rho = scalar in (-1,1)
!   nu = degree of freedom >0
!   nipol = # points used for the interpolation
!   qq = quantiles of t(nu(1)+1)
!   pp = probabilities of t(nu(1)+1)
!   pder = derivs for the interpolation
! output
!   ccdf = conditional cdf
!   cder1 = \p ccdf/\p rho
!   cder2 = \p^2 ccdf/\p rho^2
subroutine ctderivsipol(t1,t2,rho,nu, nipol,qq,pp,pder, ccdf,cder1,cder2)
  implicit none
  integer nipol,ierr ! added for interpol
  double precision qq(nipol),pp(nipol),pder(nipol),pdernew  ! added for interpol
  double precision t1,t2,rho,rho2,nu,ccdf,cder1,cder2
  double precision logpi,r2,cr,xt,tem1,tem2,const,gtem1,gtem2
  double precision lgamma   ! function in C code lgamma.c, in gnu fortran?
  !double precision pt   ! function in C code ptbeta.c
  logpi=1.1447298858494d0
  rho2 = rho*rho;
  r2 = sqrt((nu+t2*t2)/(nu+1));
  cr = sqrt(1-rho2)
  xt = (t1-rho*t2)/cr/r2;
  !ccdf = pt(xt,nu+1); 
  ! replacement for interpol
  call pchev(nipol,qq,pp,pder, 1,xt,ccdf,pdernew,ierr)
  tem1 = (t1*rho-t2)/cr**3/r2;
  tem2 = xt*xt/(nu+1);
  const = exp(lgamma(0.5*nu+1)-lgamma(0.5*nu+0.5)-0.5*logpi-0.5*log(nu+1));
  cder1 = const*(1+tem2)**(-0.5*nu-1)*tem1;
  gtem1 = (t1+2*t1*rho2-3*t2*rho)/cr**5/r2;
  gtem2 = 2*(rho*(t1*t1+t2*t2)-t1*t2*(rho2+1))/(nu+1)/cr**4/r2**2;
  cder2 = -const*(0.5*nu+1)*(1+tem2)**(-0.5*nu-2)*gtem2*tem1;
  cder2 = cder2+const*(1+tem2)**(-0.5*nu-1)*gtem1;
  return
  end

