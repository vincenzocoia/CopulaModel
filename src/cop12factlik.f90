! Frank copula log pdf only
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (-oo,oo) but could be underflow problems for th>35
! output 
!   lpdf = log pdf
subroutine lfrkpdf(u1,u2,cpar,lpdf)
  implicit none
  double precision u1,u2,cpar,lpdf
  double precision den,t0,t1,t2
  t0 = exp(-cpar);
  t1 = exp(-cpar*u1);
  t2 = exp(-cpar*u2);
  den = t1+t2-t0-t1*t2;
  lpdf = log(cpar*(1.d0-t0))-cpar*(u1+u2)-2.d0*log(abs(den));
  return
  end

! Frank copula conditional cdf only
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (-oo,oo)
! output
!   ccdf = conditional cdf
subroutine cfrkcdf(u1,u2,cpar,ccdf)
  implicit none
  double precision u1,u2,cpar,ccdf
  double precision t0,t1,t2,den,lccdf
  t0 = exp(-cpar);
  t1 = exp(-cpar*u1);
  t2 = exp(-cpar*u2);
  den = t1+t2-t0-t1*t2; 
  !lccdf = -cpar*u2+log(1.d0-t1)-log(den);
  ! correction for cpar >0 or <0
  lccdf = -cpar*u2+log((1.d0-t1)/den);
  ccdf = exp(lccdf)
  return
  end

! Gumbel copula log pdf only
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (1,oo)
! output 
!   lpdf = log pdf 
subroutine lgumpdf(u1,u2,cpar,lpdf)
  implicit none
  double precision u1,u2,cpar,lpdf
  double precision x,y,tx,ty,xd,yd,s,m,logm
  double precision den
  x = -log(u1); y = -log(u2);
  tx = log(x); ty = log(y); 
  xd = x**cpar; yd = y**cpar; s = xd+yd; m = s**(1.d0/cpar);
  den = m+cpar-1.d0; 
  logm=log(m); 
  lpdf = -m+log(den)+(1-2*cpar)*logm+(cpar-1)*(tx+ty)+x+y;
  return
  end

! Gumbel copula conditional cdf only
! inputs
!   u1,u2 = values in (0,1)
!   cpar = scalar in (1,oo)
! output
!   ccdf = conditional cdf
subroutine cgumcdf(u1,u2,cpar,ccdf)
  implicit none
  double precision u1,u2,cpar,ccdf
  double precision x,y,tx,ty,xd,yd,s,m
  double precision lccdf,logm
  x = -log(u1); y = -log(u2);
  tx = log(x); ty = log(y); 
  xd = x**cpar; yd = y**cpar; s = xd+yd; m = s**(1.d0/cpar);
  logm=log(m);
  lccdf = y-m+(1.d0-cpar)*(logm-ty)
  ccdf = exp(lccdf)
  return
  end


! t copula log pdf only
! inputs
!   t1 = qt(u1,nu) transformed to t scale
!   t2 = qt(u2,nu) transformed to t scale
!   rho = scalar in (-1,1)
!   nu = degree of freedom >0
! output 
!   lpdf = log pdf 
subroutine ltpdf(t1,t2,rho,nu,lpdf)
  implicit none
  double precision t1,t2,rho,nu,lpdf
  double precision pi,con,lgdif,coef,dt1,dt2,den
  double precision lgamma   ! function in C code lgamma.c , in gnu fortran?
  pi=3.14159265358979323846d0; ! pi
  con=0.5d0*log(pi*nu);
  lgdif=lgamma(0.5d0*(nu+1.d0))-lgamma(0.5d0*nu);
  coef = 1.d0-rho*rho; 
  den = 1.d0+(t1*t1-2.d0*rho*t1*t2+t2*t2)/(nu*coef);
  dt1 = lgdif-con-0.5d0*(nu+1.d0)*log(1.d0+t1*t1/nu);
  dt2 = lgdif-con-0.5d0*(nu+1.d0)*log(1.d0+t2*t2/nu);
  lpdf = -log(2.d0*pi)-0.5d0*log(coef)-0.5d0*(nu+2.d0)*log(den)-dt1-dt2;
  return
  end


! t copula conditional cdf only
! inputs
!   t1 = qt(u1,nu) transformed to t scale
!   t2 = qt(u2,nu) transformed to t scale
!   rho = scalar in (-1,1)
!   nu = degree of freedom >0
! output 
!   ccdf = conditional cdf 
subroutine ctcdf(t1,t2,rho,nu,ccdf)
  implicit none
  double precision t1,t2,rho,rho2,nu,ccdf
  double precision r2,cr,xt
  double precision pt   ! function in C code ptbeta.c
  rho2 = rho*rho;
  r2 = sqrt((nu+t2*t2)/(nu+1.d0));
  cr = sqrt(1.d0-rho2)
  xt = (t1-rho*t2)/cr/r2;
  ccdf = pt(xt,nu+1.d0); 
  return
  end

! BB1 copula log pdf only
! inputs
!   u1,u2 = values in (0,1)
!   param = (theta, delta), theta>0, delta>1
! output 
!   lpdf = log pdf
subroutine lbb1pdf(u1,u2,param,lpdf)
  implicit none
  double precision u1,u2,param(2),lpdf
  double precision theta,delta
  double precision t10,t20,t30,t40
  double precision m,td01,td02,s
  double precision den,coef,t1,t2,tu1,tu2
  double precision mp1,thtem,dltem,dl1,t1a,t2a,smlog
  theta=param(1); delta=param(2);
  t1 = u1**(-theta); t2 = u2**(-theta)
  t1a=t1-1.d0; t2a=t2-1.d0; smlog=log(t1a)+log(t2a)
  td01 = (t1a)**delta; td02 = (t2a)**delta
  s = td01+td02; m = s**(1.d0/delta)
  tu1 = -log(u1); tu2 = -log(u2)
  mp1=1.d0+m; 
  thtem=2.d0+1.d0/theta
  dltem=1.d0-2.d0*delta; dl1=delta-1.d0
  coef = theta*delta+1.d0
  den = theta*(dl1)+(coef)*m
  t10 = -(thtem)*log(mp1)
  t20 = (dltem)*log(m)
  t30 = log(den)
  t40 = (dl1)*(smlog)+(theta+1.d0)*(tu1+tu2)
  lpdf = t10+t20+t30+t40
  return
  end

! BB1 copula conditional cdf only
! inputs
!   u1,u2 = values in (0,1)
!   param = (theta, delta), theta>0, delta>1
! output
!   ccdf = conditional cdf
subroutine cbb1cdf(u1,u2,param,ccdf)
  implicit none
  double precision u1,u2,param(2),lcdf,ccdf
  double precision theta,delta
  double precision m,td01,td02,s
  double precision t1,t2,tu2,cf1,logm
  double precision mp1,dl1n,t1a,t2a,lt2a,lmp1
  theta=param(1); delta=param(2)
  t1 = u1**(-theta); t2 = u2**(-theta)
  t1a=t1-1.d0; t2a=t2-1.d0; 
  td01 = (t1a)**delta; td02 = (t2a)**delta
  s = td01+td02; m = s**(1.d0/delta)
  mp1=1.d0+m; lmp1=log(mp1); logm=log(m);
  tu2 = -log(u2);
  lt2a=log(t2a)
  cf1 = 1.d0+1.d0/theta;
  dl1n=1-delta
  lcdf = -cf1*lmp1+(dl1n)*(logm-lt2a)+(theta+1)*tu2;
  ccdf = exp(lcdf);
  return
  end
