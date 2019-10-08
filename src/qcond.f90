 
! C_{2|1}(v|u;cpar) for Gumbel
! input
!   u,v = values in (0,1)
!   cpar = copula parameter >1 
! output 
!   conditional cdf
double precision function pcondgumbel(v,u,cpar)
  implicit none
  double precision u,v,x,y,tem1,tem2,sm,tem,ccdf,cpar;
  x= -log(u); y= -log(v);
  tem1=x**cpar; tem2=y**cpar; sm=tem1+tem2; tem=sm**(1.d0/cpar);
  ccdf=exp(-tem);
  ccdf=ccdf*(1.+tem2/tem1)**(-1.d0+1.d0/cpar);
  ccdf=ccdf/u;
  pcondgumbel=ccdf;
  return
  end


! G(x,y)=exp(-(x^cpar+y^cpar)^(1/cpar)), cpar>=1
!   conditional P(Y>=y|x)=exp(x)*G(x,y)* (x^cpar+y^cpar)^(1/cpar-1) * x^(cpar-1)
!                           = p
!   take logs, and let z=(x^cpar+y^cpar)^(1/cpar) >=x
!   g(z) = log(p) - x + z + (cpar-1)*log(z) -(cpar-1)*log(x) =0
!   g'(z) = 1 - (cpar-1)/z
!   solve for z with NR, then solve for y=y(z,x,p,cpar)
!   y = (z^cpar - x^cpar)^(1/cpar)
!
! C_{2|1}^{-1}(p|u;cpar)
! input
!   p,u = values in (0,1)
!   cpar = copula parameter >1 
! output 
!   inverse of conditional cdf
double precision function qcondgumbel(p, u, cpar)
  implicit none
  double precision p,u,cpar,z,g,gp,x,y,con,cpar1,dif,mxdif,eps
  integer iter,mxiter

  mxiter=20; eps=1.d-6;
  x=-log(u); cpar1=cpar-1.;
  con=log(p)-x-cpar1*log(x); 
  z=x*(2.d0**(1.d0/cpar))
  mxdif=1; iter=0; 
  dif=.1d0;  ! needed in case first step leads to NaN
  do while(mxdif>eps .and. iter<mxiter)
    g=z+cpar1*log(z)+con;
    gp=1.d0+cpar1/z;
    ! if(isnan(g) || isnan(gp) || isnan(g/gp) ) { dif=dif/(-2.); }  // added for cpar>50
    !if(g==NaN .or. gp==NaN .or. g/gp==NaN ) then
    !  dif=dif/(-2.) ! added for cpar>50
    !else 
      dif=g/gp;
    !end if
    z=z-dif; iter=iter+1;
    do while(z<=x) 
      dif=dif/2.d0; z=z+dif; 
    end do
    !print *, iter, dif, z;
    mxdif=abs(dif);
  end do
  if(iter>=mxiter)  then
    !print *, "***did not converge"
    call intpr("***did not converge",19,0,0)
    !print "(4f10.6)",  p,x,cpar,z
  end if
  y=((z**cpar)-(x**cpar))**(1.d0/cpar)
  qcondgumbel=exp(-y)
  return
  end
