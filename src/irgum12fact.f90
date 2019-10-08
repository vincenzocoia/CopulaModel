! Gumbel 1-factor and 2-factor copula models for discrete item response data

! gfortran -o irgum12fact irgum12fact.f90

! sample main
!program irgum
!  implicit none
!  integer n,d,nq,i,iq,j,k,npar,ncat
!  integer, dimension(:,:), allocatable :: ydata
!  double precision, dimension(:,:), allocatable :: hess,ucuts
!  double precision, dimension(:), allocatable :: xl,wl,param,grad
!  double precision nllk
!  read *,nq
!  allocate ( xl(nq),wl(nq) )
!  do iq =1,nq   ! equidistant and equially weighted for testing
!    xl(iq)=iq/(nq+1.d0)
!    wl(iq)=1.d0/nq
!  end do
!  read *,n,d,ncat
!  print *, "n=", n, " nitems=", d, " ncat=", ncat
!  npar=d
!  allocate ( ydata(n,d), param(npar), grad(npar), hess(npar,npar), ucuts(ncat-1,d) )
!  do i=1,n
!    read *, ydata(i,:)
!  end do
!  do j=1,d
!    do k=1,ncat-1
!      ucuts(k,j)=(j+k)/(ncat+d+1.d0)
!    end do
!  end do
!  ! 1-factor
!  do j=1,d
!    param(j)=1.d0+.5d0*j
!  end do
!  print *, "1-factor"
!  print *, param
!  call irgum1fact(npar,param,d,ncat,n,ydata,ucuts,nq,wl,xl,nllk,grad,hess)
!  print *, nllk
!  print "(8f10.6)", grad
!  print "(8f10.6/)", hess
!  deallocate (param,grad,hess)
!  npar=2*d
!  allocate ( param(npar), grad(npar), hess(npar,npar) )
!  ! 2-factor
!  do j=1,d
!    param(j)=1.d0+.5d0*j
!    param(j+d)=1.d0+.1d0*j ! second factor
!  end do
!  print *, "2-factor"
!  print *, param
!  call irgum2fact(npar,param,d,ncat,n,ydata,ucuts,nq,wl,xl,nllk,grad,hess)
!  print *, nllk
!  print "(8f10.6)", grad
!  print "(8f10.6/)", hess
!  deallocate (xl, wl, ydata, param, grad, hess, ucuts)
!  stop
!  end


! 1-factor model with Gumbel, 
! inputs 
!   npar = #parameters,
!   th = parameter vector (dimension d=npar), th[j]>1
!   d = #items
!   ncat = #categories (ordinal)
!   n = sample size
!   ydata = nxd matrix of integer data, each in {0,1,...,ncat-1}
!   ucuts = (ncat-1)xd matrix of cutpoints on Uniform(0,1) scale
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes
! outputs 
!   nllk = negative log-likelihood, 
!   grad = gradient of nllk, 
!   hess = hessian of nllk
subroutine irgum1fact(npar,th,d,ncat,n,ydata,ucuts,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,ncat,nq
  integer ydata(n,d)
  double precision th(npar),wl(nq),xl(nq),ucuts(ncat-1,d)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,k
  integer, dimension(:), allocatable :: yvec
  double precision, dimension(:), allocatable :: pmf,der1,der2,fval1,integl1,grd
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision, dimension(:,:,:), allocatable :: cond,cond1,cond2
  double precision fval,integl,ww

  ! npar=d
  allocate ( yvec(d), pmf(d), der1(d), der2(d), fval1(d), integl1(d), grd(d) )
  allocate ( fval2(d,d), integl2(d,d), hss(d,d) )
  allocate ( cond(nq,ncat,d), cond1(nq,ncat,d), cond2(nq,ncat,d) )
  ! cond for cond cdf and pmf, cond1 for deriv wrt theta, cond2 for 2nd deriv
   
  !print *, "n=", n, " nitems=", d, " ncat=", ncat
  !print *, "nq=", nq
  !print "(8f10.6/)", wl
  cond=0.d0; cond1=0.d0; cond2=0.d0
  cond(:,ncat,:)=1.d0  ! condcdf is 1 for last category ncat
  ! get \p copcond(ucuts[k]|v,th[j]) / \p th[j] , k=1,...,ncat ; v=latent
  do iq=1,nq ! loop over quadrature points
    do j=1,d ! loop over items
      do k=1,ncat-1 ! loop of categories
        call irgum1derivs(xl(iq),ucuts(k,j),th(j),cond(iq,k,j),cond1(iq,k,j),cond2(iq,k,j))
      end do
    end do
  end do
  ! differencing
  do iq=1,nq 
    do j=1,d
      do k=ncat,2,-1
        cond(iq,k,j)=cond(iq,k,j)-cond(iq,k-1,j)
        cond1(iq,k,j)=cond1(iq,k,j)-cond1(iq,k-1,j)
        cond2(iq,k,j)=cond2(iq,k,j)-cond2(iq,k-1,j)
      end do
    end do
  end do
  !print "(8f10.6)", cond(1,:,d)
  !print "(8f10.6)", cond1(1,:,d)
  !print "(8f10.6)", cond2(1,:,d)
  ! looping through data
  nllk=0.d0; grad=0.d0; hess=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    yvec=ydata(i,:)+1
    !if(i==1) print "(8i3)", yvec
    integl=0.d0; integl1=0.d0; integl2=0.d0; ! initialize for integrals
    do iq=1,nq 
      ! update integrand, and derivs wrt copula parameters 
      do j=1,d
        pmf(j)=cond(iq,yvec(j),j)
        der1(j)=cond1(iq,yvec(j),j)
        der2(j)=cond2(iq,yvec(j),j)
      end do
      fval=product(pmf)
      !if(iq==1 .and. i==1) print "(4f10.6)", pmf(1),der1(1),der2(1)
      !if(iq==1 .and. i==1) print "(8f10.6)", pmf,fval
      ! fval1: vector of partial deriv wrt theta[j], j=1,...,d
      ! fval2: matrix of 2nd partial deriv wrt theta[j], theta[j2], 
      do j=1,d
        fval1(j)=fval*der1(j)/pmf(j)
        fval2(j,j)=fval*der2(j)/pmf(j)
        do j2=1,d
          if(j2 /= j) fval2(j,j2)=fval1(j)*der1(j2)/pmf(j2)
        end do
      end do
      ! update quadrature loops
      ww=wl(iq)
      integl=integl+fval*ww
      integl1=integl1+fval1*ww
      integl2=integl2+fval2*ww
    end do
    !if(i==1) print "(8f10.6/)",  integl, integl1
    ! update contribution to negative log-likelihood
    nllk=nllk-log(integl)
    grd=integl1/integl; grad=grad-grd;
    do j=1,d
      do j2=1,d
        hss(j,j2)=integl2(j,j2)/integl - grd(j)*grd(j2); 
      end do
    end do
    hess=hess-hss
  end do
  deallocate (yvec, pmf, der1, der2, fval1, fval2, integl1, integl2, grd, hss)
  deallocate (cond, cond1, cond2)
  return
  end


! 2-factor model with Gumbel
! inputs 
!   npar = #parameters,
!   param =  parameter vector (dimension 2d=npar)
!   d = #items
!   ncat = #categories (ordinal)
!   n = sample size
!   ydata = nxd matrix of integer data, each in {0,1,...,ncat-1}
!   ucuts = (ncat-1)xd matrix of cutpoints on Uniform(0,1) scale
!   nq = #quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes
! outputs 
!   nllk = negative log-likelihood, 
!   grad = gradient of nllk, 
!   hess = hessian of nllk
subroutine irgum2fact(npar,param,d,ncat,n,ydata,ucuts,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,d,n,ncat,nq
  integer ydata(n,d)
  double precision param(npar),wl(nq),xl(nq),ucuts(ncat-1,d)
  double precision nllk,grad(npar),hess(npar,npar)
  integer i,iq,j,j2,k,iq2
  integer, dimension(:), allocatable :: yvec
  double precision, dimension(:), allocatable :: pmf,tder1,tder2,fval1,integl1,grd
  double precision, dimension(:), allocatable :: th,gam,gder1,gder2,tgder
  double precision, dimension(:,:), allocatable :: fval2,integl2,hss
  double precision, dimension(:,:,:,:), allocatable :: cond,tcond1,tcond2
  double precision, dimension(:,:,:,:), allocatable :: gcond1,gcond2,tgcond2
  double precision fval,integl,ww,w2
  ! t (theta) initial for first factor, g (gamma) initial for second factor 
  double precision tccdf,tcder1,tcder2,gccdf,gcder1,gcder2,gpdf,gpder1,gpder1u

  ! npar=2*d
  allocate ( yvec(d), pmf(d), tder1(d), tder2(d), fval1(npar), integl1(npar), grd(npar) )
  allocate ( th(d), gam(d), gder1(d),gder2(d),tgder(d) )
  allocate ( fval2(npar,npar), integl2(npar,npar), hss(npar,npar) )
  allocate ( cond(nq,nq,ncat,d), tcond1(nq,nq,ncat,d), tcond2(nq,nq,ncat,d) )
  allocate ( gcond1(nq,nq,ncat,d), gcond2(nq,nq,ncat,d), tgcond2(nq,nq,ncat,d) )
  ! cond for cond cdf and pmf, cond1 for deriv wrt theta, cond2 for 2nd deriv
   
  th(1:d)=param(1:d); gam(1:d)=param((d+1):npar)
  cond=0.d0; tcond1=0.d0; tcond2=0.d0
  gcond1=0.d0; gcond2=0.d0; tgcond2=0.d0
  cond(:,:,ncat,:)=1.d0  ! condcdf is 1 for last category ncat
  ! get \p copcond(ucuts[k]|v,th[j]) / \p th[j] , k=1,...,ncat ; v=latent
  ! below might not be most efficient
  do iq=1,nq ! loop over quadrature points
  do iq2=1,nq ! loop over quadrature points
    do j=1,d ! loop over items
      do k=1,ncat-1 ! loop of categories
        call irgum1derivs(xl(iq),ucuts(k,j),th(j),tccdf,tcder1,tcder2)
        !if(iq==1 .and. iq2==1 .and. k==1) print "(4f11.6)", th(j),tccdf,tcder1,tcder2
        call irgum2derivs(xl(iq2),tccdf,gam(j),gccdf,gcder1,gcder2,gpdf,gpder1,gpder1u)
        !if(iq==1 .and. iq2==1 .and. k==1) print "(7f11.6)", gam(j),gccdf,gcder1,gcder2,gpdf,gpder1,gpder1u
        cond(iq,iq2,k,j)=gccdf
        gcond1(iq,iq2,k,j)=gcder1
        gcond2(iq,iq2,k,j)=gcder2
        tcond1(iq,iq2,k,j)=gpdf*tcder1
        tcond2(iq,iq2,k,j)=gpdf*tcder2+gpder1u*tcder1*tcder1
        tgcond2(iq,iq2,k,j)=gpder1*tcder1
      end do
    end do
  end do
  end do
  ! differencing
  do iq=1,nq 
  do iq2=1,nq 
    do j=1,d
      do k=ncat,2,-1
        cond(iq,iq2,k,j)=cond(iq,iq2,k,j)-cond(iq,iq2,k-1,j)
        tcond1(iq,iq2,k,j)=tcond1(iq,iq2,k,j)-tcond1(iq,iq2,k-1,j)
        tcond2(iq,iq2,k,j)=tcond2(iq,iq2,k,j)-tcond2(iq,iq2,k-1,j)
        gcond1(iq,iq2,k,j)=gcond1(iq,iq2,k,j)-gcond1(iq,iq2,k-1,j)
        gcond2(iq,iq2,k,j)=gcond2(iq,iq2,k,j)-gcond2(iq,iq2,k-1,j)
        tgcond2(iq,iq2,k,j)=tgcond2(iq,iq2,k,j)-tgcond2(iq,iq2,k-1,j)
      end do
    end do
  end do
  end do
  !print "(8f10.6)", cond(1,1,:,d)
  !print "(8f10.6)", tcond1(1,1,:,d)
  !print "(8f10.6)", gcond1(1,1,:,d)
  ! looping through data
  nllk=0.d0; grad=0.d0; hess=0.d0; ! initialize
  do i=1,n ! loop over rows of data set
    yvec=ydata(i,:)+1
    !if(i==1) print "(8i3)", yvec
    integl=0.d0; integl1=0.d0; integl2=0.d0; ! initialize for integrals
    do iq=1,nq 
    do iq2=1,nq 
      ! update integrand, and derivs wrt copula parameters 
      do j=1,d
        pmf(j)=cond(iq,iq2,yvec(j),j)
        tder1(j)=tcond1(iq,iq2,yvec(j),j)  ! deriv wrt th
        tder2(j)=tcond2(iq,iq2,yvec(j),j)
        gder1(j)=gcond1(iq,iq2,yvec(j),j)  ! deriv wrt gam
        gder2(j)=gcond2(iq,iq2,yvec(j),j)
        tgder(j)=tgcond2(iq,iq2,yvec(j),j)
      end do
      fval=product(pmf)
      ! fval1: vector of partial deriv wrt theta[j], j=1,...,d
      ! fval2: matrix of 2nd partial deriv wrt theta[j], theta[j2], 
      do j=1,d
        fval1(j)=fval*tder1(j)/pmf(j)
        fval1(j+d)=fval*gder1(j)/pmf(j)
        fval2(j,j)=fval*tder2(j)/pmf(j)
        fval2(j+d,j+d)=fval*gder2(j)/pmf(j)
        ! cross th , gam terms for fval2
        fval2(j,j+d)=fval* tgder(j)/pmf(j) 
        fval2(j+d,j) = fval2(j,j+d)
      end do
      do j=1,d
        do j2=1,d
          if(j2 /= j) then
            fval2(j,j2)=fval1(j)*tder1(j2)/pmf(j2)  
            fval2(j+d,j2+d)=fval1(j+d)*gder1(j2)/pmf(j2)
            fval2(j+d,j2)=fval1(j+d)*tder1(j2)/pmf(j2) 
            fval2(j,j2+d)=fval1(j)*gder1(j2)/pmf(j2) 
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
        hss(j,j2)=integl2(j,j2)/integl - grd(j)*grd(j2); 
      end do
    end do
    hess=hess-hss
  end do
  deallocate (yvec, pmf, tder1, tder2, fval1, fval2, integl1, integl2, grd, hss)
  deallocate (cond, tcond1, tcond2, gcond1, gcond2, tgcond2)
  deallocate (th, gam, gder1, gder2, tgder)
  return
  end


! Function with Gumbel condcdf and its deriv wrt theta, order 1 and 2.
! For derivatives, ccdf is function of theta and m=(x^theta+y^theta)^(1/theta)
!  conditioning is u2|u1
! inputs
!   u1,u2 = values in (0,1)
!   theta = scalar in (1,oo)
! outputs
!   ccdf = conditional cdf, 
!   cder1 = \p ccdf/\p theta
!   cder2 = \p^2 ccdf/\p theta^2
subroutine irgum1derivs(u1,u2,theta,ccdf,cder1,cder2)
  implicit none
  double precision u1,u2,theta,ccdf,cder1,cder2
  double precision x,y,tx,ty,xd,yd,s,m,logs,logm,msq,thsq,thcu
  double precision sder1,sder2,mder1,mder2
  double precision lccdf,lcder1,lcder2
  x = -log(u1); y = -log(u2);
  tx = log(x); ty = log(y); 
  xd = x**theta; yd = y**theta; s = xd+yd; m = s**(1.d0/theta);
  logs=log(s);
  thsq=theta*theta; thcu=thsq*theta
  sder1 = xd*tx+yd*ty;
  sder2 = xd*tx*tx+yd*ty*ty;
  mder1 = m*sder1/(s*theta)-m*logs/thsq;
  mder2 = -mder1*logs/thsq-2.d0*m*sder1/(s*thsq)+2.d0*m*logs/thcu;
  mder2 = mder2+sder2*m/(s*theta)+(mder1/s-m*sder1/(s*s))*sder1/theta;
  logm=log(m); msq=m*m
  lccdf = x-m+(1.d0-theta)*(logm-tx)   !condition on x
  lcder1 = -mder1+(1.d0-theta)*mder1/m-logm+tx;  !condition on x
  lcder2 = -mder2-2.d0*mder1/m+(1.d0-theta)*(mder2/m-mder1*mder1/msq);
  ccdf = exp(lccdf)
  cder1 = ccdf*lcder1
  cder2 = ccdf*(lcder1*lcder1+lcder2)
  return
  end

! Function with Gumbel condcdf and its deriv wrt theta, order 1 and 2;
!   also pdf, pdf1th= \p pdf/\p theta, pdf1v= \p pdf/ \p u2 (second arg)
! For derivatives, ccdf is function of theta and m=(x^theta+y^theta)^(1/theta)
! inputs
!   u1,u2 = values in (0,1)
!   theta = scalar in (1,oo)
! outputs
!   ccdf = conditional cdf, 
!   cder1 = \p ccdf/\p theta
!   cder2 = \p^2 ccdf/\p theta^2
!   pdf = Gumbel pdf
!   pdf1th= \p pdf/\p theta, 
!   pdf1v= \p pdf/ \p u2
subroutine irgum2derivs(u1,u2,theta,ccdf,cder1,cder2,pdf,pdf1th,pdf1v)
  implicit none
  double precision u1,u2,theta,ccdf,cder1,cder2,pdf,pdf1th,pdf1v
  double precision x,y,tx,ty,xd,yd,s,m,logs,logm,msq,thsq,thcu
  double precision sder1,sder2,mder1,mder2
  double precision lccdf,lcder1,lcder2
  double precision den,den2,mu,lpdf,lder1,lder1v
  x = -log(u1); y = -log(u2);
  tx = log(x); ty = log(y); 
  xd = x**theta; yd = y**theta; s = xd+yd; m = s**(1.d0/theta);
  logs=log(s);
  thsq=theta*theta; thcu=thsq*theta
  sder1 = xd*tx+yd*ty;
  sder2 = xd*tx*tx+yd*ty*ty;
  mder1 = m*sder1/(s*theta)-m*logs/thsq;
  mder2 = -mder1*logs/thsq-2.d0*m*sder1/(s*thsq)+2.d0*m*logs/thcu;
  mder2 = mder2+sder2*m/(s*theta)+(mder1/s-m*sder1/(s*s))*sder1/theta;
  logm=log(m); msq=m*m
  ! conditional and deriv wrt theta
  lccdf = x-m+(1.d0-theta)*(logm-tx)   !condition on x
  lcder1 = -mder1+(1.d0-theta)*mder1/m-logm+tx;  !condition on x
  lcder2 = -mder2-2.d0*mder1/m+(1.d0-theta)*(mder2/m-mder1*mder1/msq);
  ccdf = exp(lccdf)
  cder1 = ccdf*lcder1
  cder2 = ccdf*(lcder1*lcder1+lcder2)
  ! pdf and deriv wrt theta, u1
  den = m+theta-1.d0; den2=den*den
  lpdf = -m+log(den)+(1.d0-2.d0*theta)*logm+(theta-1.d0)*(tx+ty)+x+y;
  lder1 = -mder1+(mder1+1.d0)/den-2.d0*logm+(1.d0-2.d0*theta)*mder1/m+tx+ty;
  ! replace with x,y interchanged, rename lder1u to lder1v
  !mu = -m*xd/(u1*s*x);
  !lder1u = -mu+mu/den+(1-2*theta)*mu/m-(theta-1)/(u1*x)-1/u1;
  mu = -m*yd/(u2*s*y);
  lder1v = -mu+mu/den+(1.d0-2.d0*theta)*mu/m-(theta-1.d0)/(u2*y)-1.d0/u2;
  pdf=exp(lpdf); pdf1th=pdf*lder1; pdf1v=pdf*lder1v
  return
  end

