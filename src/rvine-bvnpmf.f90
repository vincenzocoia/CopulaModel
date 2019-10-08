! R-vine model with ordinal response and bivariate Gaussian pair-copulas

!program mainex
!  implicit none
!  integer d,iprint,ncateg,ii,nn,d2
!  double precision, dimension(:), allocatable :: u1vec,u2vec,parvec
!  double precision, dimension(:), allocatable :: z1vec,z2vec
!  double precision, dimension(:,:), allocatable :: parmat,zzmat
!  integer, dimension(:,:), allocatable :: A,M
!  integer ncl
!  double precision lk,nllk
!  d=4; iprint=1; d2=(d*(d-1))/2
!  allocate ( u1vec(d),u2vec(d), z1vec(d),z2vec(d), parmat(d,d), A(d,d), M(d,d) )
!  allocate ( parvec(d2), zzmat(1,2*d) )
!  u1vec = (/.1d0,.3d0,.4d0,.5d0/)
!  u2vec = (/.4d0,.6d0,.7d0,.8d0/)
!  z1vec = (/-1.2815516d0, -0.5244005d0, -0.2533471d0,  0.0d0/)
!  z2vec = (/-0.2533471d0,  0.2533471d0,  0.5244005d0,  0.8416212d0/)
!  parmat=0.d0
!  parmat(1,2:d)=0.5d0
!  parmat(2,3:d)=0.1d0
!  A=0 
!  A(1,1:4)=(/1,1,2,3/)
!  A(2,2:4)=(/2,1,2/)
!  A(3,3:4)=(/3,1/)
!  A(4,4)=4
!  M=A
!  M(2,3:4)=(/2,3/)
!  M(3,4)=3
!  call rvinepmf(d,parmat,u1vec,u2vec,A,M,iprint,lk)
!  print *,lk
!  call rvinebvnpmf(d,parmat,z1vec,z2vec,A,M,iprint,lk)
!  print *,lk
!  parvec=(/0.5d0,0.5d0,0.5d0,0.1d0,0.1d0,0.d0/)
!  ncl=1
!  zzmat(1,1:4)=z1vec; zzmat(1,5:8)=z2vec
!  call rvinediscbvnnllk(ncl,d,parvec,zzmat,A,M,nllk)
!  print *,nllk
!  deallocate ( u1vec,u2vec, z1vec, z2vec, parmat, A, M, parvec, zzmat )
!  stop
!  end

!program mainex2
!  implicit none
!  integer d,iprint,ncateg,ii,nn,d2
!  double precision, dimension(:), allocatable :: u1vec,u2vec,parvec,pr
!  double precision, dimension(:,:), allocatable :: parmat,ucuts
!  integer, dimension(:), allocatable :: jj
!  integer, dimension(:,:), allocatable :: A,M
!  double precision lk,kl
!  d=4; iprint=1; d2=(d*(d-1))/2
!  allocate ( u1vec(d),u2vec(d), parmat(d,d), A(d,d), M(d,d), jj(d) )
!  allocate ( parvec(d2) )
!  u1vec = (/.1d0,.3d0,.4d0,.5d0/)
!  u2vec = (/.4d0,.6d0,.7d0,.8d0/)
!  parmat=0.d0
!  parmat(1,2:d)=0.5d0
!  parmat(2,3:d)=0.1d0
!  A=0 
!  A(1,1:4)=(/1,1,2,3/)
!  A(2,2:4)=(/2,1,2/)
!  A(3,3:4)=(/3,1/)
!  A(4,4)=4
!  M=A
!  M(2,3:4)=(/2,3/)
!  M(3,4)=3
!  call rvinepmf(d,parmat,u1vec,u2vec,A,M,iprint,lk)
!  print *,lk
!
  ! add test for d2v
!  ncateg=3
!  nn=ncateg**d
!  do ii = 1,nn
!    call d2vplus1(d,ncateg,ii-1,jj)
!    !print *,ii, jj
!  end do
  ! add test for rvineklfn
!  allocate ( pr(nn), ucuts(ncateg+1,d) )
!  ucuts(1,:)=0.0001d0; ucuts(d,:)=0.9999d0
!  ucuts(2,:) = (/.4d0,.5d0,.4d0,.3d0/)
!  ucuts(3,:) = (/.7d0,.8d0,.6d0,.6d0/)
!  parvec=(/0.5d0,0.5d0,0.5d0,0.1d0,0.1d0,0.d0/)
!  pr=1.d0/nn
!  call rvineklfn(d,ncateg,parvec,ucuts,A,M,pr,kl)
!  print *,kl
!  deallocate ( u1vec,u2vec, parmat, A, M, jj, parvec,pr, ucuts )
!  stop
!  end

! R-vine full log-likelihood for ordinal probit
! inputs
!   nn = nrep*ncl, ncl=#clusters, 
!   d = nrep = #repeated measures per cluster
!   npred = number of predictors
!   ncateg = number of categories
!   param = vector of parameters : univariate plus R-vine partial correlations
!   A = dxd vine array with 1:d on diagonal
!   M = dxd max array computed from A
!   xmat = nn*npred matrix of covariates
!   yvec = integer-valued vector (nnx1)
! output
!   nllk = negative log-likelihood
subroutine rvineOrdwPrednllk(nn,d,npred,ncateg,param,A,M,xmat,yvec,nllk)
  implicit none
  integer nn,d,npred,ncateg
  integer A(d,d),M(d,d),yvec(nn)
  double precision param(ncateg-1+npred+ (d*(d-1))/2),xmat(nn,npred),nllk
  integer ncl,np,ii,j,ell,i,jp,iprint
  double precision, dimension(:), allocatable :: b0,bvec,rhvec,low,upp
! b0cut is vector of cutpoints
! bvec is vector of regression coefficients
! or should input be  uu1, uu2 corners of rectangles?
  double precision, dimension(:,:), allocatable :: parmat
  double precision tem,pr
  
  ncl=nn/d  ! number of clusters
  np=ncateg-1+npred+ (d*(d-1))/2
  !print *, ncl,nn,d,npred,ncateg
  allocate ( b0(ncateg+1),bvec(npred),rhvec((d*(d-1))/2),parmat(d,d) )
  allocate ( low(d), upp(d) )
  b0(1)=-6.d0; b0(ncateg+1)=6.d0
  b0(2:ncateg)=param(1:(ncateg-1))
  do ii = 2,ncateg+1
    if(b0(ii)-b0(ii-1)<=0.d0) then
      nllk=1.d10
      return
    endif
  end do
  !  intercept lower bound is b0[y]  b0[1]=-6 say 
  !  intercept upper bound is b0[y+1]  b0[ncateg+1]=6
  bvec=param(ncateg:(ncateg+npred-1))
  rhvec=param((ncateg+npred):np)
  !print *, b0
  !print *, bvec
  !print *, rhvec
  if( minval(rhvec)<=-1.d0 .or. maxval(rhvec)>=1.d0 ) then
    nllk=1.d10
    return
  endif
  parmat=0.d0
  ii=0
  do ell = 1,(d-1)
    do j = (ell+1),d
      ii=ii+1; parmat(ell,j)=rhvec(ii) 
    end do
  end do
  ! negative log-likelihood: in f90, not necessary to store zzdat
  nllk=0.d0; iprint=0
  do i=1,ncl
    ii=(i-1)*d
    do j=1,d
      tem=0.d0
      do jp=1,npred
        tem=tem+bvec(jp)*xmat(ii+j,jp)
      end do
      low(j)=tem+b0(yvec(ii+j))
      upp(j)=tem+b0(1+yvec(ii+j))
    end do
    call rvinebvnpmf(d,parmat,low,upp,A,M,iprint,pr)
    if(pr<=0) pr=1.d-15
    nllk=nllk-log(pr)
  end do
  deallocate ( b0,bvec,rhvec,parmat, low,upp )
  return
  end

! R-vine nllk for multivariate discrete, with univariate models estimated 
! inputs
!   ncl = number of clusters, 
!   d = nrep = #repeated measures per cluster
!   ncateg = number of categories
!   parvec = vector of R-vine partial correlation parameters of bivariate normal
!   zzdat = ncl x (2d) matrix with corners of rectangle, N(0,1) scale
!   A = dxd vine array with 1:d on diagonal
!   M = dxd max array computed from A
! output
!   nllk = negative log-likelihood
subroutine rvinediscbvnnllk(ncl,d,parvec,zzdat,A,M,nllk)
  implicit none
  integer ncl,d,iprint
  integer A(d,d),M(d,d)
  double precision parvec((d*(d-1))/2),zzdat(ncl,2*d),nllk
  integer ii,j,ell,i
  integer, dimension(:), allocatable :: jj1,jj2
  double precision, dimension(:,:), allocatable :: parmat
  double precision pr
  
  if(maxval(parvec)>=1. .or. minval(parvec)<= -1.) then 
    nllk=1.d10; return
  endif 
  allocate ( parmat(d,d), jj1(d), jj2(d) )
  iprint=0
  parmat=0.d0
  ii=0
  do ell = 1,(d-1)
    do j = (ell+1),d
      ii=ii+1; parmat(ell,j)=parvec(ii) 
    end do
  end do
  nllk=0.d0
  jj1=(/ (j, j=1,d) /); jj2=(/ (j, j=d+1,2*d) /)
  do i = 1,ncl
    call rvinebvnpmf(d,parmat,zzdat(i,jj1),zzdat(i,jj2),A,M,iprint,pr)
    if(pr<=0) pr=1.d-15
    nllk=nllk-log(pr)
  end do
  deallocate ( parmat, jj1, jj2 )
  return
  end


! R-vine probability for multivariate discrete 
!   d = #variables or #repeated measures per cluster
!   parmat = dxd matrix with parameters in same position as vine array A
!   A = dxd vine array with 1:d on diagonal
!   M = dxd max array computed from A
!   z1vec = lower vector of hyperrectangle
!   z2vec = upper vector of hyperrectangle on N(0,1) scale
!   iprint = print flag for intermediate steps in the probability calculations
! output
!   nllk = joint pmf or likelihood 
subroutine rvinebvnpmf(d,parmat,z1vec,z2vec,A,M,iprint,lk)
  implicit none
  integer d,iprint
  integer A(d,d),M(d,d)
  double precision parmat(d,d),z1vec(d),z2vec(d),lk
  integer d1,j,ell,ell1,n0
  double precision, dimension(:), allocatable :: upr,Fp,Fm,Vp,Vbp,Vm,Vbm
  double precision, dimension(:), allocatable :: vbpr,vpr,ovbpr,ovpr,bpp,bpm,bmp,bmm
  double precision, dimension(:,:), allocatable :: mult
  double precision tt,sp,sm
  double precision pbvncop,pnorms ! pnorms in pnorm.f90
  
  allocate ( upr(d), Fp(d), Fm(d), Vp(d), Vbp(d), Vm(d), Vbm(d))
  allocate ( vbpr(d), vpr(d), ovbpr(d), ovpr(d), bpp(d), bpm(d), bmp(d), bmm(d))
  allocate ( mult(d,d) )  ! for f_{j-l..j}, j=1,...,d, l=0,...,d-1

  d1=d-1
  upr=0.d0; Fp=0.d0; Fm=0.d0; Vp=0.d0; Vbp=0.d0; Vm=0.d0; Vbm=0.d0;
  vbpr=0.d0; vpr=0.d0;  ovbpr=0.d0; ovpr=0.d0;  
  bpp=0.d0; bpm=0.d0; bmp=0.d0; bmm=0.d0;
  mult=0.d0
  lk=0.d0
  ! univariate probabilities
  Fp=z2vec; Fm=z1vec; 
  do j=1,d
    upr(j)=pnorms(Fp(j))-pnorms(Fm(j)) 
  end do
  !if(iprint==1) then
  !  print '(10f8.4)',z2vec; print '(10f8.4)',z1vec; print '(10f8.4)',upr; 
  !endif
  mult(1,:)=upr
  ! tree 1
  !pcop=match.fun(pcopnames[1])
  n0=1
  do j = 2,d
    call pbnorm(n0,Fp(A(1,j)),Fp(j),parmat(1,j),bpp(j))
    call pbnorm(n0,Fp(A(1,j)),Fm(j),parmat(1,j),bpm(j))
    call pbnorm(n0,Fm(A(1,j)),Fp(j),parmat(1,j),bmp(j))
    call pbnorm(n0,Fm(A(1,j)),Fm(j),parmat(1,j),bmm(j))
    !bpp(j)=pbnorm(Fp(A(1,j)),Fp(j),parmat(1,j))
    !bpm(j)=pbnorm(Fp(A(1,j)),Fm(j),parmat(1,j))
    !bmp(j)=pbnorm(Fm(A(1,j)),Fp(j),parmat(1,j))
    !bmm(j)=pbnorm(Fm(A(1,j)),Fm(j),parmat(1,j))
  end do
  do j = 2,d
    mult(2,j)=bpp(j)-bpm(j)-bmp(j)+bmm(j) ! f_{j-1,j}
  end do
  !print *, mult
  lk=mult(2,2)   ! f_{1,2}
  !if(iprint==1) print *, 1,lk
  if(iprint==1) call dblepr("lk1",3,lk,1)
  ! j=d not used for vb,  j=2 not used for v?
  do j = 2,d 
    Vbp(j)=(bpp(j)-bpm(j))/upr(j) ! F_{a_{1j}|j}
    Vbm(j)=(bmp(j)-bmm(j))/upr(j)
    vbpr(j)=Vbp(j)-Vbm(j)  ! f_{a_{1j}|j}
    Vp(j)=(bpp(j)-bmp(j))/upr(A(1,j)) ! F_{j|a_{1j}}
    Vm(j)=(bpm(j)-bmm(j))/upr(A(1,j))
    vpr(j)=Vp(j)-Vm(j)    ! f_{j|a_{1j}}
  end do
  ! tree 2 and on
  ! main loop
  do ell = 2,d1 ! tree level
    ell1=ell+1
    do j = ell1,d ! Vb in  first arg, V in second arg
      !sp=ifelse(A(ell,j)<M(ell,j),Vbp(M(ell,j)),Vp(M(ell,j)))
      !sm=ifelse(A(ell,j)<M(ell,j),Vbm(M(ell,j)),Vm(M(ell,j)))
      if (A(ell,j)<M(ell,j)) then
        sp=Vbp(M(ell,j)); sm=Vbm(M(ell,j))
      else
        sp=Vp(M(ell,j)); sm=Vm(M(ell,j))
      endif
      bpp(j)=pbvncop(sp,Vp(j),parmat(ell,j))  
      ! C_{a_{ell,j},j;}
      bpm(j)=pbvncop(sp,Vm(j),parmat(ell,j))
      bmp(j)=pbvncop(sm,Vp(j),parmat(ell,j))
      bmm(j)=pbvncop(sm,Vm(j),parmat(ell,j))
    end do
    ovbpr=vbpr; ovpr=vpr
    do j = ell1,d 
      Vbp(j)=(bpp(j)-bpm(j))/ovpr(j) ! F_{a_{ell,j}|j,S} den prev vpr
      Vbm(j)=(bmp(j)-bmm(j))/ovpr(j)
      vbpr(j)=Vbp(j)-Vbm(j)  ! f_{a_{ell,j}|j,S}
      !tt=ifelse(A(ell,j)<M(ell,j),ovbpr(M(ell,j)),ovpr(M(ell,j)))
      if (A(ell,j)<M(ell,j)) then
        tt=ovbpr(M(ell,j))
      else
        tt=ovpr(M(ell,j))
      endif
      Vp(j)=(bpp(j)-bmp(j))/tt ! F_{j|a_{ell,j},S}
      Vm(j)=(bpm(j)-bmm(j))/tt
      vpr(j)=Vp(j)-Vm(j)    ! f_{j|a_{ell,j},S}
    end do
    do j = ell1,d
      mult(ell1,j)=mult(ell,M(ell,j))*vpr(j)
    end do
    lk=mult(ell1,ell1) ! ell1=d at the end
    !if(iprint==1) print *, ell,lk
     if(iprint==1) call dblepr("lk2+",4,lk,1)
  end do
  deallocate ( upr, Fp, Fm, Vp, Vbp, Vm, Vbm)
  deallocate ( vbpr, vpr, ovbpr, ovpr, bpp, bpm, bmp, bmm)
  deallocate ( mult )  
  return 
  end

! pbvncop is cdf of bivariate normal copula
! inputs
!   u,v = values in (0.1)
!   rho = correlation parameter in (-1,1)
! output
!   bivariate copula cdf
double precision function pbvncop(u,v,rho)
  implicit none
  integer n0
  double precision u,v,rho,x,y,bcdf
  double precision qnorms ! function in C, replace by version in qtnorm.f90
  if(rho==0.d0) then
    pbvncop=u*v
    return
  endif
  if(u>1-1.d-9) u=1-1.d-9
  if(v>1-1.d-9) v=1-1.d-9
  if(u<1.d-9) u=1.d-9
  if(v<1.d-9) v=1.d-9
  x=qnorms(u); y=qnorms(v)
  n0=1
  call pbnorm(n0,x,y,rho,bcdf)
  pbvncop=bcdf
  !print *,u,v,rho,pbvncop
  return
  end

! decimal to vector
! inputs
!   d = dimension, 
!   ncat = #categories, 
!   ii = integer between 1 and ncat^d 
! output
!   jj = vector of length d with elements in {1,2,...,ncat}
subroutine d2vplus1(d,ncat,ii,jj)
  implicit none
  integer d,ncat,ii,jj(d)
  integer i,tem
  tem=ii
  do i=d,1,-1
    jj(i)=mod(tem,ncat)
    tem=tem/ncat
  end do
  jj=jj+1
  ! tem%ncat=mod(tem,ncat) is the remainder from dividing ncat into tem
  !    this is tem%%ncat in R
  ! tem/ncat is integer division rounded down
  !    this is floor(tem/ncat) in R
  return
  end


! KL divergence for discrete R-vine model
! inputs
!   d = nrep = #repeated measures per cluster
!   ncateg = number of categories
!   parvec = vector of length d*(d-1)/2 for R-vine partial correlations 
!   ucuts = (ncateg+1) x d matrix of cutpts in U(0,1) scale for ordinal response
!   A = dxd vine array with 1:d on diagonal
!   M = dxd max array computed from A
!   pr = vector of joint pmf values, e.g., computed with pmfmordprobit in R
! output
!   kl = KL divergence of pr vector and R-vine pmf with Gaussian pair-copulas
subroutine rvineklfn(d,ncateg,parvec,ucuts,A,M,pr,kl)
  implicit none
  integer d,ncateg,A(d,d),M(d,d)
  double precision parvec((d*(d-1))/2),pr(ncateg**d),kl
  double precision ucuts(ncateg+1,d)
  double precision, dimension(:), allocatable :: prv,u1vec,u2vec
  double precision, dimension(:,:), allocatable :: parmat
  integer, dimension(:), allocatable :: jj
  integer nn,d2,ii,ell,j,i,iprint
  !double precision ss
  if(maxval(parvec)>=1. .or. minval(parvec)<= -1.) then 
    kl=1.d10; return
  endif 
  !  d=ncol(ucuts); ncateg=nrow(ucuts)-1
  d2=(d*(d-1))/2; nn=ncateg**d
  allocate ( prv(nn), u1vec(d), u2vec(d), parmat(d,d), jj(d) )
  parmat=0.d0
  iprint=0
  ii=0
  do ell = 1,(d-1)
    do j = (ell+1),d
      ii=ii+1; parmat(ell,j)=parvec(ii) 
    end do
  end do
  !print *, parvec
  kl=0.d0; !ss=0.d0
  do i = 1,nn
    call d2vplus1(d,ncateg,i-1,jj)
    do j = 1,d 
      u1vec(j)=ucuts(jj(j),j); u2vec(j)=ucuts(jj(j)+1,j) 
    end do
    call rvinepmf(d,parmat,u1vec,u2vec,A,M,iprint,prv(i))
   ! ? replace later with rvinebvnpmf(d,parmat,z1vec,z2vec,A,M,iprint,lk)
    !print *, i,pr(i),prv(i)
    if(prv(i)<=0) prv(i)=1.d-10
    kl=kl+pr(i)*log(pr(i)/prv(i))
    !ss=ss+prv(i)
  end do
  deallocate ( prv, u1vec, u2vec, parmat, jj )
  !print *, "ss=", ss
  return
  end

! R-vine probability for multivariate discrete  with Gaussian pair-copulas
!   d = nrep = #repeated measures per cluster
!   parmat is dxd matrix with parameters in same position as vine array A
!   A = dxd vine array with 1:d on diagonal
!   M = dxd max array computed from A
!   u1vec = lower vector of hyperrectangle
!   u2vec = upper vector of hyperrectangle on U(0,1) scale
!   iprint = print flag for intermediate steps in the probability calculations
! output
!    lk = pmf of the R-vine
subroutine rvinepmf(d,parmat,u1vec,u2vec,A,M,iprint,lk)
  implicit none
  integer d,iprint
  integer A(d,d),M(d,d)
  double precision parmat(d,d),u1vec(d),u2vec(d),lk
  integer d1,j,ell,ell1
  double precision, dimension(:), allocatable :: upr,Fp,Fm,Vp,Vbp,Vm,Vbm
  double precision, dimension(:), allocatable :: vbpr,vpr,ovbpr,ovpr,bpp,bpm,bmp,bmm
  double precision, dimension(:,:), allocatable :: mult
  double precision tt,sp,sm
  double precision pbvncop
  
  allocate ( upr(d), Fp(d), Fm(d), Vp(d), Vbp(d), Vm(d), Vbm(d))
  allocate ( vbpr(d), vpr(d), ovbpr(d), ovpr(d), bpp(d), bpm(d), bmp(d), bmm(d))
  allocate ( mult(d,d) )  ! for f_{j-l..j}, j=1,...,d, l=0,...,d-1

  !d=length(u1vec)
  d1=d-1
  upr=0.d0; Fp=0.d0; Fm=0.d0; Vp=0.d0; Vbp=0.d0; Vm=0.d0; Vbm=0.d0;
  vbpr=0.d0; vpr=0.d0;  ovbpr=0.d0; ovpr=0.d0;  
  bpp=0.d0; bpm=0.d0; bmp=0.d0; bmm=0.d0;
  mult=0.d0
  lk=0.d0
  ! univariate probabilities
  Fp=u2vec; Fm=u1vec; upr=Fp-Fm 
  !if(iprint==1) then
  !  print '(10f8.4)',u2vec; print '(10f8.4)',u1vec; print '(10f8.4)',upr; 
  !endif
  mult(1,:)=upr
  ! tree 1
  !pcop=match.fun(pcopnames[1])
  do j = 2,d
    bpp(j)=pbvncop(Fp(A(1,j)),Fp(j),parmat(1,j))
    bpm(j)=pbvncop(Fp(A(1,j)),Fm(j),parmat(1,j))
    bmp(j)=pbvncop(Fm(A(1,j)),Fp(j),parmat(1,j))
    bmm(j)=pbvncop(Fm(A(1,j)),Fm(j),parmat(1,j))
  end do
  do j = 2,d
    mult(2,j)=bpp(j)-bpm(j)-bmp(j)+bmm(j) ! f_{j-1,j}
  end do
  !print *, mult
  lk=mult(2,2)   ! f_{1,2}
  !if(iprint==1) print *, 1,lk
  if(iprint==1) call dblepr("lk1",3,lk,1)
  ! j=d not used for vb,  j=2 not used for v?
  do j = 2,d 
    Vbp(j)=(bpp(j)-bpm(j))/upr(j) ! F_{a_{1j}|j}
    Vbm(j)=(bmp(j)-bmm(j))/upr(j)
    vbpr(j)=Vbp(j)-Vbm(j)  ! f_{a_{1j}|j}
    Vp(j)=(bpp(j)-bmp(j))/upr(A(1,j)) ! F_{j|a_{1j}}
    Vm(j)=(bpm(j)-bmm(j))/upr(A(1,j))
    vpr(j)=Vp(j)-Vm(j)    ! f_{j|a_{1j}}
  end do
  ! tree 2 and on
  ! main loop
  do ell = 2,d1 ! tree level
    ell1=ell+1
    do j = ell1,d ! Vb in  first arg, V in second arg
      !sp=ifelse(A(ell,j)<M(ell,j),Vbp(M(ell,j)),Vp(M(ell,j)))
      !sm=ifelse(A(ell,j)<M(ell,j),Vbm(M(ell,j)),Vm(M(ell,j)))
      if (A(ell,j)<M(ell,j)) then
        sp=Vbp(M(ell,j)); sm=Vbm(M(ell,j))
      else
        sp=Vp(M(ell,j)); sm=Vm(M(ell,j))
      endif
      bpp(j)=pbvncop(sp,Vp(j),parmat(ell,j))  
      ! C_{a_{ell,j},j;}
      bpm(j)=pbvncop(sp,Vm(j),parmat(ell,j))
      bmp(j)=pbvncop(sm,Vp(j),parmat(ell,j))
      bmm(j)=pbvncop(sm,Vm(j),parmat(ell,j))
    end do
    ovbpr=vbpr; ovpr=vpr
    do j = ell1,d 
      Vbp(j)=(bpp(j)-bpm(j))/ovpr(j) ! F_{a_{ell,j}|j,S} den prev vpr
      Vbm(j)=(bmp(j)-bmm(j))/ovpr(j)
      vbpr(j)=Vbp(j)-Vbm(j)  ! f_{a_{ell,j}|j,S}
      !tt=ifelse(A(ell,j)<M(ell,j),ovbpr(M(ell,j)),ovpr(M(ell,j)))
      if (A(ell,j)<M(ell,j)) then
        tt=ovbpr(M(ell,j))
      else
        tt=ovpr(M(ell,j))
      endif
      Vp(j)=(bpp(j)-bmp(j))/tt ! F_{j|a_{ell,j},S}
      Vm(j)=(bpm(j)-bmm(j))/tt
      vpr(j)=Vp(j)-Vm(j)    ! f_{j|a_{ell,j},S}
    end do
    do j = ell1,d
      mult(ell1,j)=mult(ell,M(ell,j))*vpr(j)
    end do
    lk=mult(ell1,ell1) ! ell1=d at the end
    !if(iprint==1) print *, ell,lk
    if(iprint==1) call dblepr("lk2+",4,lk,1)
  end do
  !if(is.na(lk)) ! how to handle this??
  !{ lk=0.; 
  !  if(iprint) { print(parmat); print(mult); }
  !}
  !lk
  deallocate ( upr, Fp, Fm, Vp, Vbp, Vm, Vbm)
  deallocate ( vbpr, vpr, ovbpr, ovpr, bpp, bpm, bmp, bmm)
  deallocate ( mult )  
  return 
  end

