! sample main program to test subroutines
!  gfortran -o trifact mvttrifact.f90 factorderiv0.f90 lsolve.f90
!  factorderiv0.f90 is mvtfactwithgrad.f90 with intpr replaced by print (twice)

!program mainex
!  implicit none
!  integer d,mgrp,msbgrp,p,i,k,j,n
!  integer, dimension(:), allocatable :: grsize,sbgrsize
!  double precision, dimension(:), allocatable :: rh
!  double precision, dimension(:,:), allocatable :: fctmat,fctinv,grad,Robs
!  double precision, dimension(:,:), allocatable :: rhmat,amat,lgdtder,rinvder
!  double precision, dimension(:,:,:), allocatable :: jarr
!  double precision fctdet,fctldet,nllk

!  n=50
!  p=3
!  read *, mgrp,msbgrp
!  allocate (grsize(mgrp),sbgrsize(msbgrp))
!  read *, grsize
!  read *, sbgrsize
!  d=sum(grsize)  ! also sum(sbgrsize)
!  allocate (rh(d*p),fctmat(d,d),fctinv(d,d),rhmat(d,p),amat(d,p))
!  allocate (lgdtder(d,p), jarr(p,p,d), rinvder(d,d), grad(d,p), Robs(d,d))
!  do i=1,d*3
!    rh(d*3+1-i)=i/(d*p+1.d0)
!  end do
  !print "(10f8.5)", rh
!  do i=1,d
!    do j=1,d
!      Robs(i,j)=0.6d0**(abs(i-j))
!    end do
!  end do
!  print '(8f10.5)', Robs(1,:)
!  call trifct(d,mgrp,msbgrp,grsize,sbgrsize,rh(1:d),rh((d+1):(2*d)), & 
!      rh((2*d+1):(3*d)),fctmat,fctinv,fctdet,fctldet)
   !do i=1,d
   !  print '(8f10.5)',fctmat(i,:)
   !end do
!   print '(8f10.5)', fctdet
   !do i=1,d
   !  print '(8f10.5)',fctinv(i,:)
   !end do
   !other tests
!  do k=1,p
!    rhmat(:,k)=rh(((k-1)*d+1):(k*d))
!  end do
!  call pcor2load(d,p,rhmat,amat)
  !do i=1,d
  !  print '(8f10.5)',amat(i,:)
  !end do

!  call derivlogdetdloadtrif(d,mgrp,msbgrp,amat,fctinv,grsize,sbgrsize,lgdtder)
  !do i=1,d
  !  print '(5f12.8)',lgdtder(i,:)
  !end do
  !call jacobload2pcor(d,p,amat,rhmat,jarr)
  !do i=1,d
  !  print '(9f11.7)',jarr(:,:,i)
  !end do
  !call derivRinvdloadtrif(d,mgrp,msbgrp,amat,fctinv,4,1,grsize,sbgrsize,rinvder)
  !print *, "4,1"
  !do i=1,d
  !  print '(8f10.5)',rinvder(i,:)
  !end do
  !call derivRinvdloadtrif(d,mgrp,msbgrp,amat,fctinv,3,2,grsize,sbgrsize,rinvder)
  !print *, "3,2"
  !do i=1,d
  !  print '(8f10.5)',rinvder(i,:)
  !end do
  !call derivRinvdloadtrif(d,mgrp,msbgrp,amat,fctinv,4,3,grsize,sbgrsize,rinvder)
  !print *, "4,3"  ! Ok
  !do i=1,d
  !  print '(8f10.5)',rinvder(i,:)
  !end do
!  print *, "deriv of nllk"
!  call trifactnllk(d,mgrp,msbgrp,rhmat,grsize,sbgrsize,Robs,n,nllk,grad)
!  print *,nllk
!  do i=1,d
!    print '(5f13.8)',grad(i,:)
!  end do
!  deallocate (rh,fctmat,fctinv, rhmat,amat, lgdtder, jarr, rinvder,grad,Robs)
!  deallocate (grsize,sbgrsize)
!  stop
!  end

! needed for compiling this file, not for R package
! outer product of two vectors
!subroutine outer(na,nb,avec,bvec,abmat)
!  implicit none
!  integer na,nb,i
!  double precision avec(na),bvec(nb),abmat(na,nb)
!  ! na=size(avec); nb=size(bvec)
!  do i =1,na  
!    abmat(i,:)=avec(i)*bvec
!  end do
!  return
!  end

!============================================================

!computes matrix, its determinant and inverse in a tri-factor normal model
!         rh1,rh2 - *conditional* correlations (1st factor, 2nd );
! inputs  
!   d = #variables
!   mgrp = #groups
!   msbgrp = #subgroups
!   grsize = vector of group sizes, length is mgrp
!   sbgrsize = vector of subgroup sizes, length is msbgrp
!   rh1 = d-vector of correlations with common latent
!   rh2 = d-vector of partial correlations with group latent given common latent
!   rh3 = d-vector of partial correlations with subgroup latent given common 
!          and group latent
! outputs 
!   fctmat = covariance matrix in tri-factor model with parameters rh1,rh2
!   fctinv = inverse of fctmat 
!   fctdet = det(fctmat)   
!   fctldet = log ( det(fctmat) )
subroutine trifct(d,mgrp,msbgrp,grsize,sbgrsize,rh1,rh2,rh3,fctmat,fctinv,fctdet,fctldet)
  implicit none
  integer d,mgrp,msbgrp,jg,i,temi,kg,pfact,j1,j2,col
  integer grsize(mgrp),sbgrsize(msbgrp)
  double precision rh1(d),rh2(d),rh3(d),fctmat(d,d),fctinv(d,d)
  double precision fctdet,fctldet,ss,tol
  integer, dimension(:), allocatable :: nsb,nbeg,nend
  double precision, dimension(:), allocatable :: a2,a3,tem,tem2,psi2
  double precision, dimension(:,:), allocatable :: matp,prmat,tprmat,nmat,soln
  !double precision, dimension(:,:), allocatable :: fctdiag !, invpsi 
  integer, dimension(:), allocatable :: cmgr,cmst2,cmsbgr,cmst3

  tol=1.d-6
  pfact=1+mgrp+msbgrp  ! number of factors 
  ! code below works if 1<pfact<d
  allocate (nsb(mgrp),nbeg(mgrp),nend(mgrp)) ! #subgroups for each group etc
  allocate (a2(d),a3(d),tem(d),tem2(d),psi2(d))
  allocate (matp(d,pfact),prmat(d,pfact),tprmat(pfact,d),nmat(pfact,pfact)) 
  allocate (soln(pfact,d))
  allocate (cmgr(mgrp),cmst2(mgrp),cmsbgr(msbgrp),cmst3(msbgrp))
  !allocate (fctdiag(d,d) )
  temi=0
  do jg=1,mgrp
    cmst2(jg)=temi+1
    temi=temi+grsize(jg)
    cmgr(jg)=temi
  end do
  temi=0
  do kg=1,msbgrp
    cmst3(kg)=temi+1
    temi=temi+sbgrsize(kg)
    cmsbgr(kg)=temi
  end do
  nsb=0
  jg=1 
  do kg=1,msbgrp
    if(cmsbgr(kg)>cmgr(jg)) jg=jg+1
    nsb(jg)=nsb(jg)+1
  end do
  temi=0
  do jg=1,mgrp
    nbeg(jg)=temi+1
    temi=temi+nsb(jg)
    nend(jg)=temi
  end do

  a2=rh2*sqrt(1.d0-rh1**2)
  a3=rh3*sqrt(1.d0-rh1**2)*sqrt(1.d0-rh2**2)
  psi2=1.d0-rh1**2-a2**2-a3**2
  !invpsi=0.d0
  !do i=1,d
  !  invpsi(i,i)=1.d0/psi2(i) ! is this needed?
  !end do
  matp(:,1)=rh1
  ! deduce number of subgroups for each group
  ! next columns for group 1, then subgroups of group 1
  !   then group 2 and subgroups of group 2, etc
  col=1
  do jg=1,mgrp
    tem=0.d0
    tem(cmst2(jg):cmgr(jg))=a2(cmst2(jg):cmgr(jg))
    !matp(:,1+jg)=tem
    col=col+1; matp(:,col)=tem
    do kg=nbeg(jg),nend(jg)
      tem2=0.d0
      tem2(cmst3(kg):cmsbgr(kg))=a3(cmst3(kg):cmsbgr(kg))
      col=col+1; matp(:,col)=tem2
    end do
  end do
  
  fctmat=matmul(matp,transpose(matp))
  do i=1,d
    fctmat(i,i)=fctmat(i,i)+psi2(i)  ! diagonal is now 1
  end do
  do j1=1,pfact
    do j2=1,pfact
      ss=0.d0
      do i=1,d
         ss=ss+matp(i,j1)*matp(i,j2)/psi2(i)
      end do
      nmat(j1,j2)=ss
    end do
    do i=1,d
      tprmat(j1,i)=matp(i,j1)/psi2(i)
      prmat(i,j1)=matp(i,j1)/psi2(i)
    end do
  end do
  do j1=1,pfact
    nmat(j1,j1)=nmat(j1,j1)+1.d0
  end do
  !do i=1,d
  !  print "(10f8.5)", prmat(i,:)
  !end do
  !do i=1,pfact
  !  print "(12f8.5)", nmat(i,:)
  !end do
  ! nmat is pfact x pfact, prmat is pfact x d, soln is pfact x d
  call lsolve(pfact,d,nmat,tprmat,tol,soln,fctdet)
  fctdet=fctdet*product(psi2)
  fctldet=log(fctdet) 
  fctinv=-matmul(prmat,soln)
  do i=1,d
    fctinv(i,i)=fctinv(i,i)+1.d0/psi2(i)  
  end do
  deallocate (nsb,nbeg,nend)
  deallocate (a2,a3,tem,tem2,psi2)
  deallocate (matp,prmat,tprmat,nmat) 
  deallocate (cmgr,cmst2,cmsbgr,cmst3)
  return
  end

! partial derivatives of log det for tri-factor
! inputs  
!   d = #variables
!   mgrp = #groups
!   msbgrp = #subgroups
!   amat = d x 3 matrix of loadings
!   rinv = inverse of correlation matrix
!   grsize = vector of group sizes, length is mgrp
!   sbgrsize = vector of subgroup sizes, length is msbgrp
! output 
!   lgdtder = d x 3 matrix with derivative of log det wrt each loading 
subroutine derivlogdetdloadtrif(d,mgrp,msbgrp,amat,rinv,grsize,sbgrsize,lgdtder)
  implicit none
  integer d,mgrp,msbgrp,ii,jg,temi,kg
  double precision amat(d,3), rinv(d,d), lgdtder(d,3) 
  integer grsize(mgrp),sbgrsize(msbgrp)
  double precision, dimension(:), allocatable :: vec1,vec2
  integer, dimension(:), allocatable :: cmgr,cmst2,cmsbgr,cmst3
  allocate (vec1(d),vec2(d),cmgr(mgrp),cmst2(mgrp))
  allocate (cmsbgr(msbgrp),cmst3(msbgrp))
  temi=0
  do jg=1,mgrp
    cmst2(jg)=temi+1
    temi=temi+grsize(jg)
    cmgr(jg)=temi
  end do
  temi=0
  do kg=1,msbgrp
    cmst3(kg)=temi+1
    temi=temi+sbgrsize(kg)
    cmsbgr(kg)=temi
  end do
  ! common factor in col 1
  do ii=1,d
    vec1=rinv(:,ii) 
    vec2=amat(:,1); vec2(ii)=0.d0
    lgdtder(ii,1)=2.d0*sum(vec1*vec2)
  end do
  !print "(8f10.5)", lgdtder(:,1)
  ! group factor in col 2
  do jg=1,mgrp
      vec1=0.d0
    do ii=cmst2(jg),cmgr(jg)
      vec1(cmst2(jg):cmgr(jg))=rinv(cmst2(jg):cmgr(jg),ii)
      vec2(cmst2(jg):cmgr(jg))=amat(cmst2(jg):cmgr(jg),2)
      vec2(ii)=0.d0
      lgdtder(ii,2)=2.d0*sum(vec1(cmst2(jg):cmgr(jg))*vec2(cmst2(jg):cmgr(jg)))
    end do
  end do
  !print "(8f10.5)", lgdtder(:,2)
  ! subgroup factor in col 3
  do kg=1,msbgrp
      vec1=0.d0
    do ii=cmst3(kg),cmsbgr(kg)
      vec1(cmst3(kg):cmsbgr(kg))=rinv(cmst3(kg):cmsbgr(kg),ii)
      vec2(cmst3(kg):cmsbgr(kg))=amat(cmst3(kg):cmsbgr(kg),3)
      vec2(ii)=0.d0
      lgdtder(ii,3)=2.d0*sum(vec1(cmst3(kg):cmsbgr(kg))*vec2(cmst3(kg):cmsbgr(kg)))
    end do
  end do
  deallocate (vec1,vec2,cmgr,cmst2,cmsbgr,cmst3)
  return 
  end

! derivative of rinv wrt a loading parameter
! inputs  
!   d = #variables
!   mgrp = #groups
!   msbgrp = #subgroups
!   amat = d x 3 matrix of loadings
!   rinv = inverse of correlation matrix
!   rinv = inverse of correlation matrix
!   ii = index of variable. 1<=ii<=d
!   kk = index of factor. 1<=kk<=3
!   grsize = vector of group sizes, length is mgrp
!   sbgrsize = vector of subgroup sizes, length is msbgrp
! output 
!   deriv = dxd matrix with derivative of components of rinv with load[ii,kk] 
subroutine derivRinvdloadtrif(d,mgrp,msbgrp,amat,rinv,ii,kk,grsize,sbgrsize,deriv)
  implicit none
  integer d,mgrp,msbgrp,ii,kk,j1,j2,jg,temi,kg
  double precision amat(d,3), rinv(d,d), deriv(d,d)
  integer grsize(mgrp),sbgrsize(msbgrp)
  double precision, dimension(:), allocatable :: tem,vec2
  integer, dimension(:), allocatable :: cmgr,cmst2,cmsbgr,cmst3
  allocate(vec2(d),tem(d),cmgr(mgrp),cmst2(mgrp))
  allocate(cmsbgr(msbgrp),cmst3(msbgrp))
  if(kk==1) then
    vec2=amat(:,1); vec2(ii)=0.d0
    tem=matmul(rinv,vec2)
    do j1=1,d
      do j2=1,d
        deriv(j1,j2)= -rinv(j1,ii)*tem(j2)-rinv(j2,ii)*tem(j1) 
      end do
    end do
  endif
  if(kk==2) then
    temi=0
    do jg=1,mgrp
      cmst2(jg)=temi+1
      temi=temi+grsize(jg)
      cmgr(jg)=temi
    end do
    ! find group for ii
    do jg=1,mgrp 
      if(cmst2(jg)<=ii .and. ii<=cmgr(jg)) exit
    end do
    vec2=0.d0
    vec2(cmst2(jg):cmgr(jg))=amat(cmst2(jg):cmgr(jg),2); vec2(ii)=0.d0
    tem=matmul(rinv,vec2)
    do j1=1,d
      do j2=1,d
        deriv(j1,j2)= -rinv(j1,ii)*tem(j2)-rinv(j2,ii)*tem(j1) 
      end do
    end do
  endif
  if(kk==3) then
    temi=0
    do kg=1,msbgrp
      cmst3(kg)=temi+1
      temi=temi+sbgrsize(kg)
      cmsbgr(kg)=temi
    end do
    ! find group for ii
    do kg=1,msbgrp 
      if(cmst3(kg)<=ii .and. ii<=cmsbgr(kg)) exit
    end do
    vec2=0.d0
    vec2(cmst3(kg):cmsbgr(kg))=amat(cmst3(kg):cmsbgr(kg),3); vec2(ii)=0.d0
    tem=matmul(rinv,vec2)
    do j1=1,d
      do j2=1,d
        deriv(j1,j2)= -rinv(j1,ii)*tem(j2)-rinv(j2,ii)*tem(j1) 
      end do
    end do
  endif
  deallocate(vec2,tem,cmgr,cmst2,cmsbgr,cmst3)
  return
  end

! nllk and grad of tri-factor Gaussian with correlation matrix Robs
! inputs  
!   d = #variables
!   mgrp = #groups
!   msbgrp = #subgroups
!   rhmat = d x 3 matrix of correlations and partial correlations
!   grsize = vector of group sizes, length is mgrp
!   sbgrsize = vector of subgroup sizes, length is msbgrp
!   Robs = empirical correlation matrix
!   n = sample size
! outputs:
!   nllk = negative (Gaussian) log-likelihood 
!   grad = gradient of nllk
!subroutine trifctnllk(d,mgrp,rhmat,grsize,Robs,n,nllk,grad)
subroutine trifactnllk(d,mgrp,msbgrp,rhmat,grsize,sbgrsize,Robs,n,nllk,grad)
  implicit none
  integer d,mgrp,msbgrp,n,i,j,k,ell,j1,j2
  double precision rhmat(d,3), Robs(d,d), nllk, grad(d,3)
  integer grsize(mgrp),sbgrsize(msbgrp)
  double precision fctdet,fctldet,tr,traceder,lgdtderrho
  double precision, dimension(:), allocatable :: rhvec
  double precision, dimension(:,:), allocatable :: amat,fctmat,fctinv
  double precision, dimension(:,:), allocatable :: lgdtderload,rinvder
  double precision, dimension(:,:,:), allocatable :: jac,rinvderarr

  allocate (rhvec(d*3),amat(d,3),fctmat(d,d),fctinv(d,d))
  allocate (lgdtderload(d,3), rinvder(d,d), jac(3,3,d), rinvderarr(d,d,3))
  call pcor2load(d,3,rhmat,amat)
  rhvec=reshape(rhmat,(/d*3/))
  call trifct(d,mgrp,msbgrp,grsize,sbgrsize,rhvec(1:d),rhvec((d+1):(2*d)), &
      rhvec((2*d+1):(3*d)), fctmat,fctinv,fctdet,fctldet)
  tr=0.d0
  do i=1,d
    do j=1,d
      tr=tr+fctinv(i,j)*Robs(j,i)
    end do
  end do
  nllk=.5d0*n*(fctldet+tr+d*log(3.14159265358979323846d0*2.d0))
  call derivlogdetdloadtrif(d,mgrp,msbgrp,amat,fctinv,grsize,sbgrsize,lgdtderload)
  call jacobload2pcor(d,3,amat,rhmat,jac)
  ! loop through factors and variables
  ! deriv of nllk with pcor
  do i=1,d
    do k=1,3
      call derivRinvdloadtrif(d,mgrp,msbgrp,amat,fctinv,i,k,grsize,sbgrsize,rinvderarr(:,:,k))
    end do
    do k=1,3
      ! deriv of rinv wrt pcor[i,k]
      rinvder=0.d0
      lgdtderrho=0.d0
      do ell=k,3
        lgdtderrho=lgdtderrho+lgdtderload(i,ell)*jac(ell,k,i)
        rinvder=rinvder+rinvderarr(:,:,ell)*jac(ell,k,i)
      end do
      traceder=0.d0
      do j1=1,d
        do j2=1,d
          traceder=traceder+rinvder(j1,j2)*Robs(j2,j1)
        end do
      end do
      grad(i,k)=.5*n*(lgdtderrho+traceder)
    end do
  end do
  deallocate(amat,fctinv,fctmat,rhvec)
  deallocate(lgdtderload,jac,rinvder,rinvderarr)
  return
  end

! mvt with tri-factor correlation structure
! nllk and grad of tri-factor Gaussian with correlation matrix Robs
! inputs  
!   d = #variables
!   mgrp = #groups
!   msbgrp = #subgroups
!   grsize = vector of group sizes, length is mgrp
!   sbgrsize = vector of subgroup sizes, length is msbgrp
!   n = sample size
!   nu = degree of freedom parameter >0
!   rhmat = d x 3 matrix of correlations and partial correlations
!   tdata = nxd data set of t-scores
! outputs:
!   nllk = negative log-likelihood for copula
!   grad = gradient of nllk
subroutine ttrifactnllk(d,mgrp,msbgrp,grsize,sbgrsize,n,nu,rhmat,tdata,nllk,grad)
  implicit none
  integer d,mgrp,msbgrp,n,i,k,ell,j1,j2,ii
  double precision nu, rhmat(d,3), tdata(n,d), nllk, grad(d,3)
  integer grsize(mgrp),sbgrsize(msbgrp)
  double precision fctdet,fctldet,qf,qfder,gr,lgdtderrho
  double precision const,ulogpdf,nu1,nud
  double precision, dimension(:), allocatable :: rhvec,qfvec
  double precision, dimension(:,:), allocatable :: amat,fctmat,fctinv
  double precision, dimension(:,:), allocatable :: lgdtderload,rinvder
  double precision, dimension(:,:,:), allocatable :: jac,rinvderarr
  double precision lgamma   ! function in C code lgamma.c, in gnu fortran?

  allocate (rhvec(d*3),amat(d,3),fctmat(d,d),fctinv(d,d))
  allocate (lgdtderload(d,3), rinvder(d,d), jac(3,3,d), rinvderarr(d,d,3))
  allocate (qfvec(n))
  nu1=(nu+1.d0)/2.d0
  nud=(nu+d)/2.d0
  const=lgamma(nud)+(d-1.d0)*lgamma(0.5d0*nu)-d*lgamma(nu1);
  call pcor2load(d,3,rhmat,amat)
  rhvec=reshape(rhmat,(/d*3/))
  call trifct(d,mgrp,msbgrp,grsize,sbgrsize,rhvec(1:d),rhvec((d+1):(2*d)), &
      rhvec((2*d+1):(3*d)), fctmat,fctinv,fctdet,fctldet)
  nllk=0.d0
  do ii=1,n ! loop of data vectors
    qf=0.d0 ; ulogpdf=0.d0
    do j1=1,d
      do j2=1,d
        qf=qf+fctinv(j1,j2)*tdata(ii,j1)*tdata(ii,j2)
      end do
      ulogpdf=ulogpdf-nu1*log(1.d0+(tdata(ii,j1)**2)/nu)
    end do
    qfvec(ii)=nu+qf ! to use for gradient
    nllk=nllk+nud*log(1.d0+qf/nu)+ulogpdf
  end do
  nllk=nllk-n*const+n*0.5d0*fctldet
  call derivlogdetdloadtrif(d,mgrp,msbgrp,amat,fctinv,grsize,sbgrsize,lgdtderload)
  call jacobload2pcor(d,3,amat,rhmat,jac)
  ! loop through factors and variables
  ! deriv of nllk with pcor
  do i=1,d
    do k=1,3
      call derivRinvdloadtrif(d,mgrp,msbgrp,amat,fctinv,i,k,grsize,sbgrsize,rinvderarr(:,:,k))
    end do
    do k=1,3
      ! deriv of rinv wrt pcor[i,k]
      rinvder=0.d0
      ! deriv of log det wrt pcor[i,k] is lgdtderrho
      lgdtderrho=0.d0
      do ell=k,3
        lgdtderrho=lgdtderrho+lgdtderload(i,ell)*jac(ell,k,i)
        rinvder=rinvder+rinvderarr(:,:,ell)*jac(ell,k,i)
      end do
      gr=0.d0
      do ii=1,n ! loop of data vectors
        qfder=0.d0 ! for derivative of quadratic from wrt pcor[i,k]
        do j1=1,d
          do j2=1,d
            qfder=qfder+rinvder(j1,j2)*tdata(ii,j1)*tdata(ii,j2)
          end do
        end do
        gr=gr+qfder*nud/qfvec(ii)
      end do
      grad(i,k)=.5*n*lgdtderrho+gr
    end do
  end do
  deallocate(amat,fctinv,fctmat,rhvec)
  deallocate(lgdtderload,jac,rinvder,rinvderarr)
  deallocate(qfvec)
  return
  end

