
! multivariate Gaussian/t model with p-factor and bi-factor structures
! function for the gradient of the negative log-likelihood
! multivariate t(nu) nllk has U(0,1) margin;
! multivariate Gaussian nllk has N(0,1) margin 
!          (for R functions factanal.com factanal.co)
!   and converted to nllk with U(0,1) margin 
!          (in R functions mvtpfactnllk,mvtbifactnllk)

! Some details on how gradients are computed are in:
! Krupskii (2014). Structured Factor Copulas and Tail Inference.
!  PhD thesis, University of British Columbia.

! The earlier version of this code was written in R by Pavel Krupskii.
! This fortran version is more readable because the R version needed
! vectorization of 3-dimensional arrays for speed.

!program mainex
!  implicit none
!  integer d,p,i,k,j,n
!  double precision, dimension(:), allocatable :: rh
!  double precision, dimension(:,:), allocatable :: fctmat,fctinv,grad,Robs
!  double precision, dimension(:,:), allocatable :: rhmat,amat,lgdtder,rinvder
!  double precision, dimension(:,:,:), allocatable :: jarr
!  double precision fctdet,fctldet,nllk
!
!  n=50
!  read *,d,p
!  do while(p>=1) 
!    allocate (rh(d*p),fctmat(d,d),fctinv(d,d),rhmat(d,p),amat(d,p))
!    allocate (lgdtder(d,p), jarr(p,p,d), rinvder(d,d), grad(d,p), Robs(d,d))
!    do i=1,d*p
!      rh(i)=i/(d*p+1.d0)
!    end do
!    do i=1,d
!      do j=1,d
!        Robs(i,j)=0.6d0**(abs(i-j))
!      end do
!    end do
!    !print '(5f12.8)', Robs(1,:)
!    call pfct(d,p,rh,fctmat,fctinv,fctdet,fctldet)
!    do i=1,d
!      print '(8f10.6)',fctmat(i,:)
!    end do
!    print *, fctdet, fctldet
!    do i=1,d
!      print '(8f10.6)',fctinv(i,:)
!    end do
!    do k=1,p
!      rhmat(:,k)=rh(((k-1)*d+1):(k*d))
!    end do
!    call pcor2load(d,p,rhmat,amat)
!    !do i=1,d
!    !  print '(5f12.8)',amat(i,:)
!    !end do
!    call derivlogdetdload(d,p,amat,fctinv,lgdtder)
!    do i=1,d
!      print '(5f12.8)',lgdtder(i,:)
!    end do
!    call jacobload2pcor(d,p,amat,rhmat,jarr)
!    !do i=1,d
!    !  print '(9f11.7)',jarr(:,:,i)
!    !end do
!    call derivRinvdload(d,p,amat,fctinv,4,1,rinvder)
!    print *, ""
!    do i=1,d
!      print '(5f13.8)',rinvder(i,:)
!    end do
!    !call derivRinvdload(d,p,amat,fctinv,3,2,rinvder)
!    !print *, ""
!    !do i=1,d
!    !  print '(5f13.8)',rinvder(i,:)
!    !end do
!    !call derivRinvdload(d,p,amat,fctinv,2,3,rinvder)
!    !print *, ""
!    !do i=1,d
!    !  print '(5f13.8)',rinvder(i,:)
!    !end do
!    print *, "deriv of nllk"
!    call pfactnllk(d,p,rhmat,Robs,n,nllk,grad)
!    print *,nllk
!    do i=1,d
!      print '(5f13.8)',grad(i,:)
!    end do
!    deallocate (rh,fctmat,fctinv, rhmat,amat, lgdtder, jarr, rinvder,grad,Robs)
!    read *,d,p
!  end do
!  stop
!  end

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


! computes matrix, its determinant and inverse in a p-factor normal model
! inputs  
!   d = #variables
!   p = #factors 
!   rh = p*d vector of *conditional* correlations (1st factor, 2nd etc);
! outputs 
!   fctmat = covariance matrix in a p-factor model with parameters rh
!   fctdet = det(fctmat)   
!   fctinv = inverse of fctmat 
subroutine pfct(d,p,rh,fctmat,fctinv,fctdet,fctldet)
  implicit none
  integer p,d,i,j,k1,k2,k
  double precision rh(p*d),fctmat(d,d),fctinv(d,d)
  double precision fctdet,fctldet
  double precision, dimension(:), allocatable :: rhe,irhe,rh1,rh2,rh3,a2,a3,psiv
  double precision tem,tol,bdet
  double precision, dimension(:,:), allocatable :: rhmat,tmat,invtmat,tmat2,tmatd 
  double precision, dimension(:,:), allocatable :: amat,bmat, soln
  double precision det2,det3
  allocate ( rhe(d),irhe(d),rh1(d),rh2(d),rh3(d),a2(d),a3(d) ) 
  allocate ( rhmat(d,p), tmat(p,p), invtmat(p,p), tmat2(p,d),tmatd(d,d) )
  allocate ( amat(d,p), bmat(p,p), soln(p,d), psiv(d) )
  if(p==1) then
    rhe = 1.d0-rh*rh;
    irhe= 1.d0/rhe;
    call outer(d,d,rh,rh,fctmat);
    tem=sum(rh*rh/rhe) +1.d0
    fctdet = product(rhe)*tem; 
    fctldet=log(fctdet) 
    fctinv=0.d0
    do i=1,d
      fctinv(i,i)=irhe(i)
      do j=1,d
        fctinv(i,j)=fctinv(i,j)-irhe(i)*fctmat(i,j)*irhe(j)/tem
      end do
    end do
    do i=1,d
      fctmat(i,i)=1.d0
    end do
  end if
  if(p==2) then
    rh1=rh(1:d); rh2=rh((d+1):(2*d));
    a2=rh2*sqrt(1.d0-rh1**2);
    rhe = 1.d0-rh1**2-a2**2;
    irhe= 1.d0/rhe;
    rhmat(:,1)=rh1; rhmat(:,2)=a2;
    fctmat=matmul(rhmat,transpose(rhmat));
    tmat=0.d0; 
    do k1=1,p
      do k2=1,p
        tem=0.d0;
        do i=1,d
          tem=tem+rhmat(i,k1)*rhmat(i,k2)*irhe(i)
        end do
        tmat(k1,k2)=tem
      end do
      tmat(k1,k1)=tmat(k1,k1)+1
    end do
    fctdet = product(rhe)*det2(tmat); 
    fctldet=log(fctdet) 
    call inv2(tmat,invtmat)
    tmat2=matmul(invtmat,transpose(rhmat))
    tmatd=matmul(rhmat,tmat2)
    fctinv=0.d0
    do i=1,d
      fctinv(i,i)=irhe(i)
      fctmat(i,i)=1.d0
      do j=1,d
        fctinv(i,j)=fctinv(i,j)-irhe(i)*tmatd(i,j)*irhe(j)
      end do
    end do
  end if
  if(p==3) then
    rh1 = rh(1:d); rh2= rh((d+1):(2*d)); rh3= rh((2*d+1):(3*d));
    a3=rh3*sqrt(1.d0-rh1**2)*sqrt(1-rh2**2);
    a2=rh2*sqrt(1.d0-rh1**2);
    rhe = 1.d0-rh1**2-a2**2-a3**2;
    irhe= 1.d0/rhe;
    rhmat(:,1)=rh1; rhmat(:,2)=a2; rhmat(:,3)=a3
    fctmat=matmul(rhmat,transpose(rhmat));
    tmat=0.d0; 
    do k1=1,p
      do k2=1,p
        tem=0.d0;
        do i=1,d
          tem=tem+rhmat(i,k1)*rhmat(i,k2)*irhe(i)
        end do
        tmat(k1,k2)=tem
      end do
      tmat(k1,k1)=tmat(k1,k1)+1
    end do
    fctdet = product(rhe)*det3(tmat); 
    fctldet=log(fctdet) 
    call inv3(tmat,invtmat)
    tmat2=matmul(invtmat,transpose(rhmat))
    tmatd=matmul(rhmat,tmat2)
    fctinv=0.d0
    do i=1,d
      fctinv(i,i)=irhe(i)
      fctmat(i,i)=1.d0
      do j=1,d
        fctinv(i,j)=fctinv(i,j)-irhe(i)*tmatd(i,j)*irhe(j)
      end do
    end do
  endif
  ! p>=4 use linear system solver lsolve=gepp
  if(p>=4) then
    tol=1.d-8  ! tolerance for singular matrix
    rhmat=reshape(rh,(/d,p/))
    call pcor2load(d,p,rhmat,amat)
    fctmat=matmul(amat,transpose(amat))
    do i=1,d
      tem=0.d0
      do k=1,p
        tem=tem+amat(i,k)**2
      end do
      psiv(i)=1.d0-tem
      fctmat(i,i)=1.d0
    end do
    do k=1,p  ! tmat2=t(amat)%*%diag(1/psiv) is pxd
      do i=1,d
        tmat2(k,i)=amat(i,k)/psiv(i)
      end do
    end do
    bmat=matmul(tmat2,amat) ! bmat is pxp
    do k=1,p
      bmat(k,k)=bmat(k,k)+1.d0
    end do
    !do k=1,p
    !  print *,tmat2(k,:)
    !end do
    call lsolve(p,d,bmat,tmat2,tol,soln,bdet)
    !do k=1,p
    !  print *,tmat2(k,:) 
    !end do
    fctinv=matmul(transpose(tmat2),soln)
    do i=1,d
      fctinv(i,i)=fctinv(i,i)-1.d0/psiv(i)
    end do
    fctinv=-fctinv
    fctldet=log(bdet)+sum(log(psiv))
    fctdet=exp(fctldet)
  endif
  deallocate ( rhe,irhe, rh1,rh2,rh3, a2,a3, rhmat,tmat,invtmat,tmat2,tmatd )
  deallocate ( amat,bmat,psiv,soln )
  return
  end

! input: M = 2x2 matrix
! output: determinant of M
double precision function det2(M) 
  implicit none
  double precision M(2,2)
  det2=M(1,1)*M(2,2)-M(1,2)**2 
  return
  end

! input: M = 3x3 matrix
! output: determinant of M
double precision function det3(M) 
  implicit none
  double precision M(3,3)
  det3 = M(1,1)*M(2,2)*M(3,3) + 2*M(1,2)*M(1,3)*M(2,3) - M(1,3)**2 *M(2,2)  &
   - M(2,3)**2 *M(1,1) - M(1,2)**2 *M(3,3)
  return
  end

! input:  M = 2x2 matrix
! output: invM = inverse of M
subroutine inv2(M,invM)
  implicit none
  double precision M(2,2),invM(2,2),det2,dm
  dm = det2(M)
  if(abs(dm)<1.d-12)  then
    !print *,"matrix is computationally singular"; 
    call intpr("matrix is computationally singular",34,0,0)
    invM=1.d10
    return
  end if
  invM(1,2) = -M(1,2);
  invM(2,1) = invM(1,2);
  invM(1,1) = M(2,2);
  invM(2,2) = M(1,1);
  invM = invM/dm;
  return
  end

! input:  M = 3x3 matrix
! output: invM = inverse of M
subroutine inv3(M,invM)
  implicit none
  double precision M(3,3),invM(3,3),det3,dm
  dm = det3(M);
  if(abs(dm)<1.d-12)  then
    !print *,"matrix is computationally singular"; 
    call intpr("matrix is computationally singular",34,0,0)
    invM=1.d10
    return
  end if
  invM(1,1) = M(2,2)*M(3,3)-M(2,3)**2;
  invM(2,2) = M(1,1)*M(3,3)-M(1,3)**2;
  invM(3,3) = M(1,1)*M(2,2)-M(1,2)**2;
  invM(1,2) = M(1,3)*M(2,3)-M(1,2)*M(3,3);
  invM(1,3) = M(1,2)*M(2,3)-M(1,3)*M(2,2);
  invM(2,3) = M(1,2)*M(1,3)-M(2,3)*M(1,1);
  invM(2,1) = invM(1,2);
  invM(3,1) = invM(1,3);
  invM(3,2) = invM(2,3);
  invM = invM/dm;
  return
  end

! partial correlation representation to loadings for p-factor with d variables
! inputs  
!   d = #variables
!   p = #factors 
!   rhmat = d x p matrix of conditional correlations (1st factor, 2nd etc);
! output 
!   amat = d x p matrix of loadings
subroutine pcor2load(d,p,rhmat,amat)
  implicit none
  integer d,p,kf
  double precision rhmat(d,p), amat(d,p)
  double precision, dimension(:), allocatable :: tem,temc
  if(p==1) then
    amat=rhmat
    return
  end if
  ! p>=2
  allocate (tem(d),temc(d))
  amat(:,1)=rhmat(:,1)
  tem=1.d0
  do kf=2,p
    temc=sqrt(1.d0-rhmat(:,kf-1)**2)
    tem=tem*temc
    amat(:,kf)=rhmat(:,kf)*tem 
  end do
  deallocate (tem,temc)
  return
  end

! partial derivatives of log det for p-factor
! inputs  
!   d = #variables
!   p = #factors 
!   amat = d x p matrix of loadings
!   rinv = inverse of correlation matrix
! output 
!   lgdtder = d x p matrix with derivative of log det wrt each loading 
subroutine derivlogdetdload(d,p,amat,rinv,lgdtder)
  implicit none
  integer d,p,ii,kf
  double precision amat(d,p), rinv(d,d), lgdtder(d,p) 
  double precision, dimension(:), allocatable :: vec1,vec2
  allocate (vec1(d),vec2(d))
  do ii=1,d
    vec1=rinv(:,ii) 
    do kf=1,p
      vec2=amat(:,kf); vec2(ii)=0.d0
      lgdtder(ii,kf)=2.d0*sum(vec1*vec2)
    end do
  end do
  deallocate (vec1,vec2)
  return 
  end

! Jacobian of load wrt to pcor for a single row i of loading matrix
! inputs  
!   d = #variables
!   p = #factors 
!   amat = d x p matrix of loadings
!   rhmat = d x p matrix of conditional correlations (1st factor, 2nd etc);
! output 
!   jarr = 3-dimensional array, each [,,i] is lower triangular and its
!       [k,ell] position has derivative of a_{ik} wrt pcor_{i,ell}
subroutine jacobload2pcor(d,p,amat,rhmat,jarr)
  implicit none
  integer d,p,ell,k
  double precision amat(d,p), rhmat(d,p), jarr(p,p,d) 
  double precision, dimension(:), allocatable :: tem
  
  allocate (tem(d))
  tem=1.d0
  do ell=1,p
    ! problem with next line if denominator is zero
    !jarr(ell,ell,:)=amat(:,ell)/rhmat(:,ell)
    jarr(ell,ell,:)=tem
    tem=tem*sqrt(1-rhmat(:,ell)**2)
  end do
  do ell=1,(p-1)
    do k=(ell+1),p
      jarr(k,ell,:)=-rhmat(:,ell)*amat(:,k)/(1-rhmat(:,ell)**2)
    end do
  end do
  deallocate (tem)
  return 
  end

! derivative of rinv wrt a loading parameter
! inputs  
!   d = #variables
!   p = #factors 
!   amat = d x p matrix of loadings
!   rinv = inverse of correlation matrix
!   ii = index of variable. 1<=ii<=d
!   kk = index of factor. 1<=kk<=p
! output 
!   deriv = dxd matrix with derivative of components of rinv with load[ii,kk] 
subroutine derivRinvdload(d,p,amat,rinv,ii,kk,deriv)
  implicit none
  integer d,p,ii,kk,j1,j2
  double precision amat(d,p), rinv(d,d), deriv(d,d)
  double precision, dimension(:), allocatable :: tem,vec2
  allocate(vec2(d),tem(d))
  vec2=amat(:,kk); vec2(ii)=0.d0
  tem=matmul(rinv,vec2)
  do j1=1,d
    do j2=1,d
      deriv(j1,j2)= -rinv(j1,ii)*tem(j2)-rinv(j2,ii)*tem(j1) 
    end do
  end do
  deallocate(vec2,tem)
  return
  end

! nllk and grad of p-factor Gaussian with correlation matrix Robs
! inputs  
!   d = #variables
!   p = #factors 
!   rhmat = d x p matrix of partial/conditional correlations
!   Robs = empirical correlation matrix
!   n = sample size
! outputs:
!   nllk = negative (Gaussian) log-likelihood 
!   grad = gradient of nllk
!subroutine pfctnllk(d,p,rhmat,Robs,n,nllk,grad)
subroutine pfactnllk(d,p,rhmat,Robs,n,nllk,grad)
  implicit none
  integer d,p,n,i,j,k,ell,j1,j2
  double precision rhmat(d,p), Robs(d,d), nllk, grad(d,p)
  double precision fctdet,fctldet,tr,traceder,lgdtderrho
  double precision, dimension(:), allocatable :: rhvec
  double precision, dimension(:,:), allocatable :: amat,fctmat,fctinv
  double precision, dimension(:,:), allocatable :: lgdtderload,rinvder
  double precision, dimension(:,:,:), allocatable :: jac,rinvderarr

  allocate (rhvec(d*p),amat(d,p),fctmat(d,d),fctinv(d,d))
  allocate (lgdtderload(d,p), rinvder(d,d), jac(p,p,d), rinvderarr(d,d,p))
  call pcor2load(d,p,rhmat,amat)
  rhvec=reshape(rhmat,(/d*p/))
  call pfct(d,p,rhvec,fctmat,fctinv,fctdet,fctldet)
  tr=0.d0
  do i=1,d
    do j=1,d
      tr=tr+fctinv(i,j)*Robs(j,i)
    end do
  end do
  nllk=.5d0*n*(fctldet+tr+d*log(3.14159265358979323846d0*2.d0))
  call derivlogdetdload(d,p,amat,fctinv,lgdtderload)
  call jacobload2pcor(d,p,amat,rhmat,jac)
  ! loop through factors and variables
  ! deriv of nllk with pcor
  do i=1,d
    do k=1,p
      call derivRinvdload(d,p,amat,fctinv,i,k,rinvderarr(:,:,k))
    end do
    do k=1,p
      ! deriv of rinv wrt pcor[i,k]
      rinvder=0.d0
      lgdtderrho=0.d0
      do ell=k,p
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

!------------------------------------------------------------

! nllk and grad of mvt with p-factor correlation matrix 
! inputs  
!   d = #variables
!   p = #factors 
!   n = sample size
!   nu = degree of freedom parameter >0
!   rhmat = d x p matrix of partial/conditional correlations
!   tdata = nxd data set of t-scores
! outputs:
!   nllk = negative log-likelihood for copula 
!   grad = gradient of nllk
subroutine tpfactnllk(d,p,n,nu,rhmat,tdata,nllk,grad)
  implicit none
  integer d,p,n,i,k,ell,j1,j2,ii
  double precision nu, rhmat(d,p), tdata(n,d), nllk, grad(d,p)
  double precision fctdet,fctldet,qf,qfder,gr,lgdtderrho
  double precision const,ulogpdf,nu1,nud
  double precision, dimension(:), allocatable :: rhvec,qfvec
  double precision, dimension(:,:), allocatable :: amat,fctmat,fctinv
  double precision, dimension(:,:), allocatable :: lgdtderload,rinvder
  double precision, dimension(:,:,:), allocatable :: jac,rinvderarr
  double precision lgamma   ! function in C code lgamma.c, in gnu fortran?

  allocate (rhvec(d*p),amat(d,p),fctmat(d,d),fctinv(d,d))
  allocate (lgdtderload(d,p), rinvder(d,d), jac(p,p,d), rinvderarr(d,d,p))
  allocate (qfvec(n))
  nu1=(nu+1.d0)/2.d0
  nud=(nu+d)/2.d0
  const=lgamma(nud)+(d-1.d0)*lgamma(0.5d0*nu)-d*lgamma(nu1);
  call pcor2load(d,p,rhmat,amat)
  rhvec=reshape(rhmat,(/d*p/))
  call pfct(d,p,rhvec,fctmat,fctinv,fctdet,fctldet)
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
  call derivlogdetdload(d,p,amat,fctinv,lgdtderload) ! deriv of logdet
  call jacobload2pcor(d,p,amat,rhmat,jac)
  ! loop through factors and variables
  ! deriv of nllk with pcor
  do i=1,d
    do k=1,p
      call derivRinvdload(d,p,amat,fctinv,i,k,rinvderarr(:,:,k))
    end do
    do k=1,p
      ! deriv of rinv wrt pcor[i,k]
      rinvder=0.d0
      ! deriv of log det wrt pcor[i,k] is lgdtderrho
      lgdtderrho=0.d0
      do ell=k,p
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

!-------------------- extra for bifactor below-----------------------------

!computes matrix, its determinant and inverse in a bi-factor normal model
!         rh1,rh2 - *conditional* correlations (1st factor, 2nd );
! inputs  
!   d = #variables
!   mgrp = #groups
!   grsize = vector of group sizes, length is mgrp
!   rh1 = d-vector of correlations with common latent
!   rh2 = d-vector of partial correlations with group latent given common latent
! outputs 
!   fctmat = covariance matrix in bi-factor model with parameters rh1,rh2
!   fctinv = inverse of fctmat 
!   fctdet = det(fctmat)   
!   fctldet = log ( det(fctmat) )  
subroutine bifct(d,mgrp,grsize,rh1,rh2,fctmat,fctinv,fctdet,fctldet)
  implicit none
  integer d,mgrp,jg,i,temi
  integer grsize(mgrp)
  double precision rh1(d),rh2(d),fctmat(d,d),fctinv(d,d)
  double precision fctdet,fctldet,A0
  double precision, dimension(:), allocatable :: a2,A,Aj,eta,etaa,dzeta1,dzeta2,xi
  double precision, dimension(:,:), allocatable :: fctdiag,tmat 
  integer, dimension(:), allocatable :: cmgr,cmst
  allocate ( a2(d),dzeta1(d),dzeta2(d),xi(d)) 
  allocate (cmgr(mgrp),cmst(mgrp),A(mgrp),Aj(d),eta(mgrp),etaa(d))
  allocate ( fctdiag(d,d), tmat(d,d) )
  temi=0
  do jg=1,mgrp
    cmst(jg)=temi+1
    temi=temi+grsize(jg)
    cmgr(jg)=temi
  end do
  a2=rh2*sqrt(1.d0-rh1**2);
  A=0.d0; Aj=0.d0
  eta=0.d0; etaa=0.d0
  dzeta1=rh1**2/(1-rh1**2)/(1-rh2**2)
  dzeta2=a2/(1-rh1**2)/(1-rh2**2)
  fctdiag=0.d0
  do jg=1,mgrp
    fctdiag(cmst(jg):cmgr(jg),cmst(jg):cmgr(jg))=1.d0
    A(jg)=1-grsize(jg)+sum(1.d0/(1-rh2(cmst(jg):cmgr(jg))**2))
    Aj(cmst(jg):cmgr(jg))=A(jg)
    eta(jg)=sum(dzeta2(cmst(jg):cmgr(jg))*rh1(cmst(jg):cmgr(jg)))
    etaa(cmst(jg):cmgr(jg))=eta(jg)/A(jg)
  end do
  call outer(d,d,rh1,rh1,fctmat)
  call outer(d,d,a2,a2,tmat)
  fctmat=fctmat+tmat*fctdiag
  do i=1,d
    fctmat(i,i)=1.d0
  end do
  xi = dzeta1/rh1-etaa*dzeta2
  A0 = 1.d0-sum(eta**2/A)+sum(dzeta1)
  fctdet = A0*product(A)*product((1-rh1**2)*(1-rh2**2))
  call outer(d,d,xi,xi,fctinv)
  call outer(d,d,dzeta2,dzeta2,tmat)
  do i=1,d
    tmat(:,i)=tmat(:,i)/Aj
    tmat(:,i)=tmat(:,i)*fctdiag(:,i)
  end do
  fctinv= -fctinv/A0 - tmat
  do i=1,d
    fctinv(i,i)=1/(1-rh1(i)**2)/(1-rh2(i)**2) -dzeta2(i)**2/Aj(i)-xi(i)**2/A0
  end do
  deallocate (a2,dzeta1,dzeta2,xi,fctdiag,tmat )
  deallocate (cmgr,cmst,A,Aj,eta,etaa)
  fctldet=log(fctdet) 
  return
  end

! partial derivatives of log det for bi-factor
! inputs  
!   d = #variables
!   mgrp = #groups
!   amat = d x 2 matrix of loadings
!   rinv = inverse of correlation matrix
!   grsize = vector of group sizes, length is mgrp
! output 
!   lgdtder = d x 2 matrix with derivative of log det wrt each loading 
subroutine derivlogdetdloadbif(d,mgrp,amat,rinv,grsize,lgdtder)
  implicit none
  integer d,mgrp,ii,jg,temi
  double precision amat(d,2), rinv(d,d), lgdtder(d,2) 
  integer grsize(mgrp)
  double precision, dimension(:), allocatable :: vec1,vec2
  integer, dimension(:), allocatable :: cmgr,cmst
  allocate (vec1(d),vec2(d),cmgr(mgrp),cmst(mgrp))
  temi=0
  do jg=1,mgrp
    cmst(jg)=temi+1
    temi=temi+grsize(jg)
    cmgr(jg)=temi
  end do
  ! common factor in col 1
  do ii=1,d
    vec1=rinv(:,ii) 
    vec2=amat(:,1); vec2(ii)=0.d0
    lgdtder(ii,1)=2.d0*sum(vec1*vec2)
  end do
  ! group factor in col 2
  do jg=1,mgrp
      vec1=0.d0
    do ii=cmst(jg),cmgr(jg)
      vec1(cmst(jg):cmgr(jg))=rinv(cmst(jg):cmgr(jg),ii)
      vec2(cmst(jg):cmgr(jg))=amat(cmst(jg):cmgr(jg),2)
      vec2(ii)=0.d0
      lgdtder(ii,2)=2.d0*sum(vec1(cmst(jg):cmgr(jg))*vec2(cmst(jg):cmgr(jg)))
    end do
  end do
  deallocate (vec1,vec2,cmgr,cmst)
  return 
  end

! derivative of rinv wrt a loading parameter
! inputs  
!   d = #variables
!   mgrp = #groups
!   amat = d x 2 matrix of loadings
!   rinv = inverse of correlation matrix
!   rinv = inverse of correlation matrix
!   ii = index of variable. 1<=ii<=d
!   kk = index of factor. 1<=kk<=2
!   grsize = vector of group sizes, length is mgrp
! output 
!   deriv = dxd matrix with derivative of components of rinv with load[ii,kk] 
subroutine derivRinvdloadbif(d,mgrp,amat,rinv,ii,kk,grsize,deriv)
  implicit none
  integer d,mgrp,ii,kk,j1,j2,jg,temi
  double precision amat(d,2), rinv(d,d), deriv(d,d)
  integer grsize(mgrp)
  double precision, dimension(:), allocatable :: tem,vec2
  integer, dimension(:), allocatable :: cmgr,cmst
  allocate(vec2(d),tem(d),cmgr(mgrp),cmst(mgrp))
  if(kk==1) then
    vec2=amat(:,1); vec2(ii)=0.d0
    tem=matmul(rinv,vec2)
    do j1=1,d
      do j2=1,d
        deriv(j1,j2)= -rinv(j1,ii)*tem(j2)-rinv(j2,ii)*tem(j1) 
      end do
    end do
  else ! kk=2
    temi=0
    do jg=1,mgrp
      cmst(jg)=temi+1
      temi=temi+grsize(jg)
      cmgr(jg)=temi
    end do
    ! find group for ii
    do jg=1,mgrp 
      if(cmst(jg)<=ii .and. ii<=cmgr(jg)) exit
    end do
    vec2=0.d0
    vec2(cmst(jg):cmgr(jg))=amat(cmst(jg):cmgr(jg),2); vec2(ii)=0.d0
    tem=matmul(rinv,vec2)
    do j1=1,d
      do j2=1,d
        deriv(j1,j2)= -rinv(j1,ii)*tem(j2)-rinv(j2,ii)*tem(j1) 
      end do
    end do
  endif
  deallocate(vec2,tem,cmgr,cmst)
  return
  end

! nllk and grad of bi-factor Gaussian with correlation matrix Robs
! inputs  
!   d = #variables
!   mgrp = #groups
!   rhmat = d x 2 matrix of correlations and partial correlations
!   grsize = vector of group sizes, length is mgrp
!   Robs = empirical correlation matrix
!   n = sample size
! outputs:
!   nllk = negative (Gaussian) log-likelihood 
!   grad = gradient of nllk
!subroutine bifctnllk(d,mgrp,rhmat,grsize,Robs,n,nllk,grad)
subroutine bifactnllk(d,mgrp,rhmat,grsize,Robs,n,nllk,grad)
  implicit none
  integer d,mgrp,n,i,j,k,ell,j1,j2
  double precision rhmat(d,2), Robs(d,d), nllk, grad(d,2)
  integer grsize(mgrp)
  double precision fctdet,fctldet,tr,traceder,lgdtderrho
  double precision, dimension(:), allocatable :: rhvec
  double precision, dimension(:,:), allocatable :: amat,fctmat,fctinv
  double precision, dimension(:,:), allocatable :: lgdtderload,rinvder
  double precision, dimension(:,:,:), allocatable :: jac,rinvderarr

  allocate (rhvec(d*2),amat(d,2),fctmat(d,d),fctinv(d,d))
  allocate (lgdtderload(d,2), rinvder(d,d), jac(2,2,d), rinvderarr(d,d,2))
  call pcor2load(d,2,rhmat,amat)
  rhvec=reshape(rhmat,(/d*2/))
  call bifct(d,mgrp,grsize,rhvec(1:d),rhvec((d+1):(2*d)), &
      fctmat,fctinv,fctdet,fctldet)
  tr=0.d0
  do i=1,d
    do j=1,d
      tr=tr+fctinv(i,j)*Robs(j,i)
    end do
  end do
  nllk=.5d0*n*(fctldet+tr+d*log(3.14159265358979323846d0*2.d0))
  call derivlogdetdloadbif(d,mgrp,amat,fctinv,grsize,lgdtderload)
  call jacobload2pcor(d,2,amat,rhmat,jac)
  ! loop through factors and variables
  ! deriv of nllk with pcor
  do i=1,d
    do k=1,2
      call derivRinvdloadbif(d,mgrp,amat,fctinv,i,k,grsize,rinvderarr(:,:,k))
    end do
    do k=1,2
      ! deriv of rinv wrt pcor[i,k]
      rinvder=0.d0
      lgdtderrho=0.d0
      do ell=k,2
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

! mvt with bi-factor correlation structure
! nllk and grad of bi-factor Gaussian with correlation matrix Robs
! inputs  
!   d = #variables
!   mgrp = #groups
!   grsize = vector of group sizes, length is mgrp
!   n = sample size
!   nu = degree of freedom parameter >0
!   rhmat = d x 2 matrix of correlations and partial correlations
!   tdata = nxd data set of t-scores
! outputs:
!   nllk = negative log-likelihood for copula
!   grad = gradient of nllk
subroutine tbifactnllk(d,mgrp,grsize,n,nu,rhmat,tdata,nllk,grad)
  implicit none
  integer d,mgrp,n,i,k,ell,j1,j2,ii
  double precision nu, rhmat(d,2), tdata(n,d), nllk, grad(d,2)
  integer grsize(mgrp)
  double precision fctdet,fctldet,qf,qfder,gr,lgdtderrho
  double precision const,ulogpdf,nu1,nud
  double precision, dimension(:), allocatable :: rhvec,qfvec
  double precision, dimension(:,:), allocatable :: amat,fctmat,fctinv
  double precision, dimension(:,:), allocatable :: lgdtderload,rinvder
  double precision, dimension(:,:,:), allocatable :: jac,rinvderarr
  double precision lgamma   ! function in C code lgamma.c, in gnu fortran?

  allocate (rhvec(d*2),amat(d,2),fctmat(d,d),fctinv(d,d))
  allocate (lgdtderload(d,2), rinvder(d,d), jac(2,2,d), rinvderarr(d,d,2))
  allocate (qfvec(n))
  nu1=(nu+1.d0)/2.d0
  nud=(nu+d)/2.d0
  const=lgamma(nud)+(d-1.d0)*lgamma(0.5d0*nu)-d*lgamma(nu1);
  call pcor2load(d,2,rhmat,amat)
  rhvec=reshape(rhmat,(/d*2/))
  call bifct(d,mgrp,grsize,rhvec(1:d),rhvec((d+1):(2*d)), &
      fctmat,fctinv,fctdet,fctldet)
  !print *, fctldet,const
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
  call derivlogdetdloadbif(d,mgrp,amat,fctinv,grsize,lgdtderload)
  call jacobload2pcor(d,2,amat,rhmat,jac)
  ! loop through factors and variables
  ! deriv of nllk with pcor
  do i=1,d
    do k=1,2
      call derivRinvdloadbif(d,mgrp,amat,fctinv,i,k,grsize,rinvderarr(:,:,k))
    end do
    do k=1,2
      ! deriv of rinv wrt pcor[i,k]
      rinvder=0.d0
      ! deriv of log det wrt pcor[i,k] is lgdtderrho
      lgdtderrho=0.d0
      do ell=k,2
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

