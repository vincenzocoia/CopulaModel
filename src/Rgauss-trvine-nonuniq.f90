! best Gaussian truncated vines with 1,2,3,...d-2 trees
! further checks on non-uniqueness
! do the following if main program is not commented
!   and change intpr and dblepr to print statements
! gfortran -o gauss-trvine-nonuniq Rgauss-trvine-nonuniq.f90 varray.f90 pcor.f90

! gfortran -fpic -c Rgauss-trvine-nonuniq.f90 ../src/varray.f90 ../src/pcor.f90
! gfortran -shared -o trvineuniq.so Rgauss-trvine-nonuniq.o varray.o pcor.o

!program mainvarr
!  ! read in correlation matrix and pass to subroutine
!  implicit none
!  integer d,dd,i,ell,iprint,jtrunc
!  integer, dimension(:), allocatable :: bnumv
!  integer, dimension(:,:), allocatable :: permmat
!  double precision, dimension(:,:), allocatable :: r,pcmat
!  double precision, dimension(:), allocatable :: logdetmin
!  double precision eps
!  read *, d,eps,jtrunc
!  dd=(d*(d-1))/2
!  allocate ( r(d,d), bnumv(d), permmat(d,d), pcmat(dd,d), logdetmin(d) )
!  do i=1,d
!    read *, r(i,:)
!    print '(9f10.6)', r(i,:)
!  end do
!  bnumv=0
!  permmat=0
!  pcmat=0.d0
!  logdetmin=0.d0
!  iprint=1
!  call gausstrvinenonuniq(d,dd,r,eps,jtrunc,iprint,bnumv,logdetmin,permmat,pcmat)
!  print *, bnumv
!  print *, logdetmin
!  do ell =1,(d-2)
!    print *, "best for truncation level ", ell 
!    print *, permmat(:,ell)
!    print '(8f9.5)', pcmat(:,ell)
!  end do
!  deallocate ( bnumv, permmat, pcmat, logdetmin )
!  stop 
!  end

! input  
!   d = dimension <=8, 
!   dd = d*(d-1)/2
!   r = dxd correlation matrix
!   eps = tolerance for within epsilon of the best, e.g., 1.e-7
!   jtrunc = truncation level for printing extra info (e.g., 2 or 3)
!   iprint = print flag for extra intermediate steps
!  Note: some results are printed and not returned
! output best ell-truncated vine for ell=1,...,d-2 based on min logdet, with
!   bnumv(ell) = index number of the vine array for best 
!   logdetmin(ell) = log determinant of best
!   permmat(,ell) = permutation of 1:d for best
!   pcmat(.ell)   = partial correlations for best, stored as d*(d-1)/2 vector
subroutine gausstrvinenonuniq(d,dd,r,eps,jtrunc,iprint,bnumv,logdetmin, permmat,pcmat)
  implicit none
  integer d,dd,jtrunc,iprint,bnumv(d),permmat(d,d)
  double precision eps,r(d,d),logdetmin(d),pcmat(dd,d)
  integer bnum,ncond,i,j,i1,i2
  integer mtc,even,dcase,dpow,ell,inew,idiff
  integer, dimension(:), allocatable :: aperm,best
  integer, dimension(:,:), allocatable :: A1,A2,tmpind
  double precision, dimension(:), allocatable :: logdet, outmin
  double precision, dimension(:,:), allocatable :: pcm
  double precision, dimension(:,:,:), allocatable :: pc3
  integer, dimension(:,:,:), allocatable :: vineind
  ! arrays for vine index of best truncated vines at each level
  integer, dimension(:), allocatable :: nonuniq
  ! keep track on number of non-unique

  ncond=2**d-4  
  dcase=((d-2)*(d-3))/2
  dpow=2**dcase
  allocate ( aperm(d), pcm(d,d), best(d-2), nonuniq(d-2) )
  allocate ( pc3(d,d,ncond), A1(d,d),A2(d,d), logdet(d-1),outmin(d-1) )
  allocate ( vineind(d,d,d-2), tmpind(d,d) )
  nonuniq=0; tmpind=0; vineind=0  

  call allpcor(d,ncond,r,pc3)
  outmin=0.d0
  do bnum=0,(dpow-1)
    inew=0
    call vnum2array(d,bnum,A1)
    mtc=0
    loop: do 
      call nexper(d,aperm,mtc,even)
      if (aperm(d-1)<aperm(d)) then
        call varrayperm(d,A1,aperm,A2) 
        best=0;
        call trvinelogdet(d,A2,r,pc3,logdet,pcm)
        do j=1,d-2
          if(logdet(j)<outmin(j)) then
            best(j)=1
            outmin(j)=logdet(j)
          else
            if(logdet(j)-outmin(j)<eps) then
              call vine2ind(d,j,A2,tmpind)
              idiff=maxval(tmpind-vineind(:,:,j))
              if(idiff>0) idiff=1
              nonuniq(j)=nonuniq(j)+idiff
              if (j==jtrunc) then
                call intpr("within eps for ",15,j,1)
                call intpr("bnum",4,bnum,1)
                call intpr("vine array",10,0,0)
                do i=1,d
                  call intpr("next row",8, A2(i,:),d)
                end do
                call intpr("partial corr",12,0,0)
                do i=1,j
                  call dblepr("next row",8, pcm(i,1:d),d)
                end do
              end if
            end if
          end if
        end do
        ! do some printing if smaller on one of first d-2 entries
        do ell=1,(d-2)
          if(best(ell)==1) then
            inew=1
            if(iprint==1 .and. ell==jtrunc) then
              call intpr("best",4,ell,1)
              call intpr("bnum",4,bnum,1)
              call intpr("perm",4,aperm,d)
              call dblepr("logdet",6,logdet,d-1)
              call intpr("partial corr",12,0,0)
              do i=1,(d-2)
                call dblepr("next row",8, pcm(i,1:d),d)
              end do
            end if
            nonuniq(ell)=0; tmpind=0
            call vine2ind(d,ell,A2,tmpind)
            vineind(:,:,ell)=tmpind
            ! save to variables to be returned
            bnumv(ell)=bnum
            logdetmin(ell)=logdet(ell)
            permmat(:,ell)=aperm
            i1=1; i2=d-1
            do i=1,(d-1)
              pcmat(i1:i2,ell)=pcm(i,(i+1):d)
              i1=i2+1; i2=i2+(d-i-1)
            end do
          end if
        end do
      end if  ! for aperm[d-1]<aperm[d]
      if (mtc==0) then
        exit loop
      end if
    end do loop
    if(inew==1 .and. iprint==1) then
      call intpr("above bnum=",11, bnum, 1)
      call intpr("A",1, 0,0)
      do i=1,d
        call intpr("next row",8, A1(i,:),d)
      end do
      call dblepr("min after this vine array",25,outmin,d-2)
      call intpr(" ",1,0,0)
    end if
  end do ! bnum
  if(iprint==1) then
    call intpr(" ",1,0,0)
  end if
  logdetmin(d-1)=logdet(d-1)  ! this is same for all vines
  call intpr("approx nonunique for j=1,...,d-2",32,nonuniq,d-2)
  ! this only gives an idea of the non-unique (some of the items
  !  listed could be the same, also only those "near best" found after 
  !  the best are counted)
  ! To actually find the number of distinct cases within eps of the
  ! best ell-truncated vine, save the output of this program
  ! and input it to a modification of this code.
  deallocate ( pc3, A1, A2, logdet, outmin, pcm )
  deallocate ( aperm, nonuniq, vineind )
  return
  end

! vine indexing: if (i1,i2) are the conditioned variables in tree ell 
!   for ell in 1,...,j; then ind(i1,i2)=ind(i2,i1)=ell
! inputs 
!   d = dimension
!   j = index for max tree to check
!   A = dxd vine array 
! output
!   ind = dxd array where ind(i1,i2)=ind(i2,i1)=ell if i1,i2 are in tree ell
subroutine vine2ind(d,j,A,ind)
  implicit none
  integer d,j,A(d,d),ind(d,d)
  integer i1,i2,k,ell
  ! ell is index for tree
  do ell=1,j
    do k=(ell+1),d
      i1=A(ell,k) ! edge is i1,i2 | A(1:(ell-1),k)
      i2=A(k,k)
      ind(i2,i1)=ell; ind(i1,i2)=ell
    end do
  end do
  !do ell=1,d
  !  print '(10i3)', A(ell,:)
  !end do
  !do ell=1,d
  !  print *, ind(ell,:)
  !end do
  return 
  end
