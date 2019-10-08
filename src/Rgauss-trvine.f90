! best Gaussian truncated vines with 1,2,3,...d-2 trees
! version with dblpr, intpr print statements for R calls

! input  
!   d = dimension <=8, 
!   dd = d*(d-1)/2
!   r = dxd correlation matrix
!   iprint = print flag for intermediate steps
! output: best ell-truncated vine for ell=1,...,d-2 based on min logdet, with
!   bnumv(ell) = index number of the vine array for best 
!   logdetmin(ell) = log determinant of best
!   permmat(,ell) = permutation of 1:d for best
!   pcmat(.ell)   = partial correlations for best, stored as d*(d-1)/2 vector
subroutine gausstrvine(d,dd,r,iprint,bnumv,logdetmin,permmat,pcmat)
  implicit none
  integer d,dd,iprint,bnumv(d),permmat(d,d)
  double precision r(d,d),logdetmin(d),pcmat(dd,d)
  integer bnum,ncond,i,j,i1,i2
  integer mtc,even,dcase,dpow,ell,inew
  integer, dimension(:), allocatable :: aperm,best
  integer, dimension(:,:), allocatable :: A1,A2
  double precision, dimension(:), allocatable :: logdet, outmin
  double precision, dimension(:,:), allocatable :: pcm
  double precision, dimension(:,:,:), allocatable :: pc3

  ncond=2**d-4  
  dcase=((d-2)*(d-3))/2
  dpow=2**dcase
  allocate ( aperm(d), pcm(d,d), best(d-2) )
  allocate ( pc3(d,d,ncond), A1(d,d),A2(d,d), logdet(d-1),outmin(d-1) )
  
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
          end if
        end do
        ! do some printing if smaller on one of first d-2 entries
        do ell=1,(d-2)
          if(best(ell)==1) then
            inew=1
            if(iprint==1) then
              call intpr("best",4,ell,1)
              call intpr("perm",4,aperm,d)
              call dblepr("logdet",6,logdet,d-1)
              call intpr("partial corr",12,0,0)
              do i=1,(d-2)
                call dblepr("next row",8, pcm(i,1:d),d)
              end do
            end if
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
      end if ! for aperm[d-1]<aperm[d]
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

  deallocate ( pc3, A1, A2, logdet, outmin, pcm )
  deallocate ( aperm )
  return
  end

