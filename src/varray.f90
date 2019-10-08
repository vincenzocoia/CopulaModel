! conversion of some vine array functions in varray.R  to f90
! gfortran -o varray varray.f90
! gfortran -o varray varray.f90 pcor.f90

!program mainvarr
!  implicit none
!  integer d,bnum,ii,no1,first1,prev,ncond,i,j,ell1,indx
!  integer mtc,even,dfact,dcase,dpow,ell,best1,best2
!  integer, dimension(:), allocatable :: aperm,best
!  integer, dimension(:,:), allocatable :: A1,A2
!  double precision, dimension(:), allocatable :: logdet, outmax
!  double precision, dimension(:,:), allocatable :: r,pcm
!  double precision, dimension(:,:,:), allocatable :: pc3

  !d=5;
!  read *,d !,dfact
!  ncond=2**d-4  
!  dcase=((d-2)*(d-3))/2
!  dpow=2**dcase
!  allocate ( aperm(d), pcm(d,d), best(d-2) )
!  allocate ( pc3(d,d,ncond), r(d,d), A1(d,d),A2(d,d), logdet(d-1),outmax(d-1) )
!  do i=1,d
    !do j=1,d
    !   r(i,j)=0.5**(abs(i-j))
    !end do
!    read *, r(i,:)
!    print '(7f10.6)', r(i,:)
!  end do
  
!  call allpcor(d,ncond,r,pc3)
  !do i=1,d
  !  print '(5f10.6)', pc3(i,:,1)
  !end do
  !read *, bnum
!  outmax=0.d0
!  do bnum=0,(dpow-1)
  !do bnum=0,1
!    call vnum2array(d,bnum,A1)
!    print *
!    print *,"bnum=", bnum
!    do i=1,d
!      print '(7i4)', A1(i,:)
!    end do
    !outmax=0.d0 ! move this outside loop later
!    mtc=0
!    loop: do 
    ! do ii=1,dfact
    ! do ii=1,5
!      call nexper(d,aperm,mtc,even)
      !print *,aperm, even
!      call varrayperm(d,A1,aperm,A2) 
!      best=0;
!      call trvinelogdet(d,A2,r,pc3,logdet,pcm)
!      do j=1,d-2
!        if(logdet(j)<outmax(j)) then
!          best(j)=1
!          outmax(j)=logdet(j)
!        end if
!      end do
      !print '(5i3,4f10.6)', aperm,logdet
      ! do some printing if smaller on one of first d-2 entries
!      do ell=1,(d-3)
!        if(best(ell)==1) then
!          print '(a4,i1,1x,7i3)', "best",ell,aperm
!          print '(6f10.6)', logdet
!          do i=1,(d-3)
!            print '(i3,7f10.6)', i,(pcm(i,j),j=(i+1),d)
!          end do
!        end if
!      end do
      !if(best1==1) then
      !  print '(5i3,4f10.6)', aperm,logdet
      !  print '(a5,4f10.6,2x,3f10.6)', "best1",(pcm(1,j),j=2,d),(pcm(2,j),j=3,d)
      !end if
      !if(best2==1) then
      !  print '(5i3,4f10.6)', aperm,logdet
      !  print '(a5, 7f10.6)', "best2",(pcm(1,j),j=2,d),(pcm(2,j),j=3,d)
      !end if
!      if (mtc==0) then
!        exit loop
!      end if
!    end do loop
    !print *, "min for this vine array"
!    print *, "min after this vine array"
!    print '(7f10.6)', outmax
!  end do ! bnum
!  print *
!  print *,"============================================================"

!  deallocate ( pc3, r, A1, A2, logdet, outmax, pcm )
!  deallocate ( aperm )
!  stop
!  end

! vine array after permutation of indices
! inputs
!   d = dimension
!   A1 = dxd vine array
!   aperm = permutation of 1:d
! output
!   A2 = dxd vine array (after permutation of indices)
subroutine varrayperm(d,A1,aperm,A2) 
  implicit none
  integer d,A1(d,d),aperm(d),A2(d,d)
  integer i,j
  A2=A1
  do i = 1,d
   do j = i,d 
     A2(i,j)=aperm(A1(i,j))
   end do
  end do
  return 
  end

! next permutation of 1:d
! routine for generating the permutations of 1,...,n in systematic
! order, a has length n, n>=2
  subroutine nexper(n,a,mtc,even)
  implicit none
  integer n,a(n),mtc,even
  integer s,d,nm3,ia,i,i1,j,l,m
  !integer, dimension(:), allocatable :: itaken,b
  !allocate ( itaken(i), b(i) )

  ! initialize to avoid warning
  l=0
  ! interchange first two lines so that nm3 is always defined      
  nm3=n-3
  if(mtc==0) then
    a=(/(i, i=1,n)/)
    mtc=1
    even=1
    if(n.eq.1) mtc=0;
    return 
  end if
  !else 
    if(n==1) then
      a(1)=0; mtc=0; return 
    endif
    if(even==1) then
      ia=a(1); a(1)=a(2); a(2)=ia
      even=0
      !  goto 6
      if(a(n).ne.1 .or. a(1).ne.2+mod(n,2)) return
      if(n.le.3) then
        mtc=0; return
      end if
      do i=1,nm3
        if(a(i+1).ne.a(i)+1) return
      end do 
      mtc=0
    else ! even==0
      s=0
      do i1=2,n  ! error in this loop
        ia=a(i1)
        i=i1-1
        d=0
        do j=1,i
          if(a(j)>ia) d=d+1
        end do
        s=d+s
        if(d.ne.i*mod(s,2)) then
           m=mod(s+1,2)*(n+1)
           do j=1,i
             if(isign(1,a(j)-ia).ne.isign(1,a(j)-m)) then
             !!  isign(1,x)=1 if x>=0, =-1 if x<0
               m=a(j)
               l=j
             endif
           end do
           a(l)=ia
           a(i1)=m
           even=1
           return ! this is needed here
        end if
      end do
      a(1)=0; mtc=0; return 
    end if
  !end if
  return
  end


! inputs
!   d = dimension
!   b0 = binary vector of length i-3 
!   A1 = dxd vine array
!   i = column number for converting the binary representation (4<=i<=d)
! output
!   A1(,i) is replaced based on b0
subroutine vstepb(d,b0,A1,i)
  implicit none
  integer d,b0(d),A1(d,d),i
  integer ac,tem,k,ell
  integer, dimension(:), allocatable :: itaken,b
  allocate ( itaken(i), b(i) )
  itaken=0
  itaken(i)=1; itaken(i-1)=1
  b=1
  b(2:(i-2))=b0(1:(i-3))
  ac=i-2  ! active column
  A1(i,i)=i
  A1(i-1,i)=i-1
  do k = (i-2),1,-1
    if(b(k)==1) then
      tem=A1(ac,ac); itaken(tem)=1
      A1(k,i)=tem
      if(k>1) then
        ! { ac=max((1:i)[itaken==0]) }
        do ell=1,i
          if(itaken(ell)==0) ac=ell
        end do
      end if
    else
      tem=A1(k-1,ac); A1(k,i)=tem; itaken(tem)=1  
    end if
  end do
  !print *,A1(:,i)
  deallocate ( b,itaken )
  !A1[,i]
  return
  end

! decimal to binary vector
! inputs
!   n = positive integer
!   ii = integer between 0 and 2^n-1
! output
!   jj = binary vector of size n with binary representation of ii
subroutine d2b(n, ii, jj)
  implicit none
  integer n,ii, jj(n)
  integer tt,i
  tt=ii
  jj=0
  do i = n,1,-1
    jj(i)=mod(tt,2); tt=int(tt/2); 
  end do
  return 
  end

! Function with b vector and calls to vstepb
! inputs
!   d = dimension
!   bnum = integer between 0 and 2^((d-2)*(d-3)/2)-1 ;
!   2^((d-2)*(d-3)/2) is the number of d-dimensional vine arrays
!     in natural order.
!  bnum=0 is the D-vine, bnum=2^((d-2)*(d-3)/2)-1 is the C-vine
!  To get bnum from a binary matrix representation bmat,
!    bvec as a vector starts increments from the last position,
!    which is bmat[d-2,d];
!  the (d-2)*(d-3)/2 positions in bmat: [2,4], [3,4],[3,5],...,[2,d],...,[d-2,d]
! output:
!   A = dxd vine array in natural order form
subroutine vnum2array(d,bnum,A)
  implicit none
  integer d,bnum, A(d,d)
  integer dcase,dpow,i,ii
  integer, dimension(:), allocatable :: bvec,b0  
  integer, dimension(:,:), allocatable :: b

  dcase=(d-2)*(d-3)/2
  allocate ( b(d,d), b0(d), bvec(dcase) )
  dpow=2**dcase
  if(bnum<0 .or. bnum>=dpow) bnum=0
  call d2b(dcase,bnum,bvec)
  b=0; A=0
  b(1,:)=1; A(1,1)=1
  do i=2,d
    b(i,i)=1; A(i,i)=i; A(i-1,i)=i-1
  end do
  A(1,3)=1
  do i = 3,d 
    b(i-1,i)=1
  end do
  ii=0
  do i = 4,d 
    b(2:(i-2),i)=bvec((ii+1):(ii+i-3)); ii=ii+i-3 
  end do
  !do i=1,d
  !  print *,b(i,:)
  !end do
  
  ! column 4
  if(b(2,4)==1) then ! C-vine for first 4 columns 
    A(1,4)=1; A(2,4)=2 
  else 
    A(1,4)=2; A(2,4)=1 ! D-vine for first 4 columns
  end if
  ! columns 5 and higher
  do i = 5,d 
    b0(1:(i-3))=b(2:(i-2),i)  ! length i-3
    call vstepb(d,b0,A,i) ! changes to A in column i
  end do
  !print *
  !do i=1,d
  !  print *, A(i,:)
  !end do
  deallocate ( b, b0, bvec )
  return
  end

