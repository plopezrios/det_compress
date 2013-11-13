MODULE tools
 !-----------------------------------------!
 ! Auxiliary tools for main program below. !
 !-----------------------------------------!
 IMPLICIT NONE
 PRIVATE
 PUBLIC wout,i2s,errstop,check_alloc,resize_pointer,iterate_indices_multi,&
  &qsort3_int_partial_preinit,iswap1

 ! Data types.
 INTEGER,PARAMETER :: dp=kind(1.d0)

 ! Output unit (standard output).
 INTEGER,PARAMETER :: o=6

 ! Dummy type for safe overloading of WOUT routines.
 TYPE arg_separator
  LOGICAL,POINTER :: dummy=>null()
 END TYPE arg_separator

 ! Interface for WOUT routine.
 INTERFACE wout
  MODULE PROCEDURE &
   &wout_NULL,wout_C,wout_C_II,wout_C_D
 END INTERFACE wout

 ! Interface for RESIZE_POINTER routine.
 INTERFACE resize_pointer
  MODULE PROCEDURE &
   &resize_pointer_int1,resize_pointer_int2,resize_pointer_int3,&
   &resize_pointer_bool1,resize_pointer_bool2,&
   &resize_pointer_dble1,resize_pointer_dble2
 END INTERFACE resize_pointer


CONTAINS


 ! Write-out (WOUT) routines for different argument types.


 SUBROUTINE wout_NULL(unused_arg_separator,fmt)
 CHARACTER(*),INTENT(in),OPTIONAL :: fmt
 TYPE(arg_separator),INTENT(in),OPTIONAL :: unused_arg_separator
 if(present(fmt))then
  write(o,fmt)
 else
  write(o,'(a)')'' ! force empty line, unlike write(o,*) with some compilers
 endif
 END SUBROUTINE wout_NULL


 SUBROUTINE wout_C(c,unused_arg_separator,fmt)
 CHARACTER(*),INTENT(in) :: c
 CHARACTER(*),INTENT(in),OPTIONAL :: fmt
 TYPE(arg_separator),INTENT(in),OPTIONAL :: unused_arg_separator
 if(present(fmt))then
  write(o,fmt)trim(c)
 else
  write(o,'(1x,a)')trim(c)
 endif
 END SUBROUTINE wout_C


 SUBROUTINE wout_C_II(c,ii,unused_arg_separator,fmt,rfmt,adjust,rsep)
 INTEGER,INTENT(in) :: ii(:)
 LOGICAL,INTENT(in),OPTIONAL :: adjust
 CHARACTER(*),INTENT(in) :: c
 CHARACTER(*),INTENT(in),OPTIONAL :: rfmt,fmt,rsep
 TYPE(arg_separator),INTENT(in),OPTIONAL :: unused_arg_separator
 INTEGER dsize,isize,jsize
 INTEGER,PARAMETER :: n_per_line=8
 LOGICAL adjust_nums
 CHARACTER(128) tmpr,tmpr2
 adjust_nums=.false. ; if(present(adjust))adjust_nums=adjust
 if(present(fmt))then
  write(o,fmt)c,ii
 else
  dsize=size(ii,1)
  do isize=1,dsize,n_per_line
   if(present(rfmt))then
    write(tmpr,rfmt)ii(isize)
    do jsize=1,min(n_per_line-1,dsize-isize)
     write(tmpr2,rfmt)ii(isize+jsize)
     if(adjust_nums)tmpr2=adjustl(tmpr2)
     if(present(rsep))then
      tmpr=trim(tmpr)//rsep//trim(tmpr2)
     else
      if(adjust_nums)then
       tmpr=trim(tmpr)//' '//trim(tmpr2)
      else
       tmpr=trim(tmpr)//trim(tmpr2)
      endif
     endif
    enddo
   else
    write(tmpr,*)ii(isize:min(dsize,isize+n_per_line-1))
    tmpr=adjustl(tmpr)
   endif
   if(isize==1)then
    write(o,'(1x,a)')c//trim(tmpr)
   else
    write(o,'(1x,a)')trim(tmpr)
   endif
  enddo
 endif
 END SUBROUTINE wout_C_II


 SUBROUTINE wout_C_D(c,d,unused_arg_separator,fmt,rfmt,adjust)
 REAL(dp),INTENT(in) :: d
 LOGICAL,INTENT(in),OPTIONAL :: adjust
 CHARACTER(*),INTENT(in) :: c
 CHARACTER(*),INTENT(in),OPTIONAL :: rfmt,fmt
 TYPE(arg_separator),INTENT(in),OPTIONAL :: unused_arg_separator
 LOGICAL adjust_nums
 CHARACTER(128) tmpr
 adjust_nums=.false. ; if(present(adjust))adjust_nums=adjust
 if(present(fmt))then
  write(o,fmt)c,d
 else
  if(present(rfmt))then
   write(tmpr,rfmt)d
   if(adjust_nums)tmpr=adjustl(tmpr)
  else
   write(tmpr,*)d
   tmpr=adjustl(tmpr)
  endif
  write(o,'(1x,a)')c//trim(tmpr)
 endif
 END SUBROUTINE wout_C_D


 ! Integer-to-character function


 CHARACTER(20) FUNCTION i2s(n)
!----------------------------!
! Convert integer to string. !
!----------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 INTEGER i,j,k
! Initialize.
 i2s=''
! Construct string (right to left).
 i=abs(n)
 do while(i>0)
  j=i/10
  k=i-j*10 ! = mod(i,10)
  if(k<0.or.k>9)then ! catch overflows
   i2s='[overflow]'
   return
  endif
  i2s=achar(ichar('0')+k)//trim(i2s)
  i=j
 enddo ! i>0
! Deal with numbers < 1.
 if(len_trim(i2s)==0)then
  i2s='0'
 elseif(n<0)then
  i2s='-'//trim(i2s)
 endif
 END FUNCTION i2s


 ! Error handling.


 SUBROUTINE errstop(routine,error)
 !-----------------------------!
 ! Stop with an error message. !
 !-----------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: routine,error
 call wout()
 call wout('ERROR : '//trim(routine))
 call wout(trim(error))
 call wout()
 stop
 END SUBROUTINE errstop


 SUBROUTINE check_alloc(ialloc,routine,symbol)
 !----------------------------------------!
 ! Error out on unsuccessful allocations. !
 !----------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: ialloc
 CHARACTER(*),INTENT(in) :: routine,symbol
 if(ialloc==0)return
 call wout()
 if(len_trim(symbol)==0)then
  call wout('ERROR : Allocation problem in '//trim(routine)//'.')
 else
  call wout('ERROR : Allocation problem in '//trim(routine)//' ('//&
   &trim(symbol)//').')
 endif
 call wout()
 stop
 END SUBROUTINE check_alloc


 ! Pointer-resizing utilities.


 SUBROUTINE resize_pointer_int1(dims,pt,init)
 !-------------------------------------------------------!
 ! Allocate or resize a first-rank integer pointer PT to !
 ! size DIMS, keeping any exisiting data untouched. New  !
 ! elements initialized to zero unless INIT specified.   !
 !-------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: dims(1)
 INTEGER,INTENT(in),OPTIONAL :: init
 INTEGER,POINTER :: pt(:),pt_new(:)=>null()
 INTEGER,PARAMETER :: init_default=0
 INTEGER old_dims(1),ialloc
 if(.not.associated(pt))then
  allocate(pt(dims(1)),stat=ialloc)
  call check_alloc(ialloc,'RESIZE_POINTER_INT1','pt')
  if(present(init))then
   pt=init
  else
   pt=init_default
  endif
  return
 endif
 old_dims=shape(pt)
 if(all(old_dims==dims))return
 allocate(pt_new(dims(1)),stat=ialloc)
 if(any(old_dims<dims))then
  if(present(init))then
   pt_new=init
  else
   pt_new=init_default
  endif
 endif
 pt_new(1:min(old_dims(1),dims(1)))=pt(1:min(old_dims(1),dims(1)))
 deallocate(pt)
 pt=>pt_new
 nullify(pt_new)
 END SUBROUTINE resize_pointer_int1


 SUBROUTINE resize_pointer_int2(dims,pt,init)
 !--------------------------------------------------------!
 ! Allocate or resize a second-rank integer pointer PT to !
 ! size DIMS, keeping any exisiting data untouched. New   !
 ! elements initialized to zero unless INIT specified.    !
 !--------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: dims(2)
 INTEGER,INTENT(in),OPTIONAL :: init
 INTEGER,POINTER :: pt(:,:),pt_new(:,:)=>null()
 INTEGER,PARAMETER :: init_default=0
 INTEGER old_dims(2),ialloc
 if(.not.associated(pt))then
  allocate(pt(dims(1),dims(2)),stat=ialloc)
  call check_alloc(ialloc,'RESIZE_POINTER_INT2','pt')
  if(present(init))then
   pt=init
  else
   pt=init_default
  endif
  return
 endif
 old_dims=shape(pt)
 if(all(old_dims==dims))return
 allocate(pt_new(dims(1),dims(2)),stat=ialloc)
 if(any(old_dims<dims))then
  if(present(init))then
   pt_new=init
  else
   pt_new=init_default
  endif
 endif
 pt_new(1:min(old_dims(1),dims(1)),1:min(old_dims(2),dims(2)))=&
  &pt(1:min(old_dims(1),dims(1)),1:min(old_dims(2),dims(2)))
 deallocate(pt)
 pt=>pt_new
 nullify(pt_new)
 END SUBROUTINE resize_pointer_int2


 SUBROUTINE resize_pointer_int3(dims,pt,init)
 !-------------------------------------------------------!
 ! Allocate or resize a third-rank integer pointer PT to !
 ! size DIMS, keeping any existing data untouched. New   !
 ! elements initialized to zero unless INIT specified.   !
 !-------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: dims(3)
 INTEGER,INTENT(in),OPTIONAL :: init
 INTEGER,POINTER :: pt(:,:,:),pt_new(:,:,:)=>null()
 INTEGER,PARAMETER :: init_default=0
 INTEGER old_dims(3),ialloc
 if(.not.associated(pt))then
  allocate(pt(dims(1),dims(2),dims(3)),stat=ialloc)
  call check_alloc(ialloc,'RESIZE_POINTER_INT3','pt')
  if(present(init))then
   pt=init
  else
   pt=init_default
  endif
  return
 endif
 old_dims=shape(pt)
 if(all(old_dims==dims))return
 allocate(pt_new(dims(1),dims(2),dims(3)),stat=ialloc)
 if(any(old_dims<dims))then
  if(present(init))then
   pt_new=init
  else
   pt_new=init_default
  endif
 endif
 pt_new(1:min(old_dims(1),dims(1)),1:min(old_dims(2),dims(2)),&
  &1:min(old_dims(3),dims(3)))=pt(1:min(old_dims(1),dims(1)),&
  &1:min(old_dims(2),dims(2)),1:min(old_dims(3),dims(3)))
 deallocate(pt)
 pt=>pt_new
 nullify(pt_new)
 END SUBROUTINE resize_pointer_int3


 SUBROUTINE resize_pointer_bool1(dims,pt,init)
 !--------------------------------------------------------!
 ! Allocate or resize a first-rank logical pointer PT to  !
 ! size DIMS, keeping any existing data untouched. New    !
 ! elements initialized to .false. unless INIT specified. !
 !--------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: dims(1)
 LOGICAL,INTENT(in),OPTIONAL :: init
 LOGICAL,POINTER :: pt(:),pt_new(:)=>null()
 LOGICAL,PARAMETER :: init_default=.false.
 INTEGER old_dims(1),ialloc
 if(.not.associated(pt))then
  allocate(pt(dims(1)),stat=ialloc)
  call check_alloc(ialloc,'RESIZE_POINTER_BOOL1','pt')
  if(present(init))then
   pt=init
  else
   pt=init_default
  endif
  return
 endif
 old_dims=shape(pt)
 if(all(old_dims==dims))return
 allocate(pt_new(dims(1)),stat=ialloc)
 if(any(old_dims<dims))then
  if(present(init))then
   pt_new=init
  else
   pt_new=init_default
  endif
 endif
 pt_new(1:min(old_dims(1),dims(1)))=pt(1:min(old_dims(1),dims(1)))
 deallocate(pt)
 pt=>pt_new
 nullify(pt_new)
 END SUBROUTINE resize_pointer_bool1


 SUBROUTINE resize_pointer_bool2(dims,pt,init)
 !--------------------------------------------------------!
 ! Allocate or resize a second-rank logical pointer PT to !
 ! size DIMS, keeping any exisiting data untouched. New   !
 ! elements initialized to .false. unless INIT specified. !
 !--------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: dims(2)
 LOGICAL,INTENT(in),OPTIONAL :: init
 LOGICAL,POINTER :: pt(:,:),pt_new(:,:)=>null()
 LOGICAL,PARAMETER :: init_default=.false.
 INTEGER old_dims(2),ialloc
 if(.not.associated(pt))then
  allocate(pt(dims(1),dims(2)),stat=ialloc)
  call check_alloc(ialloc,'RESIZE_POINTER_BOOL2','pt')
  if(present(init))then
   pt=init
  else
   pt=init_default
  endif
  return
 endif
 old_dims=shape(pt)
 if(all(old_dims==dims))return
 allocate(pt_new(dims(1),dims(2)),stat=ialloc)
 if(any(old_dims<dims))then
  if(present(init))then
   pt_new=init
  else
   pt_new=init_default
  endif
 endif
 pt_new(1:min(old_dims(1),dims(1)),1:min(old_dims(2),dims(2)))=&
  &pt(1:min(old_dims(1),dims(1)),1:min(old_dims(2),dims(2)))
 deallocate(pt)
 pt=>pt_new
 nullify(pt_new)
 END SUBROUTINE resize_pointer_bool2


 SUBROUTINE resize_pointer_dble1(dims,pt,init)
 !------------------------------------------------------!
 ! Allocate or resize a first-rank dble pointer PT to   !
 ! size DIMS, keeping any existing data untouched. New  !
 ! elements initialized to zero unless INIT specified.  !
 !------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: dims(1)
 REAL(dp),INTENT(in),OPTIONAL :: init
 REAL(dp),POINTER :: pt(:),pt_new(:)=>null()
 REAL(dp),PARAMETER :: init_default=0.d0
 INTEGER old_dims(1),ialloc
 if(.not.associated(pt))then
  allocate(pt(dims(1)),stat=ialloc)
  call check_alloc(ialloc,'RESIZE_POINTER_DBLE1','pt')
  if(present(init))then
   pt=init
  else
   pt=init_default
  endif
  return
 endif
 old_dims=shape(pt)
 if(all(old_dims==dims))return
 allocate(pt_new(dims(1)),stat=ialloc)
 if(any(old_dims<dims))then
  if(present(init))then
   pt_new=init
  else
   pt_new=init_default
  endif
 endif
 pt_new(1:min(old_dims(1),dims(1)))=pt(1:min(old_dims(1),dims(1)))
 deallocate(pt)
 pt=>pt_new
 nullify(pt_new)
 END SUBROUTINE resize_pointer_dble1


 SUBROUTINE resize_pointer_dble2(dims,pt,init)
 !------------------------------------------------------!
 ! Allocate or resize a second-rank dble pointer PT to  !
 ! size DIMS, keeping any exisiting data untouched. New !
 ! elements initialized to zero unless INIT specified.  !
 !------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: dims(2)
 REAL(dp),INTENT(in),OPTIONAL :: init
 REAL(dp),POINTER :: pt(:,:),pt_new(:,:)=>null()
 REAL(dp),PARAMETER :: init_default=0.d0
 INTEGER old_dims(2),ialloc
 if(.not.associated(pt))then
  allocate(pt(dims(1),dims(2)),stat=ialloc)
  call check_alloc(ialloc,'RESIZE_POINTER_DBLE2','pt')
  if(present(init))then
   pt=init
  else
   pt=init_default
  endif
  return
 endif
 old_dims=shape(pt)
 if(all(old_dims==dims))return
 allocate(pt_new(dims(1),dims(2)),stat=ialloc)
 if(any(old_dims<dims))then
  if(present(init))then
   pt_new=init
  else
   pt_new=init_default
  endif
 endif
 pt_new(1:min(old_dims(1),dims(1)),1:min(old_dims(2),dims(2)))=&
  &pt(1:min(old_dims(1),dims(1)),1:min(old_dims(2),dims(2)))
 deallocate(pt)
 pt=>pt_new
 nullify(pt_new)
 END SUBROUTINE resize_pointer_dble2


 ! Index iterator.


 LOGICAL FUNCTION iterate_indices_multi(n,order,indices,is_first)
 !--------------------------------------------------------------!
 ! Provided a set of N integer indices, INDICES(1:N), with a    !
 ! range of 1..ORDER(I), generate the next set.  E.g.:          !
 !  (1,1,3) -> (1,1,4) -> (1,2,1) -> (1,2,2) ...                !
 ! The function returns a value of .false. if there are no more !
 ! index sets to loop over. If IS_FIRST=.true. on input, the    !
 ! first index set is automatically generated, and IS_FIRST     !
 ! is returned as .false. (this will cause the function to      !
 ! return .true. even if N=0).                                  !
 !--------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n,order(n)
 INTEGER,INTENT(inout) :: indices(n)
 LOGICAL,INTENT(inout) :: is_first
 INTEGER i,k
 LOGICAL is_invalid
 iterate_indices_multi=.true.
 is_invalid=n<1.or.any(order<1)
 if(is_first)then
  is_first=.false.
  ! Initialize indices to 1,1...1
  if(.not.is_invalid)indices(1:n)=1
  return
 elseif(is_invalid)then
  iterate_indices_multi=.false.
  return
 endif
 ! Increase counters
 do i=n,1,-1
  k=indices(i)+1
  if(k>order(i))cycle
  indices(i)=k
  indices(i+1:n)=1
  return
 enddo
 iterate_indices_multi=.false.
 END FUNCTION iterate_indices_multi


 ! Swap routines


 SUBROUTINE iswap1(i,j)
 !------------------------!
 ! Swap integers I and J. !
 !------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(inout) :: i,j
 INTEGER k
 k=i ; i=j ; j=k
 END SUBROUTINE iswap1


 ! Quicksort routine.


 SUBROUTINE qsort3_int_partial_preinit(n,i1,i2,x,indx)
 !--------------------------------------------------------!
 ! 3-way QuickSort (=QuickSort3) for a part of an integer !
 ! vector array with pre-initialized INDX.  I.e., we are  !
 ! rearranging INDX(I1:I2) so that X(INDX(I1:I2)) is      !
 ! sorted in ascending order.                             !
 !                                                        !
 ! QuickSort3 is the best method for integers because it  !
 ! takes repetitions into account.  Repetitions are rare  !
 ! or non-exisiting in many real cases, but are to be     !
 ! expected in integer problems.  There is an additional  !
 ! overhead due to checking for equalities, so 2-way      !
 ! QuickSort (=QuickSort) should be used when no          !
 ! repetitions are expected.                              !
 !                                                        !
 ! As with QuickSort above, we resort to Insertion sort   !
 ! when all unsorted sub-sequences are of size less than  !
 ! some threshold (M) because the lower overhead of       !
 ! insertion is favourable in small sequences.            !
 !                                                        !
 ! PLR 02.2012                                            !
 !--------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n,i1,i2,x(n)
 INTEGER,INTENT(inout) :: indx(n)
 INTEGER,PARAMETER :: m=30,stack_size=128
 INTEGER lstk(stack_size),rstk(stack_size),istk
 INTEGER l,r,i,j,k,p,q,xr

 if(i2-i1>=m)then

  istk=0 ; l=i1 ; r=i2
  do
   i=l-1 ; j=r ; xr=x(indx(r))
   p=l-1 ; q=r
   do
    i=i+1
    do while(x(indx(i))<xr)
     i=i+1
    enddo
    j=j-1
    do while(x(indx(j))>xr)
     if(j==l)exit
     j=j-1
    enddo
    if(i>=j)exit
    call iswap1(indx(i),indx(j))
    if(x(indx(i))==xr)then
     p=p+1 ; call iswap1(indx(p),indx(i))
    endif
    if(xr==x(indx(j)))then
     q=q-1 ; call iswap1(indx(j),indx(q))
    endif
   enddo
   call iswap1(indx(i),indx(r))
   j=i+1 ; i=i-1
   do k=l,p-1
    call iswap1(indx(k),indx(i)) ; i=i-1
   enddo
   do k=r-1,q+1,-1
    call iswap1(indx(j),indx(k)) ; j=j+1
   enddo
   if(r-j>=i-l.and.i-l>=m)then
    istk=istk+1
    if(istk>stack_size)call errstop('QSORT3_INT_PARTIAL_PREINIT',&
     &'Stack size too small <1>.')
    lstk(istk)=j
    rstk(istk)=r
    r=i
   elseif(i-l>r-j.and.r-j>=m)then
    istk=istk+1
    if(istk>stack_size)call errstop('QSORT3_INT_PARTIAL_PREINIT',&
     &'Stack size too small <2>.')
    lstk(istk)=l
    rstk(istk)=i
    l=j
   elseif(i-l>=m)then
    r=i
   elseif(r-j>=m)then
    l=j
   else
    if(istk<1)exit
    l=lstk(istk)
    r=rstk(istk)
    istk=istk-1
   endif
  enddo

 endif

 do i=i1+1,i2
  if(x(indx(i-1))<=x(indx(i)))cycle
  call iswap1(indx(i-1),indx(i))
  do j=i-1,i1+1,-1
   if(x(indx(j-1))<=x(indx(j)))exit
   call iswap1(indx(j-1),indx(j))
  enddo ! j
 enddo ! i

 END SUBROUTINE qsort3_int_partial_preinit


END MODULE tools
