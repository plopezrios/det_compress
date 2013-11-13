MODULE tools
 !-----------------------------------------!
 ! Auxiliary tools for main program below. !
 !-----------------------------------------!
 IMPLICIT NONE
 PRIVATE
 PUBLIC timer,timer_start,timer_end,wout,i2s,errstop,check_alloc,&
  &resize_pointer,iterate_indices_multi,qsort3_int_partial_preinit,&
  &iswap1

 ! Data types.
 INTEGER,PARAMETER :: dp=kind(1.d0)
 INTEGER,PARAMETER :: sp=kind(1.0)
 INTEGER,PARAMETER :: i64=selected_int_kind(15)

 ! Output unit (standard output).
 INTEGER,PARAMETER :: o=6

 ! Timer type.
 TYPE timed_routine
  ! This derived type contains the timing information of a particular
  ! routine in the hierarchy, as well as pointers to its parent, its
  ! first child and its next sibling.
  REAL(sp) time(3),start(3),rtime(3)
  CHARACTER(30) name
  LOGICAL error
  INTEGER(i64) hits
  TYPE(timed_routine),POINTER :: first_sublevel,next,parent
 END TYPE timed_routine

 ! Timer variables
 INTEGER :: ilevel_current=0
 REAL(sp) start_time(3),start_wtime
 INTEGER,PARAMETER :: hard_max_level=20
 TYPE(timed_routine),POINTER :: level_zero,pt_current

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


 ! Timer routines.


 SUBROUTINE timer(label,activate)
!-------------------------------------------------------------!
! Keep track of the time spent by routine LABEL, starting the !
! timer if ACTIVATE or stopping it if not. The optional flag  !
! COLLAPSE will prevent routines nested inside LABEL from     !
! being shown in the final timer report.                      !
! PLR 01.2007                                                 !
!-------------------------------------------------------------!
 IMPLICIT NONE
 LOGICAL,INTENT(in) :: activate
 CHARACTER(*),INTENT(in) :: label
 REAL(sp) current_time(3),time_diff(3)
 LOGICAL err
 TYPE(timed_routine),POINTER :: pt_search

 if(activate)then

  err=.false.
  if(ilevel_current>=hard_max_level)return
  ilevel_current=ilevel_current+1
  ! Determine whether this is a new routine or an already-known one.
  pt_search=>pt_current%first_sublevel
  if(associated(pt_search,pt_current))then
   ! Current level has no children: create first child.
   call create_timed_routine(pt_current,pt_search,label,.true.,err)
  else
   ! Search for label in tree after timed routine.
   do
    if(trim(pt_search%name)==label)exit ! found
    if(associated(pt_search%next,pt_current))then
     ! End of children from current level: add a new child.
     call create_timed_routine(pt_current,pt_search,label,.false.,err)
     exit
    endif
    pt_search=>pt_search%next ! next child
   enddo
  endif
  if(err)return ! ignore allocation errors
  ! Point pt_current at the new routine, and update its contents.
  pt_current=>pt_search
  pt_current%hits=pt_current%hits+1_i64
  call cputime(current_time)
  pt_current%start(:)=current_time(:)

 else ! deactivate

  call cputime(current_time)
  if(trim(pt_current%name)==label)then
   ! Deactivate current timer.
   time_diff(:)=current_time(:)-pt_current%start(:)
   pt_current%time(:)=pt_current%time(:)+time_diff(:)
   pt_current=>pt_current%parent
   ilevel_current=ilevel_current-1
  else
   ! This is an error. Try finding timer label in parent levels.
   pt_search=>pt_current
   do
    pt_search=>pt_search%parent
    if(associated(pt_search,level_zero))return ! timer not found - no action
    if(trim(pt_search%name)==label)exit ! timer found
   enddo
   ! Timer found, deactivate all timers in the middle and flag the error.
   do
    time_diff(:)=current_time(:)-pt_current%start(:)
    pt_current%time(:)=pt_current%time(:)+time_diff(:)
    if(associated(pt_current,pt_search))exit
    pt_current%error=.true.
    pt_current=>pt_current%parent
    ilevel_current=ilevel_current-1
   enddo
  endif ! correct deactivation or not

 endif ! activate or deactivate timer

 END SUBROUTINE timer


 SUBROUTINE create_timed_routine(my_parent,new_level,label,first,err)
 !--------------------------------------------------------------!
 ! Create a timed routine named LABEL as the last child of      !
 ! MY_PARENT, after sibling NEW_LEVEL (if not NULL), and return !
 ! a pointer to the newly created item in NEW_LEVEL on output.  !
 !--------------------------------------------------------------!
 IMPLICIT NONE
 TYPE(timed_routine),POINTER :: my_parent,new_level
 LOGICAL,INTENT(in) :: first
 LOGICAL,INTENT(out) :: err
 CHARACTER(*),INTENT(in) :: label
 INTEGER ialloc
 TYPE(timed_routine),POINTER :: pt=>null()
 err=.false.
 allocate(pt,stat=ialloc)
 if(ialloc/=0)then
  err=.true.
  return
 endif
 ! Links to the new item
 if(first)then
  my_parent%first_sublevel=>pt
 else
  new_level%next=>pt
 endif
 new_level=>pt
 nullify(pt)
 new_level%name=label
 new_level%time=0.0
 new_level%rtime=0.0
 new_level%start=0.0
 new_level%hits=0
 new_level%error=.false.
 ! Links from new item
 new_level%first_sublevel=>new_level
 new_level%next=>my_parent
 new_level%parent=>my_parent
 END SUBROUTINE create_timed_routine


 SUBROUTINE timer_start
 !------------------!
 ! Start the timers !
 !------------------!
 IMPLICIT NONE
 allocate(level_zero)
 pt_current=>level_zero
 level_zero%name='LEVEL_ZERO'
 level_zero%time=0.0
 level_zero%rtime=0.0
 level_zero%start=0.0
 level_zero%hits=1
 level_zero%error=.false.
 level_zero%first_sublevel=>level_zero
 level_zero%next=>level_zero
 level_zero%parent=>level_zero
 call cputime(start_time)
 start_wtime=walltime()
 END SUBROUTINE timer_start


 SUBROUTINE timer_end
!------------------------------------------!
! Stop the timers for the whole program.   !
!------------------------------------------!
 IMPLICIT NONE
 INTEGER level
 INTEGER(i64) hits
 REAL(sp) total_time,total_wtime,cpu_ratio_branch,cpu_ratio_tot,tv1(3),&
  &system_ratio
 LOGICAL errors
 CHARACTER(40) label
 CHARACTER(80) line
 TYPE(timed_routine),POINTER :: pt_current,pt_temp

 ! Get final time and report totals.
 call cputime(tv1)
 total_time=tv1(3)-start_time(3)
 total_wtime=walltime()-start_wtime
 call wout()
 call wout('Total CPU time : ',dble(total_time),rfmt='(f13.4)')
 call wout('Total real time: ',dble(total_wtime),rfmt='(f13.4)')
 call wout()

 ! Print header for timer table.
 level=0 ; pt_current=>level_zero ; errors=.false.
 pt_current%time(3)=total_time
 call wout('Timing information:')
 call wout('----------------------------------CPU---Sys/CPU--CPU/brn--&
  &CPU/tot--Hits--')

 ! Print all timers.
 mainprint: do
  ! Find next label to print.
  if(associated(pt_current%first_sublevel,pt_current))then
   ! No deeper level in this branch.
   if(associated(pt_current%next,pt_current%parent))then
    ! No next same-level from parent. Print unaccounted-for time for parent,
    ! then go up the 'parent' hierarchy until a 'next' pointer can be
    ! followed.
    do
     if(associated(pt_current%parent,pt_current))exit mainprint
     pt_current=>pt_current%parent ; level=level-1
     if(pt_current%rtime(3)>0.)then
      label=repeat('-',level)//'[rest]'
      system_ratio=0. ; cpu_ratio_tot=0. ; cpu_ratio_branch=0.
      if(pt_current%rtime(3)>0.)system_ratio=min(1.,&
       &max(0.,pt_current%rtime(2)/pt_current%rtime(3)))
      if(total_time>0.)cpu_ratio_tot=min(1.,max(0.,pt_current%rtime(3)/&
       &total_time))
      if(pt_current%time(3)>0.)cpu_ratio_branch=min(1.,max(0.,&
       &pt_current%rtime(3)/pt_current%time(3)))
      call format_timer_output(label,system_ratio,max(0.0,&
       &pt_current%rtime(3)),cpu_ratio_branch,cpu_ratio_tot,-1_i64,line)
      call wout(trim(line),fmt='(a)')
     endif
     if(associated(pt_current%next,pt_current%parent))cycle
     pt_current=>pt_current%next
     exit
    enddo
   else
    pt_current=>pt_current%next
   endif
  else
   ! Gather total times for children.
   tv1(:)=0. ; pt_temp=>pt_current%first_sublevel
   do
    tv1=tv1+pt_temp%time
    if(associated(pt_temp%next,pt_temp%parent))exit
    pt_temp=>pt_temp%next
   enddo
   pt_current%rtime=pt_current%time-tv1
   ! Move on to next level.
   pt_current=>pt_current%first_sublevel ; level=level+1
  endif
  system_ratio=0. ; cpu_ratio_tot=0. ; cpu_ratio_branch=0.
  if(pt_current%time(3)>0.)system_ratio=min(1.,max(0.,pt_current%time(2)/&
   &pt_current%time(3)))
  if(total_time>0.)cpu_ratio_tot=min(1.,max(0.,pt_current%time(3)/total_time))
  pt_temp=>pt_current%parent
  if(pt_temp%time(3)>0.)cpu_ratio_branch=min(1.,max(0.,pt_current%time(3)/&
   &pt_temp%time(3)))
  label=repeat('-',level-1)//trim(pt_current%name)
  if(pt_current%error)then
   label=trim(label)//'*' ; errors=.true.
  endif
  hits=pt_current%hits
  if(total_time>0.)cpu_ratio_tot=min(1.,max(0.,pt_current%time(3)/total_time))
  call format_timer_output(label,system_ratio,max(0.0,pt_current%time(3)),&
   &cpu_ratio_branch,cpu_ratio_tot,hits,line)
  call wout(trim(line),fmt='(a)')
 enddo mainprint

 ! Print footer for timer table.
 call wout(repeat('-',73))
 call wout('      CPU: CPU time (seconds)     Sys/CPU: system-to-CPU &
  &time ratio')
 call wout('   CPU/brn: CPU time (% of branch)      CPU/tot: CPU time &
  &(% of total)')
 if(errors)call wout('             *: routines without correct timer &
  &deactivation.')
 call wout(repeat('=',73))

 END SUBROUTINE timer_end


 SUBROUTINE format_timer_output(label,s_ratio,c_time,c_ratio_branch,&
  &c_ratio_tot,hits,line)
 !-----------------------------------!
 ! Produce a formatted timer string. !
 !-----------------------------------!
 IMPLICIT NONE
 INTEGER(i64),INTENT(in) :: hits
 REAL(sp),INTENT(in) :: s_ratio,c_time,c_ratio_branch,c_ratio_tot
 CHARACTER(40),INTENT(in) :: label
 CHARACTER(80),INTENT(out) :: line
 INTEGER temp_power,hits10
 INTEGER(i64) temp_hits
 REAL(sp) rhits
 CHARACTER(40) hits_string,suffix
 if(hits>=0_i64)then
  ! Divide hits by 10^(multiple of 3) to display in human-readable form.
  temp_hits=hits ; temp_power=0
  do
   if(temp_hits>=1000_i64)then
    temp_hits=temp_hits/1000_i64 ; temp_power=temp_power+3
   else
    exit
   endif
  enddo
  ! Display two digits.
  if(temp_hits<10.and.temp_power>0)then
   rhits=real(hits)*10.**(-temp_power)
   if(nint(rhits)>9)then
    hits_string='10'
   else
    hits10=int(hits/10_i64**(temp_power-1))
    hits_string=trim(i2s(hits10/10))//'.'//trim(i2s(mod(hits10,10)))
   endif
  else
   hits_string=trim(i2s(int(temp_hits)))
  endif
  select case(temp_power)
  case(0) ; suffix=''  ; case(3)  ; suffix='k' ; case(6) ; suffix='M'
  case(9) ; suffix='G' ; case(12) ; suffix='T'
  case default
   suffix='E'//trim(i2s(temp_power))
  end select
 else
  hits_string='' ; suffix=''
 endif
 write(line,"(1x,a,t30,f10.2,t42,f6.2,'%',t51,f6.2,'%',t60,f6.2,'%',t69,a3,a)")&
  &trim(label),c_time,100.*s_ratio,100.*c_ratio_branch,&
  &100.*c_ratio_tot,trim(hits_string),trim(suffix)
 END SUBROUTINE format_timer_output


 SUBROUTINE cputime(time)
 !---------------------!
 ! Interface to ETIME. !
 !---------------------!
 IMPLICIT NONE
 REAL(sp),INTENT(out) :: time(3)
 REAL(sp) user_and_system_time(2)
 time(3)=etime(user_and_system_time) ! total
 time(1)=user_and_system_time(1)     ! user
 time(2)=user_and_system_time(2)     ! system
 END SUBROUTINE cputime


 REAL(sp) FUNCTION walltime()
 !---------------------------------------!
 ! Interface to DATE_AND_TIME intrinsic. !
 !---------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,time_values(8)
 INTEGER(i64) t
 INTEGER(i64),SAVE :: time_offset
 LOGICAL,SAVE :: first_call=.true.

 ! Get current time (year, month, day, x, hour, min, sec, milliseconds)
 call date_and_time(values=time_values)

 ! Convert year to days since 1 Jan 1980.
 t=int(time_values(1)-1980,i64)*365_i64
 ! Adjust for leap years.
 do i=1980,time_values(1),4
  j=12 ; if(i==time_values(1))j=time_values(2)-1
  if(j>1.and.(mod(i,100)/=0.or.mod(i,400)==0))t=t+1_i64
 enddo ! i
 ! Convert month to days since beginning of year (potential leap day
 ! already dealt with above).
 select case(time_values(2))
 case(2)  ; t=t+31_i64  ; case(3)  ; t=t+59_i64  ; case(4)  ; t=t+90_i64
 case(5)  ; t=t+120_i64 ; case(6)  ; t=t+151_i64 ; case(7)  ; t=t+181_i64
 case(8)  ; t=t+212_i64 ; case(9)  ; t=t+243_i64 ; case(10) ; t=t+273_i64
 case(11) ; t=t+304_i64 ; case(12) ; t=t+334_i64
 end select
 ! Add day to days.
 t=t+int(time_values(3),i64)-1_i64
 ! Convert (day, hour, min, sec, milliseconds) to milliseconds since
 ! 1 Jan 1980 0:00.
 t=t*24_i64 ; t=t+time_values(5)   ! hours
 t=t*60_i64 ; t=t+time_values(6)   ! minutes
 t=t*60_i64 ; t=t+time_values(7)   ! seconds
 t=t*1000_i64 ; t=t+time_values(8) ! milliseconds

 ! Store milliseconds on first call.
 if(first_call)then
  time_offset=t
  first_call=.false.
 endif

 ! Return number of seconds - first call will always return zero.
 walltime=real(t-time_offset,sp)*1.e-3

 END FUNCTION walltime


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


 ! Error stop


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
 call traceback
 call wout()
 stop
 END SUBROUTINE errstop


 SUBROUTINE traceback
 !------------------------------------------------!
 ! Print recorded position in the code at time of !
 ! encountering an error.                         !
 !------------------------------------------------!
 IMPLICIT NONE
 TYPE(timed_routine),POINTER :: pt_temp
 call wout('Internal backtrace:')
 call wout(' Problem detected at '//trim(pt_current%name))
 pt_temp=>pt_current
 do
  if(associated(pt_temp%parent,level_zero))exit
  pt_temp=>pt_temp%parent
  call wout(' Called from '//trim(pt_temp%name))
 enddo
 nullify(pt_temp)
 call wout(' Called from MAIN')
 END SUBROUTINE traceback


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
 call traceback
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
