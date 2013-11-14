PROGRAM write_demo
  !----------------------------------------------------!
  ! Program demonstrating writing an mdet.casl file,   !
  ! which doubles as a converter from CASINO awfn.data !
  ! and correlation.data files to mdet.casl.           !
  !----------------------------------------------------!
  USE casl
  IMPLICIT NONE
  INTEGER,PARAMETER :: dp=kind(1.d0)

  ! System parameters.
  INTEGER,PARAMETER :: nspin=2
  INTEGER nemax,netot,nele(nspin)

  ! Wave function data.
  INTEGER norb,ndet
  REAL(dp),ALLOCATABLE :: detcoef(:)
  INTEGER,ALLOCATABLE :: detcoef_label(:),ecfg(:,:,:,:),orbmap(:,:,:)
  CHARACTER(80) mdet_title

  ! Read the awfn.data file.
  write(6,*)'Reading awfn.data and correlation.data.'
  call read_data

  ! Write the mdet.casl file.
  write(6,*)'Writing mdet.casl.'
  call write_mdet_casl


CONTAINS


  SUBROUTINE read_data()
  !--------------------------------------!
  ! Read awfn.data and correlation.data. !
  ! To aid clarity, we skip some error   !
  ! checks.                              !
  !--------------------------------------!
  IMPLICIT NONE
  CHARACTER(256) line
  INTEGER ierr,idum,idet,ispin,iorb,jorb,temp_ndet
  LOGICAL orb_found
  INTEGER,ALLOCATABLE :: orb_nlm(:,:)
  INTEGER,PARAMETER :: io=10

  ! Read awfn.data.
  open(io,file='awfn.data',status='old',form='formatted',iostat=ierr)
  if(ierr/=0)call errstop('READ_DATA','File awfn.data not found.')
  read(io,*)
  read(io,*)
  read(io,*)
  read(io,*)
  read(io,*)
  read(io,*)
  read(io,*)
  read(io,'(a)')line
  read(line,*,iostat=ierr)nele(1:nspin),ndet
  if(ierr/=0)call errstop('READ_DATA','The awfn.data file is not for a &
     &multi-determinant wave function.')
  if(ndet<=1)call errstop('READ_DATA','The awfn.data file is not for a &
     &multi-determinant wave function.')
  nemax=maxval(nele)
  netot=sum(nele)
  allocate(ecfg(3,nemax,nspin,ndet),detcoef(ndet),detcoef_label(ndet))
  ecfg=0 ; detcoef=0.d0 ; detcoef_label=0
  read(io,*)
  do idet=1,ndet
    do ispin=1,nspin
      do iorb=1,nele(ispin)
        read(io,*)idum,ecfg(1:3,iorb,ispin,idet)
      enddo
    enddo
  enddo
  close(io)

  ! Read correlation.data
  open(io,file='correlation.data',status='old',form='formatted',iostat=ierr)
  if(ierr/=0)call errstop('READ_DATA','File correlation.data not found.')
  do
    read(io,'(a)',iostat=ierr)line
    if(ierr/=0)call errstop('READ_DATA','The correlation.data &
       &file is not for a multideterminant wave function <1>.')
    if(trim(adjustl(line))=='START MDET')exit
  enddo
  read(io,*)
  read(io,'(a)')mdet_title
  read(io,*)
  read(io,'(a)',iostat=ierr)line
  if(ierr/=0)call errstop('READ_DATA','Problem reading correlation.data <1>.')
  if(trim(adjustl(line))/='MD')call errstop('READ_DATA','The correlation.data &
     &file is not for a multideterminant wave function <2>.')
  read(io,*,iostat=ierr)temp_ndet
  if(temp_ndet/=ndet)call errstop('READ_DATA','Mismatch between number &
     &of determinants declared in awfn.data and correlation.data.')
  do idet=1,ndet
   read(io,*)detcoef(idet),detcoef_label(idet)
  enddo
  close(io)

  ! Calculate the orbital map.
  allocate(orbmap(nemax,nspin,ndet),orb_nlm(3,ndet*netot))
  norb=0 ; orbmap=0 ; orb_nlm=0
  do idet=1,ndet
    do ispin=1,nspin
      do iorb=1,nele(ispin)
        orb_found=.false.
        do jorb=1,norb
          if(all(ecfg(1:3,iorb,ispin,idet)==orb_nlm(1:3,jorb)))then
           orb_found=.true. ; exit
          endif ! orb exists
        enddo ! jorb
        if(.not.orb_found)then ! add, assign new jorb
          norb=norb+1
          jorb=norb
          orb_nlm(1:3,norb)=ecfg(1:3,iorb,ispin,idet)
        endif ! not repeated angle
        orbmap(iorb,ispin,idet)=jorb
      enddo ! iorb
    enddo ! ispin
  enddo ! idet
  deallocate(orb_nlm)

  END SUBROUTINE read_data


  SUBROUTINE write_mdet_casl
  !-------------------------------------------------------!
  ! Write an mdet.casl file with the data we have loaded. !
  ! To aid clarity, we skip all error checks.             !
  !-------------------------------------------------------!
  IMPLICIT NONE
  CHARACTER(512) errmsg
  INTEGER idet,ispin,iorb
  call set_casl_block('mdet.casl:MDET',errmsg)
  call push_casl_context('mdet.casl:MDET')
  call set_casl_item('Title',trim(adjustl(mdet_title)),errmsg)
  call set_casl_block('Expansion',errmsg)
  call push_casl_context('Expansion')
  do idet=1,ndet
   call set_casl_block('Term '//trim(i2s(idet)),errmsg)
   call push_casl_context('Term '//trim(i2s(idet)))
   call set_casl_block('Coefficient',errmsg,prefer_inline=.true.)
   call set_casl_item('Coefficient:%u1',detcoef(idet),errmsg)
   call set_casl_item('Coefficient:Group',detcoef_label(idet),errmsg)
   do ispin=1,nspin
    if(nele(ispin)==0)cycle
    call set_casl_block('Spin '//trim(i2s(ispin)),errmsg,prefer_inline=.true.)
    call push_casl_context('Spin '//trim(i2s(ispin)))
    do iorb=1,nele(ispin)
     call set_casl_item('%u'//trim(i2s(iorb)),orbmap(iorb,ispin,idet),errmsg)
    enddo ! iorb
    call pop_casl_context() ! Spin <ispin>
   enddo ! ispin
   call pop_casl_context() ! Term <idet>
  enddo ! idet
  call pop_casl_context() ! Expansion
  call pop_casl_context() ! mdet.casl:MDET
  call write_casl('mdet.casl','mdet.casl',errmsg)
  END SUBROUTINE write_mdet_casl


 ! Tools.


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


  SUBROUTINE errstop(routine,error)
  !-----------------------------!
  ! Stop with an error message. !
  !-----------------------------!
  IMPLICIT NONE
  CHARACTER(*),INTENT(in) :: routine,error
  write(6,*)
  write(6,*)'ERROR : '//trim(routine)
  write(6,*)trim(error)
  write(6,*)
  stop
  END SUBROUTINE errstop


END PROGRAM write_demo
