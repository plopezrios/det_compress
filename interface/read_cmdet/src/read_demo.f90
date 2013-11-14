PROGRAM read_demo
  !-----------------------------------------------------!
  ! Program demonstrating reading a cmdet.casl file and !
  ! evaluating the Slater matrices in the compressed    !
  ! multideterminant wave function.                     !
  !-----------------------------------------------------!
  USE casl
  IMPLICIT NONE
  INTEGER,PARAMETER :: dp=kind(1.d0)

  ! System parameters.
  INTEGER nspin,nemax,netot
  INTEGER,ALLOCATABLE :: nele(:)

  ! Compressed multi-determinant expansion data.
  CHARACTER(256) title
  INTEGER &
   &orig_norb,orig_ndet,&
   &dedup_ndet,&
   &comp_norb,comp_ndet,comp_max_nmix
  INTEGER,ALLOCATABLE :: &
     &dedup_ndetcoef(:),dedup_idetcoef(:,:),&
     &comp_nmix(:),comp_imix(:,:),comp_imix_abs(:,:),&
     &comp_nmixcoeff(:,:),comp_imixcoeff_num(:,:,:),comp_imixcoeff_den(:,:,:),&
     &comp_ndetcoef(:),comp_idetcoef_num(:,:),comp_idetcoef_den(:,:),&
     &comp_orbmap(:,:,:)
  REAL(dp),ALLOCATABLE :: &
     &orig_detcoef(:),&
     &dedup_detcoef(:),&
     &comp_mixcoeff(:,:),comp_detcoef(:)

  ! Local variables.
  INTEGER ispin,ie,ii,idet,iorb,ipos,a,csize,i,dum1,dum2
  REAL(dp) tval
  INTEGER,ALLOCATABLE :: cseed(:)
  REAL(dp),ALLOCATABLE :: rele(:,:),orig_orbval(:,:),comp_orbval(:,:),&
     &smat(:,:,:,:)

  ! Read the cmdet.casl file.
  write(6,*)'Reading cmdet.casl file.'
  call read_cmdet_casl
  write(6,*)'Done.'
  write(6,*)
  write(6,*)'Summary:'
  write(6,*)'* Title: '//trim(title)
  write(6,'(1x,a,'//trim(i2s(nspin))//'(1x,i4))')'* Particles in each &
     &spin channel:',nele(1:nspin)
  write(6,*)'* Original expansion size       : '//trim(i2s(orig_ndet))
  write(6,*)'* Original number of orbitals   : '//trim(i2s(orig_norb))
  write(6,*)'* Deduplicated expansion size   : '//trim(i2s(dedup_ndet))
  write(6,*)'* Compressed expansion size     : '//trim(i2s(comp_ndet))
  write(6,*)'* Compressed number of orbitals : '//trim(i2s(comp_norb))
  write(6,*)

  ! Post-process data.
  call prepare_cmdet_eval

  ! Generate a random electronic configuration:
  ! Allocate configuration array.
  allocate(rele(3,netot))
  ! Initialize random number generator from timer.
  call random_seed(size=csize)
  allocate(cseed(csize))
  do i=1,csize
    call system_clock(cseed(i),dum1,dum2)
  enddo
  call random_seed(put=cseed)
  deallocate(cseed)
  ! Generate random configuration in a 1x1x1 box.
  do ii=1,netot
    do a=1,3
      call random_number(rele(a,ii))
    enddo ! a
  enddo ! ii

  ! Prepare and populate array with original orbitals.
  allocate(orig_orbval(orig_norb,netot))
  orig_orbval=0.d0
  do ii=1,netot
    do iorb=1,orig_norb
      call eval_fake_orb(iorb,rele(1:3,ii),orig_orbval(iorb,ii))
    enddo ! iorb
  enddo ! ii

  ! Compute the compressed orbitals.
  allocate(comp_orbval(comp_norb,netot))
  comp_orbval=0.d0
  do ii=1,netot
    do iorb=1,comp_norb
      tval=0.d0
      do i=1,comp_nmix(iorb)
        tval=tval+orig_orbval(comp_imix_abs(i,iorb),ii)*comp_mixcoeff(i,iorb)
      enddo ! i
      comp_orbval(iorb,ii)=tval
    enddo ! iorb
  enddo ! ii

  ! Prepare array to hold the Slater matrices in the compressed expnsion.
  allocate(smat(nemax,nemax,nspin,comp_ndet))
  smat=0.d0
  ! Populate the matrices.
  do idet=1,comp_ndet
    ii=0
    do ispin=1,nspin
      do ie=1,nele(ispin)
        ii=ii+1
        do ipos=1,nele(ispin)
          iorb=comp_orbmap(ipos,ispin,idet)
          smat(ipos,ie,ispin,idet)=comp_orbval(iorb,ii)
        enddo ! ipos
      enddo ! ie
    enddo ! ispin
  enddo ! idet

  ! NB, the total wave function is:
  ! sum_idet comp_detcoef(idet) * prod_ispin
  !   determinant_of[ smat(1:nele(ispin),1:nele(ispin),ispin,idet) ]

  ! Print the first 10 matrices to show that we have them.
  write(6,*)'First '//trim(i2s(min(comp_ndet,10)))//' Slater matrices in the &
     &compressed expansion:'
  write(6,*)'(using fake plane-wave orbitals)'
  do idet=1,min(comp_ndet,10)
    write(6,*)'* Det #'//trim(i2s(idet))//':'
    do ispin=1,nspin
      write(6,*)'  * SPIN #'//trim(i2s(ispin))//':'
      do ie=1,nele(ispin)
        write(6,'(4x,'//trim(i2s(nele(ispin)))//'(1x,es12.5))')&
           &smat(1:nele(ispin),ie,ispin,idet)
      enddo ! ie
    enddo ! ispin
  enddo ! idet


CONTAINS


  SUBROUTINE read_cmdet_casl
  !--------------------------------------------------!
  ! Read the cmdet.casl file using the CASL library. !
  ! The data is loaded in global arrays defined at   !
  ! the top of this file.                            !
  ! To aid clarity, we skip most error checks.       !
  !--------------------------------------------------!
  IMPLICIT NONE
  CHARACTER(512) errmsg
  INTEGER ierr,ispin,idet,iorb,i,j,it1,dedup_max_ndetcoef,&
     &comp_max_nmixcoeff,comp_max_ndetcoef
  LOGICAL exists,is_block

  ! Instruct the CASL module to load the cmdet.casl file.  If non-blank
  ! the errmsg string will contain information on syntax errors in the
  ! CASL file on return.
  call read_casl('cmdet.casl',errmsg)
  if(len_trim(errmsg)>0)call errstop('READ_CMDET_CASL',trim(errmsg))

  ! If the cmdet.casl file was successfully read, block 'cmdet.casl:CMDET'
  ! will be present in the CASL tree.
  call query_casl_item('cmdet.casl:CMDET',exists=exists)
  if(.not.exists)call errstop('READ_CMDET_CASL','File cmdet.casl does not &
     &exist or does not contain a CMDET block.')

  ! Move the CASL context to cmdet.casl:CMDET so that we don't have to
  ! prepend "cmdet.casl:CMDET" to every request.
  call push_casl_context('cmdet.casl:CMDET')

  ! Get title if present.
  call get_casl_item('Title',title,ierr)
  if(ierr/=0)title='No title.'

  ! Get system parameters from original expansion block.
  ! The number of original determinants is the number of children of Original
  ! expansion.
  call query_casl_item('Original expansion',nchildren=orig_ndet)
  call push_casl_context('Original expansion')
  ! The number of spin channels is the number of children of Term 1
  call query_casl_item('Term 1',nchildren=it1)
  nspin=it1-1
  allocate(nele(nspin))
  nele=0
  do ispin=1,nspin
    ! The number of electrons in each spin channel is number of children of
    ! Term 1:Spin <ispin>
    call query_casl_item('Term 1:Spin '//trim(i2s(ispin)),&
     &nchildren=nele(ispin))
  enddo ! ispin
  nemax=maxval(nele)
  netot=sum(nele)
  ! Get the original expansion coefficients.
  allocate(orig_detcoef(orig_ndet))
  orig_detcoef=0.d0
  do idet=1,orig_ndet
   call get_casl_item('Term '//trim(i2s(idet)),orig_detcoef(idet),ierr)
  enddo ! idet
  ! Return to previous context.
  call pop_casl_context() ! out of "Original expansion"

  ! Load deduplication data.  This block contains blocks of the form:
  !  c_1: [ 1, 5, -7 ]
  ! meaning c_1 = k_1 + k_5 - k_7  (c = dedup. coeff., k = orig. coeff.)
  call query_casl_item('Deduplicated coefficients',nchildren=dedup_ndet)
  call push_casl_context('Deduplicated coefficients')
  ! Do a first pass to gather the maximum number of original coefficients
  ! to be combined.
  allocate(dedup_ndetcoef(dedup_ndet))
  dedup_ndetcoef=0
  do idet=1,dedup_ndet
   call query_casl_item('c_'//trim(i2s(idet)),nchildren=dedup_ndetcoef(idet))
  enddo ! idet
  ! Allocate data arrays once maximum size is known.
  dedup_max_ndetcoef=maxval(dedup_ndetcoef)
  allocate(dedup_idetcoef(dedup_max_ndetcoef,dedup_ndet))
  dedup_idetcoef=0
  ! Do a second pass to read the data.
  do idet=1,dedup_ndet
    call push_casl_context('c_'//trim(i2s(idet)))
    do i=1,dedup_ndetcoef(idet)
      call get_casl_item('%u'//trim(i2s(i)),dedup_idetcoef(i,idet),ierr)
    enddo ! i
    call pop_casl_context() ! out of "c_<idet>"
  enddo ! idet
  ! Return to previous context.
  call pop_casl_context() ! out of "Deduplicated expansion"

  ! Move into Compressed expansion block.
  call push_casl_context('Compressed expansion')

  ! Load compressed orbital pool.  This block contains blocks of the form:
  !   Orbital 1:
  !     Component 1: [ 5, Num: [ 7, -10 ], Den: [ 8, 6 ] ]
  !     Component 2: [ 9 ]
  ! meaning Phi_1 = c_7*(-c_10) / (c_8*c_6) * phi_5 + phi_9
  ! (where Phi = comp. orb., phi = orig. orb., c = dedup. coeff.)
  call query_casl_item('Orbital pool',nchildren=comp_norb)
  call push_casl_context('Orbital pool')
  ! Do a first pass to gather maxiumum array sizes.
  allocate(comp_nmix(comp_norb))
  comp_nmix=0
  comp_max_nmixcoeff=0
  do iorb=1,comp_norb
    call query_casl_item('Orbital '//trim(i2s(iorb)),nchildren=comp_nmix(iorb))
    call push_casl_context('Orbital '//trim(i2s(iorb)))
    do i=1,comp_nmix(iorb)
      call query_casl_item('Component '//trim(i2s(i))//':Num',exists=exists,&
         &is_block=is_block,nchildren=it1)
      if(exists.and.is_block)comp_max_nmixcoeff=max(comp_max_nmixcoeff,it1)
    enddo ! i
    call pop_casl_context() ! out of "Orbital <iorb>"
  enddo ! iorb
  ! Allocate data arrays once maximum size is known.
  comp_max_nmix=maxval(comp_nmix)
  allocate(comp_imix(comp_max_nmix,comp_norb),&
     &comp_nmixcoeff(comp_max_nmix,comp_norb),&
     &comp_imixcoeff_num(comp_max_nmixcoeff,comp_max_nmix,comp_norb),&
     &comp_imixcoeff_den(comp_max_nmixcoeff,comp_max_nmix,comp_norb))
  comp_imix=0
  comp_nmixcoeff=0
  comp_imixcoeff_num=0
  comp_imixcoeff_den=0
  ! Do a second pass to read the data.
  do iorb=1,comp_norb
    call push_casl_context('Orbital '//trim(i2s(iorb)))
    do i=1,comp_nmix(iorb)
      call push_casl_context('Component '//trim(i2s(i)))
      call get_casl_item('%u1',comp_imix(i,iorb),ierr)
      call query_casl_item('Num',exists=exists,is_block=is_block,nchildren=it1)
      if(exists.and.is_block)comp_nmixcoeff(i,iorb)=it1
      do j=1,comp_nmixcoeff(i,iorb)
        call get_casl_item('Num:%u'//trim(i2s(j)),&
           &comp_imixcoeff_num(j,i,iorb),ierr)
        call get_casl_item('Den:%u'//trim(i2s(j)),&
           &comp_imixcoeff_den(j,i,iorb),ierr)
      enddo ! j
      call pop_casl_context() ! out of "Component <i>"
    enddo ! i
    call pop_casl_context() ! out of "Orbital <iorb>"
  enddo ! i
  call pop_casl_context() ! out of "Orbital pool"

  ! Take a moment to infer the number of orbitals in the original orbital pool.
  orig_norb=maxval(comp_imix)

  ! Load compressed determinant data.  This block contains blocks of the form:
  !   Term 1:
  !     Coefficient:
  !       Num: [ 2, 6, 8 ]
  !       Den: [ -1, 3, 5 ]
  !     Spin 1: [ 4, 7, 8 ]
  !     Spin 2: [ 3 ]
  ! meaning:
  ! * C_1 = c_2*c_6*c_8 / ((-c_1)*c_3*c_5)
  ! * det_1_up = [ Phi_4, Phi_7, Phi_8 ]
  ! * det_1_down = [ Phi_3 ]
  ! (where C = comp. det. coeff., c = dedup. det. coeff., det_1_x = comp. det.
  ! of spin x, Phi = comp. orb.)
  call query_casl_item('Expansion',nchildren=comp_ndet)
  call push_casl_context('Expansion')
  ! Do a first pass to gather maximum array sizes.
  allocate(comp_ndetcoef(comp_ndet))
  comp_ndetcoef=0
  do idet=1,comp_ndet
    call query_casl_item('Term '//trim(i2s(idet))//':Coefficient:Num',&
       &exists=exists,is_block=is_block,nchildren=it1)
    if(exists.and.is_block)comp_ndetcoef(idet)=it1
  enddo ! idet
  ! Allocate data arrays once maximum size is known.
  comp_max_ndetcoef=maxval(comp_ndetcoef)
  allocate(comp_idetcoef_num(comp_max_ndetcoef,comp_ndet),&
     &comp_idetcoef_den(comp_max_ndetcoef,comp_ndet))
  comp_idetcoef_num=0
  comp_idetcoef_den=0
  allocate(comp_orbmap(nemax,nspin,comp_ndet))
  comp_orbmap=0
  ! Do a second pass to read the data.
  do idet=1,comp_ndet
    call push_casl_context('Term '//trim(i2s(idet)))
    call push_casl_context('Coefficient')
    do i=1,comp_ndetcoef(idet)
      call get_casl_item('Num:%u'//trim(i2s(i)),comp_idetcoef_num(i,idet),ierr)
      if(i<comp_ndetcoef(idet))then
        call get_casl_item('Den:%u'//trim(i2s(i)),comp_idetcoef_den(i,idet),&
           &ierr)
      endif
    enddo ! i
    call pop_casl_context() ! out of "Coefficient"
    do ispin=1,nspin
      call push_casl_context('Spin '//trim(i2s(ispin)))
      do iorb=1,nele(ispin)
        call get_casl_item('%u'//trim(i2s(iorb)),comp_orbmap(iorb,ispin,idet),&
           &ierr)
      enddo ! iorb
      call pop_casl_context() ! out of "Spin <ispin>"
    enddo ! ispin
    call pop_casl_context() ! out of "Term <idet>"
  enddo ! idet
  call pop_casl_context() ! out of "Expansion"

  ! Exit remaining contexts.
  call pop_casl_context() ! out of "Compressed expansion"
  call pop_casl_context() ! out of "cmdet.casl:CMDET"

  ! Free up memory used by CASL structure.
  call delete_casl_item('cmdet.casl')

  END SUBROUTINE read_cmdet_casl


  SUBROUTINE prepare_cmdet_eval
  !---------------------------------------------------------!
  ! Compute quantities required to calculate the compressed !
  ! expansion at run-time.  This only needs to be done once !
  ! after reading the compressed expansion data, and every  !
  ! time the determinant coefficients DETCOEF change.       !
  !---------------------------------------------------------!
  IMPLICIT NONE
  INTEGER idet,iorb,i,j,it1,it2
  REAL(dp) t1,t2,tcoef

  ! Compute value of deduplicated determinant coefficients.
  allocate(dedup_detcoef(dedup_ndet))
  do idet=1,dedup_ndet
    tcoef=0.d0
    do i=1,dedup_ndetcoef(idet)
      ! DEDUP_IDETCOEF absorbs sign changes, so deal with that here.
      it1=dedup_idetcoef(i,idet)
      t1=0.d0
      if(it1>0)then
        t1=orig_detcoef(it1)
      elseif(it1<0)then
        t1=-orig_detcoef(-it1)
      endif
      tcoef=tcoef+t1
    enddo ! i
    dedup_detcoef(idet)=tcoef
  enddo ! idet

  ! Compute data to transform orbitals.
  allocate(comp_imix_abs(comp_max_nmix,comp_norb),&
     &comp_mixcoeff(comp_max_nmix,comp_norb))
  do iorb=1,comp_norb
    do i=1,comp_nmix(iorb)
      comp_mixcoeff(i,iorb)=1.d0
      tcoef=1.d0
      do j=1,comp_nmixcoeff(i,iorb)
        ! COMP_IMIXCOEFF_NUM absorbs sign changes, so deal with that here.
        it1=comp_imixcoeff_num(j,i,iorb)
        t1=1.d0
        if(it1>0)then
          t1=dedup_detcoef(it1)
        elseif(it1<0)then
          t1=-dedup_detcoef(-it1)
        endif
        ! COMP_IMIXCOEFF_DEN absorbs sign changes, so deal with that here.
        it2=comp_imixcoeff_den(j,i,iorb)
        t2=1.d0
        if(it2>0)then
          t2=dedup_detcoef(it2)
        elseif(it2<0)then
          t2=-dedup_detcoef(-it2)
        endif
        if(t2/=0.d0)then
          tcoef=tcoef*(t1/t2)
        else
          tcoef=0.d0
        endif
      enddo ! j
      ! COMP_IMIX absorbs sign changes, so deal with that here.
      if(comp_imix(i,iorb)<0)then
        comp_mixcoeff(i,iorb)=-tcoef
        comp_imix_abs(i,iorb)=-comp_imix(i,iorb)
      else
        comp_mixcoeff(i,iorb)=tcoef
        comp_imix_abs(i,iorb)=comp_imix(i,iorb)
      endif
    enddo ! i
  enddo ! iorb

  ! Compute compressed expansion coefficients.
  allocate(comp_detcoef(comp_ndet))
  do idet=1,comp_ndet
    tcoef=1.d0
    do i=1,comp_ndetcoef(idet)
      ! COMP_IDETCOEF_NUM absorbs sign changes, so deal with that here.
      it1=comp_idetcoef_num(i,idet)
      t1=1.d0
      if(it1>0)then
        t1=dedup_detcoef(it1)
      elseif(it1<0)then
        t1=-dedup_detcoef(-it1)
      endif
      ! COMP_IDETCOEF_DEN absorbs sign changes, so deal with that here.
      it2=comp_idetcoef_den(i,idet)
      t2=1.d0
      if(it2>0)then
        t2=dedup_detcoef(it2)
      elseif(it2<0)then
        t2=-dedup_detcoef(-it2)
      endif
      if(t2/=0.d0)then
        tcoef=tcoef*(t1/t2)
      else
        tcoef=0.d0
      endif
    enddo ! i
    comp_detcoef(idet)=tcoef
  enddo ! idet

  END SUBROUTINE prepare_cmdet_eval


  SUBROUTINE eval_fake_orb(iorb,rvec,orbval)
  !-----------------------------------------------------!
  ! Evaluate the value ORBVAL of the IORB-th orbital at !
  ! position RVEC(1:3).                                 !
  ! Note that we don't know what the orbitals in the    !
  ! original expansion were, since cmdet.casl does not  !
  ! hold that information, so we evaluate plane waves.  !
  !-----------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: iorb
  REAL(dp),INTENT(in) :: rvec(3)
  REAL(dp),INTENT(inout) :: orbval
  REAL(dp),PARAMETER :: twopi=3.14159265358979324d0*2.d0
  INTEGER a,ik
  ! We use the orbitals:
  ! cos(2*pi*x), cos(2*pi*y), cos(2*pi*z),
  ! cos(4*pi*x), cos(4*pi*y), cos(4*pi_z), ...
  ! i.e., some weirdly excited free fluid wave function.
  ! Get plane wave direction.
  a=mod(iorb-1,3)+1
  ! Get k vector size.
  ik=(iorb-1)/3+1
  ! Compute orbital.
  orbval=cos(twopi*dble(ik)*rvec(a))
  END SUBROUTINE eval_fake_orb


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


END PROGRAM read_demo
