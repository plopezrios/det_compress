PROGRAM det_compress
 !----------------------------------------------------------------!
 ! Compress a multideterminant expansion into a shorter expansion !
 ! by combining determinants together.                            !
 !                                                                !
 ! PLR 08.2013                                                    !
 !----------------------------------------------------------------!
 USE tools
 IMPLICIT NONE

 ! Variable declarations

 ! Data types.
 INTEGER,PARAMETER :: dp=kind(1.d0)

 ! System parameters (global variables).
 INTEGER nspin,nemax,netot
 INTEGER,ALLOCATABLE :: nele(:)

 ! Original expansion type
 TYPE original
  ! Number of original determinants and orbitals
  INTEGER ndet,norb
  ! Which orbital goes in each determinant row
  INTEGER,POINTER :: orbmap(:,:,:)=>null()
  ! Which determinant group (CSF, typically) each determinant belongs in
  INTEGER,POINTER :: detcoef_label(:)=>null()
  ! Determinant coefficients
  REAL(dp),POINTER :: detcoef(:)=>null()
 END TYPE original

 ! Deduplicated expansion type
 TYPE deduplicated
  ! Number of de-duplicated determinant.
  INTEGER ndet
  ! Which orbital goes in each de-duplicated determinant row
  INTEGER,POINTER :: orbmap(:,:,:)=>null()
  ! Which original determinant coefficients need to be added together to
  ! obtain the de-duplicated determinant coefficients.
  INTEGER,POINTER :: ndetcoef(:)=>null(),idetcoef(:,:)=>null()
  ! Effective determinant group after de-duplication
  INTEGER,POINTER :: idetlabel_eff(:)=>null()
 END TYPE deduplicated

 ! Orbital pool type (sub-type of the compressed and operation_set types).
 ! This type describes how to construct compressed orbitals from original
 ! orbitals and de-duplicated determinant coefficients.
 ! NB, we allocate in blocks to reduce the frequency of re-allocations,
 ! thus the *_alloc version of some array sizes.
 TYPE orbital_pool
  ! Number of orbitals.
  INTEGER :: norb=0
  INTEGER :: norb_alloc=0 ! allocation size
  ! Number and indices of original orbitals used in the compressed orbital.
  INTEGER,POINTER :: nmix(:)=>null(),imix(:,:)=>null()
  ! How to construct the coefficients of the linear combination from
  ! de-duplicated determinant coefficients.
  INTEGER,POINTER :: nmixcoeff(:,:)=>null(),imixcoeff_num(:,:,:)=>null(),&
   &imixcoeff_den(:,:,:)=>null()
  INTEGER :: max_nmix=0,max_nmixcoeff=0
  INTEGER :: nmix_alloc=0,nmixcoeff_alloc=0 ! allocation sizes
 END TYPE orbital_pool
 ! Block allocation deltas: we reallocate the arrays in the ORBITAL_POOL
 ! type by increasing their dimensions by the following amounts:
 INTEGER,PARAMETER :: orbpool_realloc_norb=200,orbpool_realloc_nmix=10,&
  &orbpool_realloc_nmixcoeff=5

 ! Parametrized expansion type (sub-type of the compressed and operation_set
 ! types).
 ! This type describes how to construct the compressed determinants and
 ! compressed determinant coefficients from de-duplicated determinant
 ! coefficients.  The values of %ORBMAP refer to the orbitals in a
 ! compressed orbital pool; the rest of this type is independent from
 ! the orbital pool.
 ! NB, we allocate in blocks to reduce the frequency of re-allocations,
 ! thus the *_alloc version of some array sizes.
 TYPE parametrized_expansion
  ! Number of compressed determinants and of compressed orbitals
  INTEGER :: ndet=0
  INTEGER :: ndet_alloc=0 ! allocation size
  ! Which compressed orbital goes in each compressed determinant row
  INTEGER,POINTER :: orbmap(:,:,:)=>null()
  ! How to construct compressed determinant coefficients from de-duplicated
  ! determinant coefficients.
  INTEGER,POINTER :: ndetcoef(:)=>null(),idetcoef_num(:,:)=>null(),&
   &idetcoef_den(:,:)=>null()
  INTEGER :: max_ndetcoef=1
  INTEGER :: ndetcoef_alloc=1 ! allocation size
 END TYPE parametrized_expansion
 ! Block allocation deltas: we reallocate the arrays in the
 ! PARAMETRIZED_EXPANSION type by increasing their dimensions by the
 ! following amounts:
 INTEGER,PARAMETER :: parexp_realloc_ndet=200,parexp_realloc_ndetcoef=5

 ! Expansion instances
 TYPE(original),POINTER :: orig=>null()
 TYPE(deduplicated),POINTER :: dedup=>null()
 TYPE(parametrized_expansion),POINTER :: comp=>null()
 TYPE(orbital_pool),POINTER :: orbpool=>null()

 ! List of operations.  This type describes the elements of W^(0) and
 ! U^(n) [n>0], in the notation of the paper, as well as the elements of
 ! the auxiliary set U^(0)* which contains the elements of U^(0) that are
 ! required to construct U^(1).
 ! NB, we allocate in blocks to reduce the frequency of re-allocations,
 ! thus the *_alloc version of some array sizes.
 TYPE operation_set
  ! Number of operations in set
  INTEGER :: nop=0
  INTEGER :: nop_alloc=0 ! allocation size
  ! Maximum number of determinants in an operation
  INTEGER :: max_ndet=0
  INTEGER :: ndet_alloc=0 ! allocation size
  ! Description of operation:
  ! Number of deduplicated determinants involved in each operation
  INTEGER,POINTER :: ndet(:)=>null()
  ! List of deduplicated determinants involved in each operation
  INTEGER,POINTER :: idet(:,:)=>null()
  ! Which spin determinant is combined in each operation - unused in
  ! U^(0)* and U^(n).
  INTEGER,POINTER :: ispin(:)=>null()
  ! Row index of orbital to be combined for each determinant in list
  ! - unused in U^(0)* and U^(n).
  INTEGER,POINTER :: iorb(:,:)=>null()
 END TYPE
 ! Block allocation deltas: we reallocate the arrays in the OPERATION_SET
 ! type by increasing their dimensions by the following amounts:
 INTEGER,PARAMETER :: opset_realloc_nop=10000,opset_realloc_ndet=50

 ! Title read in from mdet.casl.
 CHARACTER(256) title

 ! String to read user input.
 CHARACTER(80) char_80
 INTEGER ierr

 ! Numerical tolerances to catch various types of zero:
 ! Numbers smaller than this in absolute value are zero:
 REAL(dp),PARAMETER :: tol_zero=1.d-100
 ! Numbers differing by less than this relative amount are equal:
 REAL(dp),PARAMETER :: tol_relative=1.d-10
 ! Numbers differing by less than this absolute amount are equal:
 REAL(dp),PARAMETER :: tol_absolute=1.d-110

 ! Maximum number of quick real-valued LPSOLVE iterations.
 INTEGER,PARAMETER :: LP_SOLVE_REAL_MAX_NTRY=10

 ! Input variables:
 LOGICAL :: DO_COMPRESS=.true.
 LOGICAL :: USE_LPSOLVE=.true.
 LOGICAL :: UNIFIED_ITERATION=.true.
 LOGICAL :: IGNORE_COEFF_LABELS=.false.
 ! NB, modes in terms of DO_COMPRESS / USE_LPSOLVE / UNIFIED_ITERATION:
 ! (a) De-duplicate: F F F
 ! (b) Quick:        T F F
 ! (c) Good:         T T F
 ! (d) Best:         T T T
 ! Note that USE_LPSOLVE=F and UNIFIED_ITERATION=T are incompatible in
 ! the present implementation.

 ! Local variables
 INTEGER iter,prev_ndet

 ! Main program

 ! Print title and license info
 call wout('Multi-determinant compressor')
 call wout('============================')
 call wout()
 call wout('This program uses the lp_solve library, licensed under the LGPL.')
 call wout('Refer to the README files that accompany the source for details.')
 call wout()

 ! Load data and report
 call read_mdet_casl(orig)
 call wout('System:')
 call wout(' Title: '//trim(adjustl(title)))
 if(nspin==1)then
  call wout(' Number of electrons   : '//trim(i2s(netot))//' (same-spin)')
 elseif(nspin==2)then
  call wout(' Number of electrons   : '//trim(i2s(netot))//' ('//&
   &trim(i2s(nele(1)))//' up, '//trim(i2s(nele(2)))//' down)')
 else
  call wout(' Number of electrons   : '//trim(i2s(netot))//' (of '//&
   &trim(i2s(nspin))//' different spins)')
 endif
 call wout(' Number of determinants: '//trim(i2s(orig%ndet)))
 call wout(' Number of CSFs        : '//trim(i2s(maxval(orig%detcoef_label))))
 call wout(' Number of orbitals    : '//trim(i2s(orig%norb)))
 call wout()

 ! User interaction: pick operational level.
 call wout('Select operational level:')
 call wout('(a) de-duplicate')
 call wout('(b) quick (de-duplicate + greedy + simple iterative method')
 call wout('(c) good (de-duplicate + LPSOLVE + simple iterative method')
 call wout('(d) best (de-duplicate + LPSOLVE + unified iteration method')
 call wout()
 call wout('** TYPE A, B, C, OR D, THEN PRESS ENTER **')
 do
  read(5,*,iostat=ierr)char_80
  call wout()
  if(ierr/=0)then
   call wout('Quitting.')
   stop
  endif
  select case(trim(adjustl(char_80)))
  case('a','A')
   DO_COMPRESS=.false. ; USE_LPSOLVE=.false. ; UNIFIED_ITERATION=.false.
   call wout('Using (a) de-duplication')
  case('b','B')
   USE_LPSOLVE=.false. ; UNIFIED_ITERATION=.false.
   call wout('Using (b) quick')
  case('c','C')
   UNIFIED_ITERATION=.false.
   call wout('Using (c) good')
  case('d','D')
   call wout('Using (d) best')
  case default
   call wout('Wrong option, try again.')
   cycle
  end select
  call wout()
  exit
 enddo

 ! De-duplication
 call wout('De-duplicating expansion:')
 call deduplication(orig,dedup,comp,orbpool)
 call wout(' After de-duplication: '//trim(i2s(comp%ndet))//&
  &' determinants')
 call wout()

 ! Compression
 if(DO_COMPRESS)then

  call wout('Compressing expansion:')

  if(UNIFIED_ITERATION)then ! no need to loop

   prev_ndet=comp%ndet
   ! Compress
   call compress(orig,dedup,comp,orbpool)
   ! Report successful iteration
   if(comp%ndet<prev_ndet)call wout(' After compression: '//&
    &trim(i2s(comp%ndet))//' determinants, '//trim(i2s(orbpool%norb))//&
    &' orbitals')

  else ! .not.UNIFIED_ITERATION -> need to loop

   ! Loop over iterations
   iter=0
   do
    iter=iter+1
    prev_ndet=comp%ndet
    ! Compress
    call compress(orig,dedup,comp,orbpool)
    ! Exit if compression did not yield any (further) gains
    if(comp%ndet==prev_ndet)exit
    ! Report successful iteration
    if(iter==1)then
     call wout(' After 1 compression iteration: '//trim(i2s(comp%ndet))//&
      &' determinants, '//trim(i2s(orbpool%norb))//' orbitals')
    else ! iter>1
     call wout(' After '//trim(i2s(iter))//' compression iterations: '//&
      &trim(i2s(comp%ndet))//' determinants, '//trim(i2s(orbpool%norb))//&
      &' orbitals')
    endif ! iter==1 or not
   enddo ! iter

  endif ! UNIFIED_ITERATION or not

 else ! .not.DO_COMPRESS

  call wout('Skipping compression.')

 endif ! DO_COMPRESS or not

 call wout()

 ! Test that the compressed expansion expands to the original expansion.
 call wout('Testing compressed expansion:')
 select case(test_compression(orig,dedup,comp,orbpool))
 case(0,1)
  call wout(' Success: compressed expansion expands back to de-duplicated &
   &expansion.')
 case(2)
  call errstop('DET_COMPRESS','Compression failed: coefficient value &
   &mismatch.')
 case(3)
  call errstop('DET_COMPRESS','Compression failed: coefficient sign &
   &mismatch.')
 case(4)
  call errstop('DET_COMPRESS','Compression failed: coefficients diverge.')
 case(5)
  call errstop('DET_COMPRESS','Compression failed: reconstruction is &
   &larger than original.')
 case(6)
  call errstop('DET_COMPRESS','Compression failed: reconstruction is &
   &shorter than original.')
 case(7)
  call errstop('DET_COMPRESS','Compression failed: reconstruction and &
   &original have different terms.')
 case default
  call errstop('DET_COMPRESS','Compression failed: unknown error code')
 end select
 call wout()

 ! Test that proportionality constraints are obeyed.
 if(.not.IGNORE_COEFF_LABELS.and.maxval(orig%detcoef_label)>1)then
  call wout('Testing coefficient proportionality:')
  call test_coeff_proportionality(orig,dedup,comp,orbpool)
  call wout(' Success: compressed expansion respects coefficient &
   &proportionality.')
  call wout()
 endif

 ! Write compressed expansion
 call write_cmdet_casl(orig,dedup,comp,orbpool)


CONTAINS


 SUBROUTINE read_mdet_casl(orig)
 !-----------------------------------------------------------!
 ! Read the mdet.casl file containing the expansion details. !
 !-----------------------------------------------------------!
 USE casl
 IMPLICIT NONE
 TYPE(original),POINTER :: orig
 LOGICAL exists,is_block
 INTEGER ierr,ialloc,ispin,idet,iorb
 CHARACTER(512) errmsg

 ! Load file
 call read_casl('mdet.casl',errmsg)
 if(len_trim(errmsg)>0)call errstop('READ_MDET_CASL',trim(errmsg))
 call query_casl_item('mdet.casl:MDET',exists=exists)
 if(.not.exists)call errstop('READ_MDET_CASL','File mdet.casl not found or &
  &does not contain MDET block.')
 if(.not.push_casl_context('mdet.casl:MDET'))&
  &call errstop('READ_MDET_CASL','Could not push MDET onto the CASL stack.')

 ! Get title and basic system parameters
 call get_casl_item('Title',title,ierr=ierr)
 if(ierr/=0)title='No title.'

 ! Get basic system parameters
 allocate(orig)

 ! Get into expansion block
 call query_casl_item('Expansion',exists=exists,is_block=is_block,&
  &nchildren=orig%ndet)
 if(.not.exists.or..not.is_block)call errstop('READ_MDET_CASL',&
  &'Could not find "Expansion" block.')
 if(.not.push_casl_context('Expansion'))call errstop('READ_MDET_CASL',&
  &'Could not push Expansion onto the CASL stack.')

 ! Get system size from first term
 nspin=0
 do
  call query_casl_item('Term 1:Spin '//trim(i2s(nspin+1)),exists=exists,&
   &is_block=is_block)
  if(.not.exists.or..not.is_block)exit
  nspin=nspin+1
 enddo
 if(nspin==0)call errstop('READ_MDET_CASL','Could not find number of spin &
  &channels in the system, possibly "Term 1" block is incomplete.')
 allocate(nele(nspin),stat=ialloc)
 call check_alloc(ialloc,'READ_MDET_CASL','nele')
 nele=0
 do ispin=1,nspin
  call query_casl_item('Term 1:Spin '//trim(i2s(ispin)),nchildren=nele(ispin))
 enddo ! ispin
 nemax=maxval(nele)
 netot=sum(nele)

 ! Allocate orbmap
 allocate(orig%orbmap(nemax,nspin,orig%ndet),orig%detcoef(orig%ndet),&
  &orig%detcoef_label(orig%ndet),stat=ialloc)
 call check_alloc(ialloc,'READ_MDET_CASL','orbmap, detcoef')
 orig%orbmap=0
 orig%detcoef=0.d0
 orig%detcoef_label=0

 ! Read expansion coefficients and orbmap
 do idet=1,orig%ndet
  if(.not.push_casl_context('Term '//trim(i2s(idet))))call errstop&
   &('READ_MDET_CASL','Could not push Term <idet> onto the CASL stack.')
  call get_casl_item('Coefficient:%u1',orig%detcoef(idet),ierr=ierr)
  if(ierr/=0)call errstop('READ_MDET_CASL','Could not read detcoef(<idet>).')
  call get_casl_item('Coefficient:Group',orig%detcoef_label(idet),ierr=ierr)
  if(ierr/=0)call errstop('READ_MDET_CASL','Could not read &
   &detcoef_label(<idet>).')
  do ispin=1,nspin
   if(nele(ispin)==0)cycle
   if(.not.push_casl_context('Spin '//trim(i2s(ispin))))call errstop&
    &('READ_MDET_CASL','Could not push Spin <ispin> onto the CASL stack.')
   do iorb=1,nele(ispin)
    call get_casl_item('%u'//trim(i2s(iorb)),orig%orbmap(iorb,ispin,idet),&
     &ierr=ierr)
    if(ierr/=0)call errstop('READ_MDET_CASL','Could not read orbmap(...).')
   enddo ! iorb
   if(.not.pop_casl_context())call errstop('READ_MDET_CASL',&
    &'Could not pop out of Spin <ispin> in the CASL stack.')
  enddo ! ispin
  if(.not.pop_casl_context())call errstop('READ_MDET_CASL',&
   &'Could not pop out of Term <idet> in the CASL stack.')
 enddo ! idet
 if(.not.pop_casl_context())call errstop('READ_MDET_CASL',&
  &'Could not pop out of Expansion in the CASL stack.')
 if(.not.pop_casl_context())call errstop('READ_MDET_CASL',&
  &'Could not pop out of MDET in the CASL stack.')
 if(pop_casl_context())call errstop('READ_MDET_CASL',&
  &'Got lost in CASL context changes.')

 ! Get number of orbitals
 orig%norb=maxval(orig%orbmap)

 ! Reclaim CASL memory and clean up.
 call delete_casl_item('mdet.casl')

 END SUBROUTINE read_mdet_casl


 SUBROUTINE deduplication(orig,dedup,comp,orbpool)
 !--------------------------------------------------------------------!
 ! Remove redundant terms from the initial expansion, combining their !
 ! coefficients.  While at it, construct the data structures to hold  !
 ! the compressed expansion.                                          !
 !--------------------------------------------------------------------!
 IMPLICIT NONE
 TYPE(original),POINTER :: orig
 TYPE(deduplicated),POINTER :: dedup
 TYPE(parametrized_expansion),POINTER :: comp
 TYPE(orbital_pool),POINTER :: orbpool
 ! Local variables
 INTEGER idet,jdet,ispin,iorb,jorb,i,j,k,iswp,jeff,jeff0,nm,&
  &max_ndetcoef,ndetcoef_alloc,ndetcoef_alloc_delta
 LOGICAL det_active(orig%ndet)
 REAL(dp) t1,t2
 ! Work pointers.
 INTEGER,POINTER :: tmp_idet_eff(:)=>null(),dedup_idetlabel(:,:)=>null()
 ! Sections of DEDUP_IDETLABEL for faster access.
 INTEGER,POINTER :: dedup_idetlabel_ndet(:)=>null(),&
  &dedup_idetlabel_idet(:)=>null(),dedup_idetlabel_jdet(:)=>null()
 ! Sections of components of DEDUP for faster access.
 INTEGER,POINTER :: dedup_orbmap_idet(:,:)=>null(),&
  &dedup_idetcoef_1_idet=>null(),dedup_idetcoef_ndet(:)=>null(),&
  &dedup_idetcoef_idet(:)=>null(),dedup_idetcoef_jdet(:)=>null()

 ! The following makes all of our examples do a single allocation without
 ! massively overestimating the size of the array.
 ndetcoef_alloc_delta=nint(sqrt(sqrt(dble(orig%ndet))))
 if(ndetcoef_alloc_delta<5)ndetcoef_alloc_delta=5
 if(ndetcoef_alloc_delta>20)ndetcoef_alloc_delta=20

 ! Initialize variables and pointers.
 det_active=.true.
 ndetcoef_alloc=ndetcoef_alloc_delta
 max_ndetcoef=0
 allocate(dedup)
 dedup%ndet=0
 call resize_pointer((/orig%ndet/),dedup%ndetcoef)
 call resize_pointer((/ndetcoef_alloc,orig%ndet/),dedup%idetcoef)
 call resize_pointer((/nemax,nspin,orig%ndet/),dedup%orbmap)
 if(.not.IGNORE_COEFF_LABELS)call resize_pointer(&
  &(/ndetcoef_alloc,orig%ndet/),dedup_idetlabel)

 ! Copy and sort the original orbmap.  We use dedup%orbmap and
 ! dedup%idetcoef(1,:) as temporary storage space.  Subsequent access
 ! to these arrays is safe: values are read for term idet and written
 ! for term jdet where jdet <= idet, avoiding mix-ups of temporary and
 ! final data.  The final arrays are reallocated to the correct size
 ! and any remaining temporary data are discarded.
 do idet=1,orig%ndet
  dedup_orbmap_idet=>dedup%orbmap(:,:,idet)
  dedup_orbmap_idet=orig%orbmap(:,:,idet)
  dedup_idetcoef_1_idet=>dedup%idetcoef(1,idet)
  dedup_idetcoef_1_idet=idet
  do ispin=1,nspin
   do iorb=1,nele(ispin)
    iswp=iorb
    do jorb=iorb+1,nele(ispin)
     if(dedup_orbmap_idet(iswp,ispin)>dedup_orbmap_idet(jorb,ispin))iswp=jorb
    enddo ! jorb
    if(iswp>iorb)then
     call iswap1(dedup_orbmap_idet(iswp,ispin),dedup_orbmap_idet(iorb,ispin))
     dedup_idetcoef_1_idet=-dedup_idetcoef_1_idet
    endif
   enddo ! iorb
  enddo ! ispin
 enddo ! idet
 ! Clean shortcut pointers.
 nullify(dedup_orbmap_idet,dedup_idetcoef_1_idet)

 ! Merge terms with the same determinants and different coefficients.
 do idet=1,orig%ndet
  if(.not.det_active(idet))cycle
  dedup%ndet=dedup%ndet+1
  dedup%ndetcoef(dedup%ndet)=1
  max_ndetcoef=max(1,max_ndetcoef)
  dedup_orbmap_idet=>dedup%orbmap(:,:,idet)
  dedup%idetcoef(1,dedup%ndet)=dedup%idetcoef(1,idet) ! see sorting code above
  dedup%orbmap(:,:,dedup%ndet)=dedup_orbmap_idet(:,:) ! see sorting code above
  if(.not.IGNORE_COEFF_LABELS)dedup_idetlabel(1,dedup%ndet)=&
   &orig%detcoef_label(idet)
  do jdet=idet+1,orig%ndet
   if(.not.det_active(jdet))cycle
   if(all(dedup_orbmap_idet(:,:)==dedup%orbmap(:,:,jdet)))then
    dedup%ndetcoef(dedup%ndet)=dedup%ndetcoef(dedup%ndet)+1
    max_ndetcoef=max(dedup%ndetcoef(dedup%ndet),max_ndetcoef)
    if(max_ndetcoef>ndetcoef_alloc)then
     ndetcoef_alloc=ndetcoef_alloc+ndetcoef_alloc_delta
     call resize_pointer((/ndetcoef_alloc,orig%ndet/),dedup%idetcoef)
     if(.not.IGNORE_COEFF_LABELS)&
      &call resize_pointer((/ndetcoef_alloc,orig%ndet/),dedup_idetlabel)
    endif
    dedup%idetcoef(dedup%ndetcoef(dedup%ndet),dedup%ndet)=&
     &dedup%idetcoef(1,jdet) ! see sorting code above
    if(.not.IGNORE_COEFF_LABELS)&
     &dedup_idetlabel(dedup%ndetcoef(dedup%ndet),dedup%ndet)=&
     &orig%detcoef_label(jdet)
    det_active(jdet)=.false.
   endif
  enddo ! jdet
  ! Sort list of coefficients by group, then by index
  if(.not.IGNORE_COEFF_LABELS)then
   dedup_idetlabel_ndet=>dedup_idetlabel(:,dedup%ndet)
   dedup_idetcoef_ndet=>dedup%idetcoef(:,dedup%ndet)
   do i=1,dedup%ndetcoef(dedup%ndet)
    iswp=i
    do j=i+1,dedup%ndetcoef(dedup%ndet)
     if(dedup_idetlabel_ndet(iswp)>dedup_idetlabel_ndet(j))then
      iswp=j
     elseif(dedup_idetlabel_ndet(iswp)==dedup_idetlabel_ndet(j))then
      if(abs(dedup_idetcoef_ndet(iswp))>abs(dedup_idetcoef_ndet(j)))iswp=j
     endif
    enddo ! j
    if(iswp>i)then
     call iswap1(dedup_idetlabel_ndet(iswp),dedup_idetlabel_ndet(i))
     call iswap1(dedup_idetcoef_ndet(iswp),dedup_idetcoef_ndet(i))
    endif
   enddo ! i
  endif ! .not.IGNORE_COEFF_LABELS
 enddo ! idet
 ! Clean shortcut pointers.
 nullify(dedup_orbmap_idet,dedup_idetcoef_ndet)

 ! Trim pointers to final size
 call resize_pointer((/dedup%ndet/),dedup%ndetcoef)
 call resize_pointer((/max_ndetcoef,dedup%ndet/),dedup%idetcoef)
 call resize_pointer((/nemax,nspin,dedup%ndet/),dedup%orbmap)

 ! Create list of "effective" coefficient labels
 call resize_pointer((/dedup%ndet/),dedup%idetlabel_eff)

 if(.not.IGNORE_COEFF_LABELS)then

  call resize_pointer((/dedup%ndet/),tmp_idet_eff)
  jeff0=maxval(orig%detcoef_label)
  jeff=jeff0
  do idet=1,dedup%ndet
   dedup_idetcoef_idet=>dedup%idetcoef(:,idet)
   dedup_idetlabel_idet=>dedup_idetlabel(:,idet)
   dedup%idetlabel_eff(idet)=dedup_idetlabel_idet(1)
   nm=dedup%ndetcoef(idet)
   if(nm==1)cycle
   if(all(dedup_idetlabel_idet(2:nm)==dedup_idetlabel_idet(1)))cycle
   ! Deduplicated coefficient needs an "effective" label
   do j=jeff0+1,jeff
    jdet=tmp_idet_eff(j-jeff0)
    dedup_idetcoef_jdet=>dedup%idetcoef(:,jdet)
    dedup_idetlabel_jdet=>dedup_idetlabel(:,jdet)
    nm=dedup%ndetcoef(jdet)
    if(dedup%ndetcoef(idet)/=nm)cycle
    if(any(dedup_idetlabel_idet(1:nm)/=dedup_idetlabel_jdet(1:nm)))cycle
    ! Check proportionality
    t1=orig%detcoef(abs(dedup_idetcoef_idet(1)))/&
     &orig%detcoef(abs(dedup_idetcoef_jdet(1)))
    if(dedup_idetcoef_idet(1)<0.neqv.dedup_idetcoef_jdet(1)<0)t1=-t1
    do k=2,nm
     t2=orig%detcoef(abs(dedup_idetcoef_idet(k)))/&
      &orig%detcoef(abs(dedup_idetcoef_jdet(k)))
     if(dedup_idetcoef_idet(k)<0.neqv.dedup_idetcoef_jdet(k)<0)t2=-t2
     if(.not.compare_numbers(t1,t2))exit
    enddo ! k
    if(k<=nm)cycle
    exit
   enddo ! jdet
   dedup%idetlabel_eff(idet)=j
   if(j>jeff)then ! create new "effective" label
    jeff=jeff+1
    tmp_idet_eff(jeff-jeff0)=idet
   endif
  enddo ! idet

  ! Clean up.
  deallocate(tmp_idet_eff,dedup_idetlabel)
  nullify(tmp_idet_eff,dedup_idetlabel)
  nullify(dedup_idetcoef_idet,dedup_idetcoef_jdet,dedup_idetlabel_idet,&
   &dedup_idetlabel_jdet)

 else ! IGNORE_COEFF_LABELS

  ! Make all effective labels equal.
  dedup%idetlabel_eff(1:dedup%ndet)=1

 endif ! .not.IGNORE_COEFF_LABELS or IGNORE_COEFF_LABELS

 ! Construct initial compressed expansion (=deduplicated).
 allocate(comp)
 comp%ndet=dedup%ndet   ;  comp%ndet_alloc=dedup%ndet
 call resize_pointer((/dedup%ndet/),comp%ndetcoef,init=1)
 call resize_pointer((/1,dedup%ndet/),comp%idetcoef_num)
 call resize_pointer((/1,dedup%ndet/),comp%idetcoef_den)
 call resize_pointer((/nemax,nspin,dedup%ndet/),comp%orbmap)
 comp%orbmap(:,:,:)=dedup%orbmap(:,:,:)
 do idet=1,dedup%ndet
  comp%idetcoef_num(1,idet)=idet
 enddo ! idet

 ! Construct initial orbpool (=original orbitals).
 allocate(orbpool)
 orbpool%norb=orig%norb ; orbpool%norb_alloc=orig%norb
 orbpool%max_nmix=1     ; orbpool%nmix_alloc=1
 call resize_pointer((/orig%norb/),orbpool%nmix,init=1)
 call resize_pointer((/1,orig%norb/),orbpool%imix)
 call resize_pointer((/1,orig%norb/),orbpool%nmixcoeff)
 call resize_pointer((/0,1,orig%norb/),orbpool%imixcoeff_num)
 call resize_pointer((/0,1,orig%norb/),orbpool%imixcoeff_den)
 do iorb=1,orig%norb
  orbpool%imix(1,iorb)=iorb
 enddo ! iorb

 END SUBROUTINE deduplication


 SUBROUTINE compress(orig,dedup,comp,orbpool)
 !-------------------------------------------------------!
 ! Compress the de-duplicated multideterminant expansion !
 ! by combining determinants together.                   !
 !-------------------------------------------------------!
 IMPLICIT NONE
 TYPE(original),POINTER :: orig
 TYPE(deduplicated),POINTER :: dedup
 TYPE(parametrized_expansion),POINTER :: comp
 TYPE(orbital_pool),POINTER :: orbpool
 ! PAREXP and OPSET instances.
 TYPE(parametrized_expansion),POINTER :: new=>null(),comp_u0star=>null(),&
  &comp_un=>null()
 TYPE(operation_set),POINTER :: opset_pairwise=>null(),opset_w0=>null(),&
  &opset_u0star=>null(),opset_un=>null()
 ! Variables for pairwise comparison and partitioning.
 INTEGER maxtag,tag_idet,itag,tag_jdet,ntag_replace,ireplace
 INTEGER dettag(comp%ndet),tag2part(comp%ndet),detpartition(comp%ndet),&
  &ndet_in_partition(comp%ndet),itag_replace(comp%ndet)
 ! Variables for greedy algorithm.
 INTEGER greedy_indx(netot*comp%ndet)
 ! Variables for constructing and solving linear program.
 INTEGER itry,ivar,ieqn
 LOGICAL anyvarfail
 INTEGER det2eqn_map(comp%ndet),det2eqn_map_rec(comp%ndet)
 LOGICAL,POINTER :: varfail(:)=>null()
 ! Flag to determine whether determinants have been compressed or not.
 LOGICAL detactive(comp%ndet)
 ! Variables for LP_SOLVE
 INTEGER,POINTER :: lp_vartype(:)=>null(),lp_ineq_sign(:)=>null()
 REAL(dp),POINTER :: lp_array(:,:)=>null(),lp_rhs(:)=>null(),&
  &lp_targetcoeff(:)=>null(),lp_solution(:)=>null()
 INTEGER lp_nvar,lp_neqn,lp_target_opt_sign,lp_solution_status
 ! Miscellaneous local variables.
 INTEGER idet,jdet,iorb,jorb,ispin,iop,jop,nop0,nop0_prev,ilevel,npart,ipart,&
  &ndetcoef1
 INTEGER,POINTER :: op_idet(:)=>null()
 ! Sections of components of COMP for faster access.
 INTEGER,POINTER :: comp_orbmap_idet(:,:)=>null(),comp_orbmap_jdet(:,:)=>null()

 ! Interface to wrapper around LP_SOLVE
 INTERFACE
  SUBROUTINE lp_wrapper(neqn,nvar,vartype,lparray,ineq_sign,lprhs,targetcoeff,&
   &target_opt_sign,solution,solution_status)
  INTEGER,INTENT(in) :: neqn,nvar,vartype(nvar),ineq_sign(neqn),target_opt_sign
  REAL(kind(1.d0)),INTENT(in) :: lparray(nvar,neqn),lprhs(neqn),&
   &targetcoeff(nvar)
  REAL(kind(1.d0)),INTENT(out) :: solution(neqn)
  INTEGER,INTENT(out) :: solution_status
  END SUBROUTINE lp_wrapper
 END INTERFACE

 ! Form all possible sets of compressible terms
 dettag=0
 maxtag=0
 do idet=1,comp%ndet
  comp_orbmap_idet=>comp%orbmap(:,:,idet)
  if(dettag(idet)>0)then
   tag_idet=dettag(idet)
  else
   tag_idet=maxtag+1
  endif
  ntag_replace=1
  itag_replace(1)=tag_idet
  ! Build list of determinants which IDET can be combined with
  do jdet=idet+1,comp%ndet
   comp_orbmap_jdet=>comp%orbmap(:,:,jdet)
   ! Find whether IDET and JDET can be merged.
   call compare_pair(comp_orbmap_idet,comp_orbmap_jdet,iorb,jorb,ispin)
   if(ispin/=0)then
    ! They can, so add an operation to the pairwise operation set, add
    ! JDET to the list, and replace IDET's tentative tag with JDET's tag if
    ! it is a lower number.
    call add_op(opset_pairwise,2,(/idet,jdet/),(/iorb,jorb/),ispin)
    tag_jdet=dettag(jdet)
    dettag(jdet)=tag_idet
    if(tag_jdet>0)then
     ! Replace IDET's tentative tag with JDET's tag if the latter is
     ! a lower number.
     tag_idet=min(tag_idet,tag_jdet)
     ! See if JDET's tag is in the list of tags; if not, add it.
     do ireplace=1,ntag_replace
      if(itag_replace(ireplace)==tag_jdet)exit
     enddo ! ireplace
     if(ireplace>ntag_replace)then
      ntag_replace=ntag_replace+1
      itag_replace(ntag_replace)=tag_jdet
     endif ! new tag to replace
    endif
   endif
  enddo ! jdet
  ! Set IDET's tag.
  dettag(idet)=tag_idet
  ! Propagate connections by replacing all listed tags with the final tag we've
  ! found for IDET.
  if(ntag_replace>1)then
   do jdet=1,comp%ndet
    do ireplace=1,ntag_replace
     if(itag_replace(ireplace)==dettag(jdet))then
      dettag(jdet)=tag_idet
      exit
     endif
    enddo ! ireplace
   enddo ! jdet
  endif ! any tags to replace
  ! Update maximum tag value
  maxtag=max(maxtag,tag_idet)
 enddo ! idet

 ! Clean up.
 nullify(comp_orbmap_idet,comp_orbmap_jdet)

 ! Return if no pairwise operations found.
 if(.not.associated(opset_pairwise))return

 ! Construct map from tags to partitions
 tag2part=0
 npart=0
 do idet=1,comp%ndet
  itag=dettag(idet)
  if(itag>0)then
   if(tag2part(itag)==0)then
    npart=npart+1
    tag2part(itag)=npart
   endif
  endif
 enddo ! idet
 ! Assign definitive tags
 detpartition=0
 ndet_in_partition(1:npart)=0
 do idet=1,comp%ndet
  itag=dettag(idet)
  if(itag>0)then
   ipart=tag2part(itag)
   detpartition(idet)=ipart
   if(ipart>0)ndet_in_partition(ipart)=ndet_in_partition(ipart)+1
  endif
 enddo ! idet

 ! Report maximum partition size.
 call wout(' '//trim(i2s(npart))//' partitions of maximum size '//&
  &trim(i2s(maxval(ndet_in_partition(1:npart)))))

 ! Flag all determinants in COMP as active (we use this at the end to
 ! re-add uncompressible determinants to the final expansion).
 detactive=.true.

 ! Initialize NEW expansion.
 allocate(new)

 ! Loop over partitions
 do ipart=1,npart

  if(ndet_in_partition(ipart)<2)cycle

  ! Generate compression operations
  do iop=1,opset_pairwise%nop
   if(opset_pairwise%ndet(iop)<2)cycle
   idet=opset_pairwise%idet(1,iop)
   if(detpartition(idet)/=ipart)cycle
   jdet=opset_pairwise%idet(2,iop)
   iorb=opset_pairwise%iorb(1,iop)
   jorb=opset_pairwise%iorb(2,iop)
   ispin=opset_pairwise%ispin(iop)
   call merge_op(opset_w0,idet,jdet,iorb,jorb,ispin)
  enddo ! iop
  if(.not.associated(opset_w0))cycle
  if(opset_w0%nop<1)cycle

  ! Proceed through recursion levels for this partition.
  if(UNIFIED_ITERATION)then
   ! Prepare for recursion level 0 -> 1
   call construct_u0star(orig,dedup,comp,opset_w0,opset_u0star,comp_u0star,&
    &orbpool)
   if(associated(opset_u0star))then
    if(opset_u0star%nop>0)then
     ! Recursion level 0 -> 1
     call construct_un(1,opset_u0star%nop,dedup,opset_u0star,comp_u0star,&
      &opset_un,comp_un,orbpool)
     if(associated(opset_un))then
      if(opset_un%nop>0)then
       nop0=0
        do ilevel=2,netot
         ! Recursion level ILEVEL-1 -> ILEVEL
         nop0_prev=nop0
         nop0=opset_un%nop
         call construct_un(nop0_prev+1,nop0,dedup,opset_un,comp_un,opset_un,&
          &comp_un,orbpool)
         if(opset_un%nop==nop0)exit
        enddo ! ilevel
      endif ! recursion 1->2 possible (2)
     endif ! recursion 1->2 possible (1)
    endif ! recursion 0->1 possible (2)
   endif ! recursion 0->1 possible (1)
  endif ! UNIFIED_ITERATION

  ! Solve linear program
  if(USE_LPSOLVE)then

   ! Generate system of equations
   ieqn=0
   do idet=1,comp%ndet
    if(detpartition(idet)==ipart)then
     ieqn=ieqn+1 ; det2eqn_map(idet)=ieqn
    endif
   enddo ! idet
   if(UNIFIED_ITERATION)then
    do idet=1,comp%ndet
     if(detpartition(idet)==ipart)then
      ieqn=ieqn+1 ; det2eqn_map_rec(idet)=ieqn
     endif
    enddo ! idet
   endif
   ! Problem size
   lp_nvar=opset_w0%nop
   if(associated(opset_un))lp_nvar=lp_nvar+opset_un%nop
   lp_neqn=ieqn
   ! Allocate LP problem arrays
   call resize_pointer((/lp_nvar,lp_neqn/),lp_array)
   call resize_pointer((/lp_neqn/),lp_rhs)
   call resize_pointer((/lp_nvar/),lp_vartype)
   call resize_pointer((/lp_neqn/),lp_ineq_sign)
   call resize_pointer((/lp_nvar/),lp_targetcoeff)
   call resize_pointer((/lp_nvar/),lp_solution)
   ! Populate equations
   do iop=1,opset_w0%nop
    lp_targetcoeff(iop)=1.d0
    do jdet=1,opset_w0%ndet(iop)
     idet=opset_w0%idet(jdet,iop)
     ieqn=det2eqn_map(idet)
     lp_array(iop,ieqn)=1.d0
     lp_rhs(ieqn)=1.d0
     lp_ineq_sign(ieqn)=1
    enddo
   enddo
   if(associated(opset_un))then
    do iop=1,opset_un%nop
     lp_targetcoeff(iop+opset_w0%nop)=1.d0
     do jdet=1,opset_un%ndet(iop)
      ! 1st inequation
      idet=opset_un%idet(jdet,iop)
      ieqn=det2eqn_map(idet)
      lp_array(iop+opset_w0%nop,ieqn)=1.d0
      lp_rhs(ieqn)=1.d0
      lp_ineq_sign(ieqn)=1
      ! 2nd inequation
      ieqn=det2eqn_map_rec(idet)
      lp_array(iop+opset_w0%nop,ieqn)=1.d0
      lp_rhs(ieqn)=1.d0
      lp_ineq_sign(ieqn)=-1
     enddo
    enddo
   endif
   lp_target_opt_sign=-1 ! minimize

   ! Solve the linear program by calling LPSOLVE.
   call resize_pointer((/lp_nvar/),varfail)
   itry=1
   do
    call lp_wrapper(lp_neqn,lp_nvar,lp_vartype,lp_array,lp_ineq_sign,lp_rhs,&
     &lp_targetcoeff,lp_target_opt_sign,lp_solution,lp_solution_status)
    if(lp_solution_status/=0)call errstop('COMPRESS_EXPANSION',&
     &'LP_SOLVE complains and don''t know what to do..')
    ! Check the real-valued problem gives a binary-valued solution.
    anyvarfail=.false.
    do ivar=1,lp_nvar
     if(abs(lp_solution(ivar))>1.d-7.and.abs(lp_solution(ivar)-1.d0)>1.d-7)then
      varfail(ivar)=.true.
      anyvarfail=.true.
     else
      varfail(ivar)=.false.
     endif
    enddo ! ivar
    ! If it does, this is the solution.
    if(.not.anyvarfail)exit
    ! If not, flag the problematic variables as binary and try again.
    itry=itry+1
    if(itry<=LP_SOLVE_REAL_MAX_NTRY)then
     where(varfail)lp_vartype=2
    else
     ! Too many attempts, flag all as binary and solve the exact problem.
     lp_vartype(:)=2
    endif
    lp_solution=0.d0
    lp_solution_status=0
   enddo ! itry
   deallocate(varfail)
   nullify(varfail)

   ! Apply operations contained in solution
   if(lp_solution_status==0)then
    if(UNIFIED_ITERATION.and.associated(opset_un))then
     ! Apply operations in U^(n) first.
     do iop=1,opset_un%nop
      if(nint(lp_solution(iop+opset_w0%nop))==1)then
       ! Apply operation by simply copying the IDET-th orbmap from the
       ! subclass of OPSET.
       call apply_op_un(iop,comp_un,new)
       ! Remove determinants involved in operation from elements in W^(0).
       op_idet=>opset_un%idet(1:opset_un%ndet(iop),iop)
       call update_opset(op_idet,opset_w0,detactive)
      endif
     enddo
    endif
    ! Apply operations in W^(0)
    do iop=1,opset_w0%nop
     if(nint(lp_solution(iop))==1.and.opset_w0%ndet(iop)>1)then
      ! Apply the operation.
      call apply_op(opset_w0%ndet(iop),&
       &opset_w0%idet(1:opset_w0%ndet(iop),iop),&
       &opset_w0%ispin(iop),opset_w0%iorb(1:opset_w0%ndet(iop),iop),orig,&
       &dedup,comp,new,orbpool)
      ! Remove determinants involved in operation from elements in W^(0).
      op_idet=>opset_w0%idet(1:opset_w0%ndet(iop),iop)
      call update_opset(op_idet,opset_w0,detactive)
     endif
    enddo ! iop
   endif

   ! Clean up pointers used for the LPSOLVE call for this partition
   deallocate(lp_vartype,lp_array,lp_rhs,lp_ineq_sign,lp_targetcoeff,&
    &lp_solution)
   nullify(lp_vartype,lp_array,lp_rhs,lp_ineq_sign,lp_targetcoeff,&
    &lp_solution)

  else ! .not.USE_LPSOLVE

   ! Apply the greedy algorithm.
   ! NB, by construction this will not work with the unified iteration
   ! algorithm.  Besides using opset_un below, one would need the ability
   ! to invalidate operations in U^(n), which we don't have at present.

   ! Initialize pointers, and set operation index to 1..n
   do iop=1,opset_w0%nop
    greedy_indx(iop)=iop
   enddo ! i

   ! Loop over operations
   do jop=opset_w0%nop,1,-1
    ! Sort remaining JOP operations by number of determinants.  This operation
    ! is costs ~ O[JOP log(JOP)], but for JOP<OPSET_NOP the list should be
    ! almost sorted and should be more like ~ O[JOP].
    call qsort3_int_partial_preinit(opset_w0%nop,1,jop,opset_w0%ndet,&
     &greedy_indx)
    ! Pick last (largest) operation in remaining set
    iop=greedy_indx(jop)
    ! Quit if it's a no-op - we are done
    if(opset_w0%ndet(iop)<2)exit
    ! Apply the operation.
    call apply_op(opset_w0%ndet(iop),&
     &opset_w0%idet(1:opset_w0%ndet(iop),iop),&
     &opset_w0%ispin(iop),opset_w0%iorb(1:opset_w0%ndet(iop),iop),orig,&
     &dedup,comp,new,orbpool)
    ! Remove determinants involved in operation from elements in W^(0).
    op_idet=>opset_w0%idet(1:opset_w0%ndet(iop),iop)
    call update_opset(op_idet,opset_w0,detactive)
   enddo ! jop

  endif ! USE_LPSOLVE or not

  ! Clean up.
  nullify(op_idet)
  ! Avoid subsequent reallocations of OPSET_*, instead zero counters.
  opset_w0%nop=0
  opset_w0%max_ndet=0
  if(associated(opset_u0star))then
   opset_u0star%nop=0
   opset_u0star%max_ndet=0
   comp_u0star%ndet=0
   comp_u0star%max_ndetcoef=1
  endif
  if(associated(opset_un))then
   opset_un%nop=0
   opset_un%max_ndet=0
   comp_un%ndet=0
   comp_un%max_ndetcoef=1
  endif

 enddo ! ipart

 ! Free memory from operation sets, etc.
 call deallocate_opset(opset_w0)
 call deallocate_opset(opset_u0star)
 call deallocate_parexp(comp_u0star)
 call deallocate_opset(opset_un)
 call deallocate_parexp(comp_un)

 ! Add any remaining active determinants to expansion.
 do idet=1,comp%ndet
  if(.not.detactive(idet))cycle
  new%ndet=new%ndet+1
  call parexp_make_room(new)
  ! Set detcoef to that of first determinant in the operation.
  ndetcoef1=comp%ndetcoef(idet)
  new%ndetcoef(new%ndet)=ndetcoef1
  new%idetcoef_num(1:ndetcoef1,new%ndet)=&
   &comp%idetcoef_num(1:ndetcoef1,idet)
  new%idetcoef_den(1:ndetcoef1,new%ndet)=&
   &comp%idetcoef_den(1:ndetcoef1,idet)
  ! Copy the orbmap.
  new%orbmap(1:nemax,1:nspin,new%ndet)=comp%orbmap(1:nemax,1:nspin,idet)
 enddo ! idet

 ! Purge unused orbitals from orbpool.
 call purge_orbpool(orbpool,new)

 ! Clean up pointers.
 call deallocate_opset(opset_pairwise)

 ! Point comp at new.
 deallocate(comp%orbmap,comp%ndetcoef,comp%idetcoef_num,comp%idetcoef_den)
 deallocate(comp)
 comp=>new
 nullify(new)

 END SUBROUTINE compress


 SUBROUTINE compare_pair(comp_orbmap_idet,comp_orbmap_jdet,diff_iorb,&
  &diff_jorb,diff_ispin)
 !---------------------------------------------------------!
 ! Compare terms IDET and JDET, and if they differ only by !
 ! one element, return its location on each determinant in !
 ! DIFF_IORB and DIFF_JORB, and in which spin determinant  !
 ! this mismatch occurs in DIFF_ISPIN, else return         !
 ! DIFF_ISPIN=0.                                           !
 !---------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,POINTER :: comp_orbmap_idet(:,:),comp_orbmap_jdet(:,:)
 INTEGER,INTENT(out) :: diff_iorb,diff_jorb,diff_ispin
 INTEGER i,iorb,ispin
 ! Sections of COMP_ORBMAP_* for faster access.
 INTEGER,POINTER :: comp_orbmap_idet_ispin(:)=>null(),&
  &comp_orbmap_jdet_ispin(:)=>null()


 ! Initialize to "equal determinants".
 diff_iorb=0 ; diff_jorb=0 ; diff_ispin=-1

 ! Loop over the NSPIN determinants.
 do ispin=1,nspin
  comp_orbmap_idet_ispin=>comp_orbmap_idet(:,ispin)
  comp_orbmap_jdet_ispin=>comp_orbmap_jdet(:,ispin)
  ! Compare orbitals in same position.
  do iorb=1,nele(ispin)
   if(comp_orbmap_idet_ispin(iorb)/=comp_orbmap_jdet_ispin(iorb))exit
  enddo ! iorb
  ! If there is no difference for this spin move on to next.
  if(iorb>nele(ispin))cycle
  ! If there already was a difference, the terms cannot be mixed.
  if(diff_ispin>0)then
   diff_ispin=0
   return
  endif
  ! Flag the difference.
  diff_ispin=ispin
  ! Determine how to match up rest of determinant.  NB, orbital indices
  ! are sorted.
  if(iorb==nele(ispin))then
   ! E.g.  [1,2,3,4,5] and [1,2,3,4,6] -> diff at 5th position, no more checks
   diff_iorb=iorb
   diff_jorb=iorb
   i=iorb
  elseif(comp_orbmap_idet_ispin(iorb+1)==comp_orbmap_jdet_ispin(iorb+1))then
   ! E.g. [1,2,4,5,6] and [1,3,4,5,6] -> diff at 2nd position, check (4:5)
   diff_iorb=iorb
   diff_jorb=iorb
   i=iorb+1
  elseif(comp_orbmap_idet_ispin(iorb+1)==comp_orbmap_jdet_ispin(iorb))then
   ! E.g. [1,2,3,4,6] and [1,3,4,5,6] -> diff at 2nd position for first
   ! det -> diff at 4th position for second det -> check (5:5)
   diff_iorb=iorb
   do i=iorb,nele(ispin)-1
    if(comp_orbmap_idet_ispin(i+1)/=comp_orbmap_jdet_ispin(i))exit
   enddo ! i
   diff_jorb=i
  elseif(comp_orbmap_idet_ispin(iorb)==comp_orbmap_jdet_ispin(iorb+1))then
   ! E.g. [1,3,4,5,6] and [1,2,3,4,6] -> diff at 2nd position for second
   ! det -> diff at 4th position for first det -> check (5:5)
   diff_jorb=iorb
   do i=iorb,nele(ispin)-1
    if(comp_orbmap_idet_ispin(i)/=comp_orbmap_jdet_ispin(i+1))exit
   enddo ! i
   diff_iorb=i
  else
   ! E.g.  [1,2,3,4,5] and [1,3,5,6,7] -> more than one diff
   diff_ispin=0
   return
  endif
  ! Run rest of check.
  do iorb=i+1,nele(ispin)
   if(comp_orbmap_idet_ispin(iorb)/=comp_orbmap_jdet_ispin(iorb))then
    diff_ispin=0
    return
   endif
  enddo ! iorb
 enddo ! ispin

 ! Return if determinants are equal.
 if(diff_ispin<1)return

 ! Compute relative sign of permutation and flag it in diff_jorb.
 if(mod(abs(diff_iorb-diff_jorb),2)==1)diff_jorb=-diff_jorb

 END SUBROUTINE compare_pair


 SUBROUTINE construct_u0star(orig,dedup,comp,opset_w0,opset_u0star,&
  &comp_u0star,orbpool)
 !----------------------------------------------------------!
 ! Compare operations IOP and JOP in OPSET belonging in set !
 ! W^(0) to see if their results can be combined to produce !
 ! operations in U^(1).  If so, create the necessary        !
 ! operations in U^(0)*.  The construction of U^(1) is not  !
 ! done here; it's deferred so that the generic procedure   !
 ! to construct U^(n) from U^(n-1) can be used (using       !
 ! U^(0)* as input).                                        !
 !                                                          !
 ! NB, U^(0)* is a partial enumeration set containing only  !
 ! those elements in U^(0) that are required to construct   !
 ! U^(1).  Here we may add unnecessary terms to U^(0)*, but !
 ! these get rejected in CONSTRUCT_UN.                      !
 !----------------------------------------------------------!
 IMPLICIT NONE
 TYPE(original),POINTER :: orig
 TYPE(deduplicated),POINTER :: dedup
 TYPE(parametrized_expansion),POINTER :: comp,comp_u0star
 TYPE(operation_set),POINTER :: opset_w0,opset_u0star
 TYPE(orbital_pool),POINTER :: orbpool
 INTEGER iop,jop,ispin,iorb,jorb,idet,jdet,n,ie,je,diff_iorb,diff_jorb,&
  &diff_ispin,i,j,k,l,nn,temp_n,temp_n_label,id,jd,ik,jk,kk,lk,nk,kop
 INTEGER &
  &temp_idet_iop(opset_w0%max_ndet),temp_idet_jop(opset_w0%max_ndet),&
  &temp_iorb_iop(opset_w0%max_ndet),temp_iorb_jop(opset_w0%max_ndet),&
  &temp_iorb_common(opset_w0%max_ndet),&
  &nmixcoeff1(opset_w0%max_ndet),&
  &imixcoeff_num1(2*opset_w0%max_ndet,opset_w0%max_ndet),&
  &imixcoeff_den1(2*opset_w0%max_ndet,opset_w0%max_ndet),&
  &mixlabel1(opset_w0%max_ndet),&
  &temp_num(opset_w0%max_ndet),temp_den(opset_w0%max_ndet),&
  &temp_num_label(opset_w0%max_ndet),temp_den_label(opset_w0%max_ndet),&
  &label_orb_indx(opset_w0%max_ndet),label_iorb(opset_w0%max_ndet),&
  &order(opset_w0%max_ndet),indices(opset_w0%max_ndet),&
  &list_idet_iop(opset_w0%max_ndet),list_iorb_iop(opset_w0%max_ndet),&
  &list_idet_jop(opset_w0%max_ndet),list_iorb_jop(opset_w0%max_ndet)
 INTEGER,POINTER :: temp_orbmap_iop(:,:)=>null(),temp_orbmap_jop(:,:)=>null()
 LOGICAL lmixsign1(opset_w0%max_ndet)
 LOGICAL flip_sign,ldum,is_first,add_iop,add_jop
 REAL(dp) rquot_num,rquot_den

 ! Allocate work pointers.
 call resize_pointer((/nemax,nspin/),temp_orbmap_iop)
 call resize_pointer((/nemax,nspin/),temp_orbmap_jop)

 ! Loop over operations IOP.
 do iop=1,opset_w0%nop

  ! We are analyzing W^(0), where compressed determinants consist entirely
  ! of original orbitals except for one which is compressed.  Construct the
  ! orbmap for the result of operation IOP, replacing the single compressed
  ! orbital in it with -1 as the last element.
  ispin=opset_w0%ispin(iop)
  iorb=abs(opset_w0%iorb(1,iop))
  temp_orbmap_iop(:,:)=comp%orbmap(:,:,opset_w0%idet(1,iop))
  if(iorb<nele(ispin))temp_orbmap_iop(iorb:nele(ispin)-1,ispin)=&
   &temp_orbmap_iop(iorb+1:nele(ispin),ispin)
  temp_orbmap_iop(nele(ispin),ispin)=-1

  ! Loop over operations JOP to compare IOP against.
  do jop=iop+1,opset_w0%nop

   ! Check both operations act on the same spin; if not, cycle.
   if(ispin/=opset_w0%ispin(jop))cycle

   ! Construct the orbmap for the result of operation JOP, replacing the
   ! single compressed orbital in it with -1, effectively assuming that
   ! the compressed orbital is equal to that in the result of operation IOP.
   ! Compare the two orbmaps and if they do not differ by a single element,
   ! cycle.
   jorb=abs(opset_w0%iorb(1,jop))
   temp_orbmap_jop(:,:)=comp%orbmap(:,:,opset_w0%idet(1,jop))
   if(jorb<nele(ispin))temp_orbmap_jop(jorb:nele(ispin)-1,ispin)=&
    &temp_orbmap_jop(jorb+1:nele(ispin),ispin)
   temp_orbmap_jop(nele(ispin),ispin)=-1
   call compare_pair(temp_orbmap_iop,temp_orbmap_jop,diff_iorb,diff_jorb,&
    &diff_ispin)
   if(diff_ispin==0)cycle

   ! Now look at the compressed orbital in each of the two compressed
   ! determinants and see which original orbitals they have in common,
   ! and cycle if they have less than two orbitals in common.
   n=0
   do idet=1,opset_w0%ndet(iop)
    ie=opset_w0%iorb(idet,iop)
    id=opset_w0%idet(idet,iop)
    iorb=comp%orbmap(abs(ie),ispin,id)
    do jdet=1,opset_w0%ndet(jop)
     je=opset_w0%iorb(jdet,jop)
     jd=opset_w0%idet(jdet,jop)
     jorb=comp%orbmap(abs(je),ispin,jd)
     if(iorb==jorb)then
      n=n+1
      temp_iorb_common(n)=iorb
      temp_idet_iop(n)=id
      temp_idet_jop(n)=jd
      temp_iorb_iop(n)=ie
      temp_iorb_jop(n)=je
      exit
     endif
    enddo ! jdet
   enddo ! idet
   ! Cycle if there are less than two common orbitals.
   if(n<2)cycle

   ! We now inspect the set of common orbitals and generate a list of all
   ! subsets of it where the orbital coefficients are proportional.
   ! For the above list of common original orbitals, compute the ratios
   ! of each of their respective coefficients.
   do i=1,n
    call divide_quot(&
     &comp%ndetcoef(temp_idet_iop(i)),comp%idetcoef_num(:,temp_idet_iop(i)),&
     &comp%idetcoef_den(:,temp_idet_iop(i)),&
     &comp%ndetcoef(temp_idet_jop(i)),comp%idetcoef_num(:,temp_idet_jop(i)),&
     &comp%idetcoef_den(:,temp_idet_jop(i)),&
     &nmixcoeff1(i),imixcoeff_num1(1,i),imixcoeff_den1(1,i),lmixsign1(i))
   enddo ! i

   ! FIXME - the code in add_to_orbpool might be more efficient, consider
   ! using it instead of the following.

   ! Figure out if two or more of the above coefficients are equal, and assign
   ! a label to each set of equal coefficients.
   do i=1,n
    mixlabel1(i)=i
   enddo ! i
   do i=1,n
    do j=i+1,n
     ! Skip if already flagged as the same.
     if(mixlabel1(i)==mixlabel1(j))cycle
     ! Compare coefficients i & j.  If they are the same, change the label
     ! of j to match that of i.
     call divide_quot(&
      &nmixcoeff1(i),imixcoeff_num1(1,i),imixcoeff_den1(1,i),&
      &nmixcoeff1(j),imixcoeff_num1(1,j),imixcoeff_den1(1,j),&
      &temp_n,temp_num,temp_den,flip_sign)
     if(lmixsign1(i))flip_sign=.not.flip_sign
     if(lmixsign1(j))flip_sign=.not.flip_sign
     if(temp_n>0)then
      ! Translate coefficients to labels.
      temp_n_label=temp_n
      do l=1,temp_n
       k=temp_num(l)
       if(k==0)then
        temp_num_label(l)=0
       elseif(k>0)then
        temp_num_label(l)=dedup%idetlabel_eff(k)
       else
        temp_num_label(l)=-dedup%idetlabel_eff(-k)
       endif
       k=temp_den(l)
       if(k==0)then
        temp_den_label(l)=0
       elseif(k>0)then
        temp_den_label(l)=dedup%idetlabel_eff(k)
       else
        temp_den_label(l)=-dedup%idetlabel_eff(-k)
       endif
      enddo ! l
      ! Simplify, and discard if labels do not cancel out (i.e., if
      ! simplification leaves two or more different labels).
      call simplify_quot(temp_n_label,temp_num_label,temp_den_label,ldum)
      if(temp_n_label>0)then
       if(count(temp_num_label(1:temp_n_label)>0)+&
        &count(temp_den_label(1:temp_n_label)>0)>1)cycle
      endif
      ! Finally, check that the numerical value of the quotient is one.
      rquot_num=1.d0
      rquot_den=1.d0
      do l=1,temp_n
       rquot_num=rquot_num*eval_dedup_detcoef(orig,dedup,temp_num(l))
       rquot_den=rquot_den*eval_dedup_detcoef(orig,dedup,temp_den(l))
      enddo ! l
      if(flip_sign)rquot_num=-rquot_num
      if(.not.compare_numbers(rquot_num,rquot_den))cycle
     endif
     where(mixlabel1(:)==mixlabel1(j))mixlabel1(:)=mixlabel1(i)
    enddo ! j
   enddo ! i

   ! Count how many sets of size >1 there are.
   nn=0
   do i=1,n
    if(count(mixlabel1==i)<2)then
     where(mixlabel1==i)mixlabel1=0
    else
     nn=nn+1
     if(i/=nn)then
      where(mixlabel1==i)mixlabel1=nn
     endif
    endif
   enddo ! i
   ! Cycle if there are no such sets.
   if(nn==0)cycle

   ! Create all required elements in U^(0)* (as operations on de-duplicated
   ! determinants).
   ! Loop over labels
   do i=1,nn

    ! Gather orbitals in this label
    k=0
    do j=1,n
     if(mixlabel1(j)==i)then
      k=k+1
      label_orb_indx(k)=j
      label_iorb(k)=temp_iorb_common(j)
     endif
    enddo ! j

    ! Loop over all possible combinations of more than one orbital from the
    ! k orbitals
    order(1:k)=2
    is_first=.true.
    indices(1:k)=0
    do while(iterate_indices_multi(k,order,indices,is_first))
     ! Make packed list of the NK values of LABEL_ORB_INDX pointing at the
     ! NK orbitals we are going to mix.
     nk=0
     do j=1,k
      if(indices(j)==1)cycle
      nk=nk+1
      kk=label_orb_indx(j)
      list_idet_iop(nk)=temp_idet_iop(kk)
      list_iorb_iop(nk)=temp_iorb_iop(kk)
      list_idet_jop(nk)=temp_idet_jop(kk)
      list_iorb_jop(nk)=temp_iorb_jop(kk)
     enddo ! j
     if(nk<2)cycle ! need at least two orbitals
     ! Sort operations by IDET.
     do ik=1,nk
      kk=ik ; lk=ik
      do jk=ik+1,nk
       if(list_idet_iop(jk)<list_idet_iop(kk))kk=jk
       if(list_idet_jop(jk)<list_idet_jop(lk))lk=jk
      enddo ! jk
      if(kk>ik)then
       call iswap1(list_idet_iop(ik),list_idet_iop(kk))
       call iswap1(list_iorb_iop(ik),list_iorb_iop(kk))
      endif
      if(lk>ik)then
       call iswap1(list_idet_jop(ik),list_idet_jop(lk))
       call iswap1(list_iorb_jop(ik),list_iorb_jop(lk))
      endif
     enddo ! ik
     ! Check whether operations already exist.
     add_iop=.true. ; add_jop=.true.
     if(associated(opset_u0star))then
      do kop=1,opset_u0star%nop
       if(opset_u0star%ndet(kop)/=nk)cycle
       if(add_iop)then
        if(all(opset_u0star%idet(1:nk,kop)==list_idet_iop(1:nk)))&
         &add_iop=.false.
       endif
       if(add_jop)then
        if(all(opset_u0star%idet(1:nk,kop)==list_idet_jop(1:nk)))&
         &add_jop=.false.
       endif
       if(.not.add_iop.and..not.add_jop)exit
      enddo ! kop
     endif
     if(.not.add_iop.and..not.add_jop)cycle
     ! Create operations in U^(0)*.
     if(add_iop)then
      call add_op(opset_u0star,nk,list_idet_iop)
      call apply_op(nk,list_idet_iop,ispin,list_iorb_iop,orig,dedup,comp,&
       &comp_u0star,orbpool)
     endif
     if(add_jop)then
      call add_op(opset_u0star,nk,list_idet_jop)
      call apply_op(nk,list_idet_jop,ispin,list_iorb_jop,orig,dedup,comp,&
       &comp_u0star,orbpool)
     endif
    enddo ! indices

   enddo ! i

  enddo ! jop

 enddo ! iop

 ! Clean up.
 deallocate(temp_orbmap_iop,temp_orbmap_jop)
 nullify(temp_orbmap_iop,temp_orbmap_jop)

 END SUBROUTINE construct_u0star


 SUBROUTINE construct_un(iop1,iop2,dedup,opset_in,comp_in,opset_out,comp_out,&
  &orbpool)
 !--------------------------------------------------------!
 ! Construct U^(n) in OPSET_OUT from U^(n-1) in OPSET_IN, !
 ! using the orbmaps in COMP_IN and adding new orbmaps to !
 ! COMP_OUT and orbitals to ORBPOOL.                      !
 !--------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: iop1,iop2
 TYPE(deduplicated),POINTER :: dedup
 TYPE(operation_set),POINTER :: opset_in,opset_out
 TYPE(parametrized_expansion),POINTER :: comp_in,comp_out
 TYPE(orbital_pool),POINTER :: orbpool
 INTEGER nop0,iop,jop,iorb,jorb,ispin,idet,jdet,ielem,nelem,nop,ndet,&
  &ncombine,icombine,orig_norb
 INTEGER elem_list(iop2-iop1+1),combine_jop(iop2-iop1+1),&
  &combine_ispin(iop2-iop1+1),combine_iorb(iop2-iop1+1),&
  &combine_jorb(iop2-iop1+1),order(iop2-iop1+1),indices(iop2-iop1+1),&
  &iop_vector(iop2-iop1+1),iorb_vector(iop2-iop1+1),&
  &idet_vector(dedup%ndet)
 LOGICAL is_first,mask(dedup%ndet),aliased_opsets
 ! Sections of components of OPSET_IN for faster access.
 INTEGER,POINTER :: &
  &opset_in_ndet_iop=>null(),opset_in_idet_iop(:)=>null(),&
  &opset_in_ndet_jop=>null(),opset_in_idet_jop(:)=>null()
 ! Sections of components of COMP_IN for faster access.
 INTEGER,POINTER :: &
  &comp_in_orbmap_iop(:,:)=>null(),comp_in_orbmap_jop(:,:)=>null()

 ! Deal with different pointers to the same targets in the arguments.
 aliased_opsets=associated(opset_in,opset_out)

 ! Initialize mask and pointers to components of TYPEd pointers.
 mask(1:dedup%ndet)=.false.
 orig_norb=orig%norb

 ! Record current number of operations in OPSET_OUT.
 nop0=0
 if(associated(opset_out))nop0=opset_out%nop

 ! Loop over pairs of operations in U^(n-1).
 do iop=iop1,iop2-1

  ! Initialize vector storing combinations of operations.
  ncombine=1
  combine_jop(ncombine)=iop
  ! Point pointers.
  opset_in_ndet_iop=>opset_in%ndet(iop)
  opset_in_idet_iop=>opset_in%idet(:,iop)
  comp_in_orbmap_iop=>comp_in%orbmap(:,:,iop)

  ! Check whether operation IOP combines with operation JOP.
  do jop=iop+1,iop2

   ! Point pointers.
   opset_in_ndet_jop=>opset_in%ndet(jop)

   ! Check if IOP and JOP operate on the same number of determinants,
   ! which is a necessary condition for them to be combinable.
   if(opset_in_ndet_iop/=opset_in_ndet_jop)cycle

   ! Point pointers.
   opset_in_idet_jop=>opset_in%idet(:,jop)
   comp_in_orbmap_jop=>comp_in%orbmap(:,:,jop)

   ! Check if IOP and JOP operate on any of the same de-duplicated
   ! determinants, in which case the operations are incompatible.
   do idet=1,opset_in_ndet_iop
    do jdet=1,opset_in_ndet_jop
     if(opset_in_idet_iop(idet)==opset_in_idet_jop(jdet))exit
    enddo ! jdet
    if(jdet<=opset_in_ndet_jop)exit
   enddo ! idet
   if(idet<=opset_in_ndet_iop)cycle

   ! IOP and JOP should yield compressed determinants with:
   ! * n compressed orbitals which are equal.
   ! * NETOT - n - 1 original orbitals which are equal.
   ! * 1 original orbital which is different.
   ! It suffices to check that just one orbital differs and that it is an
   ! original orbital.
   call compare_pair(comp_in_orbmap_iop,comp_in_orbmap_jop,iorb,jorb,ispin)
   if(ispin/=0)then
    ! Check that this is indeed an orbital from the original set, else
    ! we will generate redundant operations already in U^(n-1) and
    ! place them into U^(n).
    if(comp_in%orbmap(abs(iorb),ispin,iop)>=orig_norb.or.&
     &comp_in%orbmap(abs(jorb),ispin,jop)>=orig_norb)cycle
    ! Record which spin and orbital we are supposed to mix.
    ncombine=ncombine+1
    combine_jop(ncombine)=jop
    combine_ispin(ncombine)=ispin
    combine_iorb(ncombine)=iorb
    combine_jorb(ncombine)=jorb
   endif

  enddo ! jop

  ! Add all possible combinations of just-found operations to OPSET.
  do ispin=1,nspin
   do iorb=1,nele(ispin)
    nelem=1
    elem_list(nelem)=1
    combine_ispin(1)=ispin
    combine_iorb(1)=iorb
    combine_jorb(1)=iorb
    do icombine=2,ncombine
     if(combine_ispin(icombine)==ispin.and.&
      &abs(combine_iorb(icombine))==iorb)then
      nelem=nelem+1
      elem_list(nelem)=icombine
     endif
    enddo ! jop
    if(nelem<2)cycle

    ! Loop over combinations of the found operations
    order(1:nelem)=2
    is_first=.false.
    indices(1)=2 ; indices(2:nelem)=1
    do while(iterate_indices_multi(nelem,order,indices,is_first))
     ! Check for determinants being used multiple times.  Note no element
     ! clashes with the first since we've imposed that in the loop above,
     ! so we start checking at the second element.
     check_loop: do ielem=2,nelem
      if(indices(ielem)==1)cycle
      jop=combine_jop(elem_list(ielem))
      ! Point pointers.
      opset_in_ndet_jop=>opset_in%ndet(jop)
      opset_in_idet_jop=>opset_in%idet(:,jop)
      ! Flag used determinants, exit if already used.
      do idet=1,opset_in_ndet_jop
       if(mask(opset_in_idet_jop(idet)))exit check_loop
       mask(opset_in_idet_jop(idet))=.true.
      enddo ! idet
     enddo check_loop ! ielem
     ! Undo the MASK setting so we don't need to set the whole MASK to .false.
     ! again - it's expensive because it's a very large vector.
     uncheck_loop: do ielem=2,nelem
      if(indices(ielem)==1)cycle
      jop=combine_jop(elem_list(ielem))
      ! Point pointers.
      opset_in_ndet_jop=>opset_in%ndet(jop)
      opset_in_idet_jop=>opset_in%idet(:,jop)
      ! Unflag used determinants, exit if already unflagged.
      do idet=1,opset_in_ndet_jop
       if(.not.mask(opset_in_idet_jop(idet)))exit uncheck_loop
       mask(opset_in_idet_jop(idet))=.false.
      enddo ! idet
     enddo uncheck_loop
     if(ielem<=nelem)then
      ! Determinants have been used multiple times, so skip all combinations
      ! of operations that include the sub-combination we have just checked
      ! by setting the indices to the right of the offending index to 2, and
      ! move on to the next index set.
      if(ielem<nelem)indices(ielem+1:nelem)=2
      cycle
     endif
     ! The combination of operations is valid, so add a new operation
     nop=0
     do ielem=1,nelem
      if(indices(ielem)==1)cycle
      nop=nop+1
      jop=combine_jop(elem_list(ielem))
      iop_vector(nop)=jop
      iorb_vector(nop)=combine_jorb(elem_list(ielem))
     enddo ! ielem
     ! Compute IDET_VECTOR.
     ndet=nop*opset_in%ndet(iop_vector(1))
     call eval_idet_from_iop(opset_in,nop,iop_vector,ndet,idet_vector)
     ! Check operation is not repeated.
     if(associated(opset_out))then
      do jop=nop0+1,opset_out%nop
       if(opset_in%ndet(jop)/=ndet)cycle
       if(all(idet_vector(1:ndet)==opset_in%idet(1:ndet,jop)))exit
      enddo ! jop
      if(jop<=opset_out%nop)cycle
     endif
     ! Add operation.
     call add_op(opset_out,ndet,idet_vector(1:ndet))
     ! Construct the orbmap for the result of this operation.
     call apply_op(nop,iop_vector(1:nop),ispin,iorb_vector(1:nop),&
      &orig,dedup,comp_in,comp_out,orbpool)
    enddo ! indices

   enddo ! iorb
  enddo ! ispin

 enddo ! iop

 ! Clean up.
 nullify(opset_in_ndet_iop,opset_in_ndet_jop,opset_in_idet_iop,&
  &opset_in_idet_jop,comp_in_orbmap_iop,comp_in_orbmap_jop)

 END SUBROUTINE construct_un


 SUBROUTINE add_op(opset,ndet,idet_vector,iorb_vector,ispin)
 !----------------------------------------------------------!
 ! Add a new operation to OPSET between terms IDET and JDET !
 ! where the differing orbitals are in determinant ISPIN at !
 ! positions IORB and JORB, respectively, where a factor    !
 ! FACTOR is required to do a potential merge.              !
 ! This routine will allocate/enlarge OPSET as needed.      !
 !----------------------------------------------------------!
 IMPLICIT NONE
 TYPE(operation_set),POINTER :: opset
 INTEGER,INTENT(in) :: ndet,idet_vector(ndet)
 INTEGER,INTENT(in),OPTIONAL :: iorb_vector(ndet),ispin
 INTEGER ialloc
 if(.not.associated(opset))then
  allocate(opset,stat=ialloc)
  call check_alloc(ialloc,'OPSET_MAKE_ROOM','contents')
  opset%max_ndet=ndet ; opset%ndet_alloc=0
  opset%nop=1         ; opset%nop_alloc=0
 else
  opset%max_ndet=max(opset%max_ndet,ndet)
  opset%nop=opset%nop+1
 endif
 call opset_make_room(opset,detail=present(iorb_vector))
 opset%ndet(opset%nop)=ndet
 opset%idet(1:ndet,opset%nop)=idet_vector(1:ndet)
 if(present(ispin))opset%ispin(opset%nop)=ispin
 if(present(iorb_vector))opset%iorb(1:ndet,opset%nop)=iorb_vector(1:ndet)
 END SUBROUTINE add_op


 SUBROUTINE eval_idet_from_iop(opset,nop,iop_vector,ndet,idet_vector)
 !--------------------------------------------------------------------!
 ! Compute the list of de-duplicated determinants IDET_VECTOR(1:NDET) !
 ! involved in a U^(n) operation combining the NOP U^(n-1) operations !
 ! IOP_VECTOR(1:NOP) as listed in OPSET.                              !
 !--------------------------------------------------------------------!
 IMPLICIT NONE
 TYPE(operation_set),POINTER :: opset
 INTEGER,INTENT(in) :: nop,iop_vector(nop),ndet
 INTEGER,INTENT(inout) :: idet_vector(ndet)
 INTEGER nn,ipos,i,j,k
 nn=ndet/nop
 ipos=0
 do i=1,nop
  idet_vector(ipos+1:ipos+nn)=opset%idet(1:nn,iop_vector(i))
  ipos=ipos+nn
 enddo ! i
 ! Sort
 do i=1,ndet
  k=i
  do j=i+1,ndet
   if(idet_vector(j)<idet_vector(k))k=j
  enddo ! k
  if(k>i)call iswap1(idet_vector(i),idet_vector(k))
 enddo ! i
 END SUBROUTINE eval_idet_from_iop


 SUBROUTINE merge_op(opset,idet,jdet,iorb,jorb,ispin)
 !----------------------------------------------------------!
 ! Given the outcome of the comparison between determinants !
 ! IDET,JDET (given by IORB,JORB,ISPIN), either:            !
 ! - enlarge *the* previously-existing compatible operation !
 !   (if it exists) with one of IDET or JDET, else          !
 ! - add a new operation with IDET and JDET                 !
 !----------------------------------------------------------!
 IMPLICIT NONE
 TYPE(operation_set),POINTER :: opset
 INTEGER,INTENT(in) :: idet,jdet,iorb,jorb,ispin
 INTEGER iop,where_idet,where_jdet,id
 ! Sections of components of OPSET for faster access.
 INTEGER,POINTER :: opset_idet_iop(:)=>null(),opset_iorb_iop(:)=>null()

 if(ispin==0)return

 ! Find whether this particular operation has already been found for other
 ! terms including either idet or jdet.
 if(associated(opset))then
  do iop=1,opset%nop
   if(ispin/=opset%ispin(iop))cycle
   opset_idet_iop=>opset%idet(:,iop)
   opset_iorb_iop=>opset%iorb(:,iop)
   where_idet=0 ; where_jdet=0
   do id=1,opset%ndet(iop)
    if(opset_idet_iop(id)==idet)then
     if(opset_iorb_iop(id)==iorb)then
      where_idet=id
     elseif(opset_iorb_iop(id)==-iorb)then
      where_idet=-id
     endif
    endif
    if(opset_idet_iop(id)==jdet)then
     if(opset_iorb_iop(id)==jorb)then
      where_jdet=id
     elseif(opset_iorb_iop(id)==-jorb)then
      where_jdet=-id
     endif
    endif
   enddo
   if(where_idet/=0.and.where_jdet/=0)then ! pairwise operation already in
    return
   elseif(where_jdet>0)then
    call extend_op(opset,iop,idet,iorb)
    nullify(opset_idet_iop,opset_iorb_iop)
    return
   elseif(where_jdet<0)then
    call extend_op(opset,iop,idet,-iorb)
    nullify(opset_idet_iop,opset_iorb_iop)
    return
   elseif(where_idet>0)then
    call extend_op(opset,iop,jdet,jorb)
    nullify(opset_idet_iop,opset_iorb_iop)
    return
   elseif(where_idet<0)then
    call extend_op(opset,iop,jdet,-jorb)
    nullify(opset_idet_iop,opset_iorb_iop)
    return
   endif
  enddo
 endif

 ! Not a previous operation, so add new one.
 call add_op(opset,2,(/idet,jdet/),(/iorb,jorb/),ispin)

 ! Clean up
 nullify(opset_idet_iop,opset_iorb_iop)

 END SUBROUTINE merge_op


 SUBROUTINE extend_op(opset,iop,idet,iorb)
 !--------------------------------------------------------------!
 ! Add a new determinant IDET to the IOP-th operation in OPSET, !
 ! where the differing orbital is IORB (which includes the sign !
 ! required for the merge).                                     !
 !--------------------------------------------------------------!
 IMPLICIT NONE
 TYPE(operation_set),POINTER :: opset
 INTEGER,INTENT(in) :: iop,idet,iorb
 opset%ndet(iop)=opset%ndet(iop)+1
 opset%max_ndet=max(opset%max_ndet,opset%ndet(iop))
 call opset_make_room(opset)
 opset%idet(opset%ndet(iop),iop)=idet
 opset%iorb(opset%ndet(iop),iop)=iorb
 END SUBROUTINE extend_op


 SUBROUTINE apply_op(op_ndet,op_idet,op_ispin,op_iorb,orig,dedup,old_parexp,&
  &new_parexp,orbpool)
 !--------------------------------------------------------!
 ! Apply operation IOP in OPSET using the determinants in !
 ! OLD_PAREXP and place the result in the relevant        !
 ! of NEW_PAREXP.  This routine can be used for           !
 ! operations in W^(0) or U^(n), but for the latter the   !
 ! determinant lists in OLD_PAREXP need to refer to       !
 ! (temporary) ORBMAPs describing the result of           !
 ! operations in U^(n).                                   !
 ! The assumption is made that the new compressed orbital !
 ! is going to be built from original orbitals (so that   !
 ! nmix=1).  This should always be the case.              !
 !--------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: op_ndet,op_idet(op_ndet),op_iorb(op_ndet),op_ispin
 TYPE(original),POINTER :: orig
 TYPE(deduplicated),POINTER :: dedup
 TYPE(parametrized_expansion),POINTER :: old_parexp,new_parexp
 TYPE(orbital_pool),POINTER :: orbpool
 ! Local variables
 INTEGER,POINTER :: new_orbmap(:,:)=>null()
 INTEGER imix1(op_ndet),nmixcoeff1(op_ndet),&
  &imixcoeff_num1(2*old_parexp%max_ndetcoef,op_ndet),&
  &imixcoeff_den1(2*old_parexp%max_ndetcoef,op_ndet),&
  &idetcoef_num1(2*old_parexp%max_ndetcoef+orbpool%max_nmixcoeff),&
  &idetcoef_den1(2*old_parexp%max_ndetcoef+orbpool%max_nmixcoeff)
 INTEGER nm,i,j,iorb,jorb,korb,lorb,idet1,ndetcoef1,iswp,ialloc
 LOGICAL flip_sign,flip_sign_orbmap,aliased_parexp

 ! Quick return.
 if(op_ndet==0)return

 ! Check if we need to be careful with the reallocation of new_parexp.
 aliased_parexp=associated(old_parexp,new_parexp)

 ! Add a determinant and initialize pointers.
 if(.not.associated(new_parexp))then
  allocate(new_parexp,stat=ialloc)
  call check_alloc(ialloc,'APPLY_OP','new_parexp')
 endif
 new_parexp%ndet=new_parexp%ndet+1
 new_parexp%max_ndetcoef=max(new_parexp%max_ndetcoef,old_parexp%max_ndetcoef)
 call parexp_make_room(new_parexp)
 if(aliased_parexp)old_parexp=>new_parexp

 ! Initialize determinant coefficient to that of first determinant in
 ! operation.
 idet1=op_idet(1)
 new_parexp%ndetcoef(new_parexp%ndet)=old_parexp%ndetcoef(idet1)
 new_parexp%idetcoef_num(1:old_parexp%ndetcoef(idet1),new_parexp%ndet)=&
  &old_parexp%idetcoef_num(1:old_parexp%ndetcoef(idet1),idet1)
 new_parexp%idetcoef_den(1:old_parexp%ndetcoef(idet1),new_parexp%ndet)=&
  &old_parexp%idetcoef_den(1:old_parexp%ndetcoef(idet1),idet1)
 if(op_iorb(1)<0)new_parexp%idetcoef_num(1,new_parexp%ndet)=&
  &-new_parexp%idetcoef_num(1,new_parexp%ndet)

 ! Gather orbital indices into imix1 while computing size of imixcoeff_*.
 do i=1,op_ndet
  jorb=old_parexp%orbmap(abs(op_iorb(i)),op_ispin,op_idet(i))
  imix1(i)=orbpool%imix(1,jorb)
  if(op_iorb(i)<0)imix1(i)=-imix1(i)
 enddo ! i

 ! Construct list of mixing coefficient ratios.
 do i=1,op_ndet
  jorb=old_parexp%orbmap(abs(op_iorb(i)),op_ispin,op_idet(i))
  call divide_quot(&
   &old_parexp%ndetcoef(op_idet(i)),old_parexp%idetcoef_num(:,op_idet(i)),&
   &old_parexp%idetcoef_den(:,op_idet(i)),&
   &old_parexp%ndetcoef(idet1),old_parexp%idetcoef_num(:,idet1),&
   &old_parexp%idetcoef_den(:,idet1),&
   &nmixcoeff1(i),imixcoeff_num1(1,i),imixcoeff_den1(1,i),flip_sign)
  if(flip_sign)imix1(i)=-imix1(i)
 enddo ! i

 ! Sort mixing data by imix1.
 do i=1,op_ndet
  iswp=i
  do j=i+1,op_ndet
   if(abs(imix1(iswp))>abs(imix1(j)))iswp=j
  enddo ! j
  if(iswp>i)then
   call iswap1(imix1(iswp),imix1(i))
   call iswap1(nmixcoeff1(iswp),nmixcoeff1(i))
   do j=1,max(nmixcoeff1(iswp),nmixcoeff1(i))
    call iswap1(imixcoeff_num1(j,iswp),imixcoeff_num1(j,i))
    call iswap1(imixcoeff_den1(j,iswp),imixcoeff_den1(j,i))
   enddo ! j
  endif
 enddo ! i

 ! Add orbital to pool of orbitals (as new or as existing).
 call add_to_orbpool(op_ndet,2*old_parexp%max_ndetcoef,imix1,nmixcoeff1,&
  &imixcoeff_num1,imixcoeff_den1,orig,dedup,orbpool,korb,ndetcoef1,&
  &idetcoef_num1,idetcoef_den1)

 ! Store determinant coefficient in new%idetcoef_*.
 nm=new_parexp%ndetcoef(new_parexp%ndet)
 new_parexp%ndetcoef(new_parexp%ndet)=nm+ndetcoef1
 new_parexp%max_ndetcoef=max(new_parexp%max_ndetcoef,&
  &new_parexp%ndetcoef(new_parexp%ndet))
 call parexp_make_room(new_parexp)
 if(aliased_parexp)old_parexp=>new_parexp
 new_parexp%idetcoef_num(nm+1:nm+ndetcoef1,new_parexp%ndet)=&
  &idetcoef_num1(1:ndetcoef1)
 new_parexp%idetcoef_den(nm+1:nm+ndetcoef1,new_parexp%ndet)=&
  &idetcoef_den1(1:ndetcoef1)
 if(korb<0)new_parexp%idetcoef_num(1,new_parexp%ndet)=&
  &-new_parexp%idetcoef_num(1,new_parexp%ndet)

 ! Produce orbmap.
 new_orbmap=>new_parexp%orbmap(:,:,new_parexp%ndet)
 new_orbmap(:,:)=old_parexp%orbmap(:,:,op_idet(1))
 new_orbmap(abs(op_iorb(1)),op_ispin)=abs(korb)

 ! Sort orbmap.
 flip_sign_orbmap=.false.
 do iorb=1,nele(op_ispin)
  lorb=iorb
  do jorb=iorb+1,nele(op_ispin)
   if(new_orbmap(jorb,op_ispin)<new_orbmap(lorb,op_ispin))lorb=jorb
  enddo ! jorb
  if(lorb/=iorb)then
   call iswap1(new_orbmap(iorb,op_ispin),new_orbmap(lorb,op_ispin))
   flip_sign_orbmap=.not.flip_sign_orbmap
  endif ! must swap entries
 enddo ! iorb
 nullify(new_orbmap)

 ! Simplify determinant coefficient.
 call simplify_quot(new_parexp%ndetcoef(new_parexp%ndet),&
  &new_parexp%idetcoef_num(:,new_parexp%ndet),&
  &new_parexp%idetcoef_den(:,new_parexp%ndet),flip_sign)
 if(flip_sign.neqv.flip_sign_orbmap)&
  &new_parexp%idetcoef_num(1,new_parexp%ndet)=&
  &-new_parexp%idetcoef_num(1,new_parexp%ndet)

 ! Recalculate max_ndetcoef after simplification.
 new_parexp%max_ndetcoef=max(new_parexp%max_ndetcoef,&
  &new_parexp%ndetcoef(new_parexp%ndet))

 END SUBROUTINE apply_op


 SUBROUTINE apply_op_un(iop,parexp_un,new_parexp)
 !---------------------------------------------------------------!
 ! Transfer the result of applying operation IOP in OPSET to the !
 ! OLD expansion into the NEW expansion, where operation IOP     !
 ! belongs in U^(n) and its resulting ORBMAP has already been    !
 ! calculated and is stored in OPSET%PAREXP as determinant IDET. !
 !---------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: iop
 TYPE(parametrized_expansion),POINTER :: parexp_un,new_parexp
 ! Local variables
 LOGICAL aliased_parexp

 ! Check if we need to be careful with the reallocation of new_parexp.
 aliased_parexp=associated(parexp_un,new_parexp)

 ! Add a determinant and initialize pointers.
 new_parexp%ndet=new_parexp%ndet+1
 new_parexp%max_ndetcoef=max(new_parexp%max_ndetcoef,parexp_un%ndetcoef(iop))
 call parexp_make_room(new_parexp)
 if(aliased_parexp)parexp_un=>new_parexp

 ! For this operation we already have the deterinant and orbital data stored.
 ! Copy the determinant data.
 new_parexp%orbmap(:,:,new_parexp%ndet)=parexp_un%orbmap(:,:,iop)
 new_parexp%ndetcoef(new_parexp%ndet)=parexp_un%ndetcoef(iop)
 call parexp_make_room(new_parexp)
 if(aliased_parexp)parexp_un=>new_parexp
 new_parexp%idetcoef_num(1:new_parexp%ndetcoef(new_parexp%ndet),&
  &new_parexp%ndet)=parexp_un%idetcoef_num(1:parexp_un%ndetcoef(iop),iop)
 new_parexp%idetcoef_den(1:new_parexp%ndetcoef(new_parexp%ndet),&
  &new_parexp%ndet)=parexp_un%idetcoef_den(1:parexp_un%ndetcoef(iop),iop)

 END SUBROUTINE apply_op_un


 SUBROUTINE add_to_orbpool(n,nn,imix1,nmixcoeff1,imixcoeff_num1,&
  &imixcoeff_den1,orig,dedup,orbpool,iorb,ndetcoef1,idetcoef_num1,&
  &idetcoef_den1)
 !---------------------------------------------------------!
 ! Given a pool of orbitals, and a new orbital, either add !
 ! it to the pool if it's not already in it, or return the !
 ! index of the existing orbital in the pool.  Also return !
 ! the determinant prefactors returned by the comparison.  !
 !---------------------------------------------------------!
 IMPLICIT NONE
 ! INTENT(in)
 INTEGER,INTENT(in) :: n,nn
 INTEGER,INTENT(in) :: imix1(n),nmixcoeff1(n),imixcoeff_num1(nn,n),&
  &imixcoeff_den1(nn,n)
 TYPE(original),POINTER :: orig
 TYPE(deduplicated),POINTER :: dedup
 ! INTENT(inout)
 INTEGER,INTENT(inout) :: iorb
 INTEGER,INTENT(inout) :: ndetcoef1
 INTEGER,INTENT(inout) :: idetcoef_num1(*),idetcoef_den1(*)
 TYPE(orbital_pool),POINTER :: orbpool
 ! Local variables
 INTEGER i,j,k,jorb
 LOGICAL same,ldum
 INTEGER nquot,nquot1,nquot_label,nquot_label1
 INTEGER &
  &iquot_num(nn+orbpool%max_nmixcoeff),&
  &iquot_den(nn+orbpool%max_nmixcoeff),&
  &iquot_num1(nn+orbpool%max_nmixcoeff),&
  &iquot_den1(nn+orbpool%max_nmixcoeff),&
  &iquot_num_label(nn+orbpool%max_nmixcoeff),&
  &iquot_den_label(nn+orbpool%max_nmixcoeff),&
  &iquot_num_label1(nn+orbpool%max_nmixcoeff),&
  &iquot_den_label1(nn+orbpool%max_nmixcoeff)
 LOGICAL quot_flip_sign,quot_flip_sign1
 REAL(dp) rquot_num,rquot_den,rquot_num1,rquot_den1

 ! Initialize output
 ndetcoef1=0

 ! Check if this combination is already present
 iorb=0
 do jorb=1,orbpool%norb
  ! Check that the same orbitals appear in both compressed orbitals.
  if(orbpool%nmix(jorb)/=n)cycle
  if(any(abs(orbpool%imix(1:n,jorb))/=abs(imix1(1:n))))cycle

  ! For each orbital, construct quotient between coefficients and simplify it.
  do i=1,n
   call divide_quot(&
    &nmixcoeff1(i),imixcoeff_num1(1,i),imixcoeff_den1(1,i),&
    &orbpool%nmixcoeff(i,jorb),orbpool%imixcoeff_num(:,i,jorb),&
    &orbpool%imixcoeff_den(:,i,jorb),&
    &nquot,iquot_num,iquot_den,quot_flip_sign)
   if(orbpool%imix(i,jorb)<0.neqv.imix1(i)<0)quot_flip_sign=&
    &.not.quot_flip_sign
   ! Translate coefficients to labels.
   nquot_label=nquot
   do j=1,nquot
    k=iquot_num(j)
    if(k==0)then
     iquot_num_label(j)=0
    elseif(k>0)then
     iquot_num_label(j)=dedup%idetlabel_eff(k)
    else
     iquot_num_label(j)=-dedup%idetlabel_eff(-k)
    endif
    k=iquot_den(j)
    if(k==0)then
     iquot_den_label(j)=0
    elseif(k>0)then
     iquot_den_label(j)=dedup%idetlabel_eff(k)
    else
     iquot_den_label(j)=-dedup%idetlabel_eff(-k)
    endif
   enddo ! j
   ! Simplify and sort labels (labels have no signs, so ignore flip_sign).
   call simplify_quot(nquot_label,iquot_num_label(:),iquot_den_label(:),ldum)
   if(i==1)then
    ! Store first set of labels.
    nquot1=nquot
    iquot_num1(1:nquot)=iquot_num(1:nquot)
    iquot_den1(1:nquot)=iquot_den(1:nquot)
    quot_flip_sign1=quot_flip_sign
    nquot_label1=nquot_label
    iquot_num_label1(1:nquot_label)=iquot_num_label(1:nquot_label)
    iquot_den_label1(1:nquot_label)=iquot_den_label(1:nquot_label)
   else
    ! Compare to first set of labels.
    if(nquot_label/=nquot_label1)exit
    if(any(iquot_num_label(1:nquot_label)/=&
     &iquot_num_label1(1:nquot_label)))exit
    if(any(iquot_den_label(1:nquot_label)/=&
     &iquot_den_label1(1:nquot_label)))exit
   endif
   ! Compute numerical value of quotient.
   rquot_num=1.d0
   rquot_den=1.d0
   do j=1,nquot
    rquot_num=rquot_num*eval_dedup_detcoef(orig,dedup,iquot_num(j))
    rquot_den=rquot_den*eval_dedup_detcoef(orig,dedup,iquot_den(j))
   enddo ! j
   if(i==1)then
    ! Store numerical value of first quotient.
    rquot_num1=rquot_num
    rquot_den1=rquot_den
   else
    ! Compare to numerical value of first quotient.
    if(quot_flip_sign1.eqv.quot_flip_sign)then
     same=compare_numbers(rquot_num1*rquot_den,rquot_num*rquot_den1)
    else
     same=compare_numbers(rquot_num1*rquot_den,-rquot_num*rquot_den1)
    endif
    if(.not.same)exit
   endif
  enddo ! i
  if(i<=n)cycle
  iorb=jorb ; exit

 enddo ! jorb

 if(iorb==0)then

  ! This orbital is new, so add it to the pool.
  iorb=orbpool%norb+1
  orbpool%norb=iorb
  orbpool%max_nmix=max(orbpool%max_nmix,n)
  if(n>0)orbpool%max_nmixcoeff=max(orbpool%max_nmixcoeff,&
   &maxval(nmixcoeff1))
  call orbpool_make_room(orbpool)
  orbpool%nmix(orbpool%norb)=n
  do i=1,n
   orbpool%imix(i,orbpool%norb)=imix1(i)
   orbpool%nmixcoeff(i,orbpool%norb)=nmixcoeff1(i)
   if(nmixcoeff1(i)>0)then
    orbpool%imixcoeff_num(1:nmixcoeff1(i),i,orbpool%norb)=&
     &imixcoeff_num1(1:nmixcoeff1(i),i)
    orbpool%imixcoeff_den(1:nmixcoeff1(i),i,orbpool%norb)=&
     &imixcoeff_den1(1:nmixcoeff1(i),i)
   endif
  enddo ! i

  ! Set unused output variables.
  ndetcoef1=0

 else

  ! This orbital is old, in which case we may need to modify the detcoef.
  ndetcoef1=nquot1
  if(ndetcoef1>0)then
   idetcoef_num1(1:ndetcoef1)=iquot_num1(1:ndetcoef1)
   idetcoef_den1(1:ndetcoef1)=iquot_den1(1:ndetcoef1)
  endif
  if(quot_flip_sign1)iorb=-iorb

 endif

 END SUBROUTINE add_to_orbpool


 SUBROUTINE purge_orbpool(orbpool,parexp)
 !--------------------------------------------------------------!
 ! Eliminate orbitals from ORBPOOL that are not used in PAREXP, !
 ! updating PAREXP%ORBMAP accordingly.                          !
 !--------------------------------------------------------------!
 IMPLICIT NONE
 TYPE(orbital_pool),POINTER :: orbpool
 TYPE(parametrized_expansion),POINTER :: parexp
 TYPE(orbital_pool),POINTER :: new_orbpool=>null()
 INTEGER idet,ispin,iorb,jorb,i,norb,nmix,max_nmix,nmixcoeff,&
  &max_nmixcoeff,ialloc,remap(orbpool%norb)
 LOGICAL is_used(orbpool%norb)

 ! Flag presence of orbitals in orbmap.
 is_used=.false.
 do idet=1,parexp%ndet
  do ispin=1,nspin
   do iorb=1,nele(ispin)
    is_used(parexp%orbmap(iorb,ispin,idet))=.true.
   enddo
  enddo ! ispin
 enddo ! idet

 ! Analyze used orbitals.
 remap=0
 norb=0 ; max_nmix=0 ; max_nmixcoeff=0
 do iorb=1,orbpool%norb
  if(.not.is_used(iorb))cycle
  norb=norb+1
  remap(iorb)=norb
  nmix=orbpool%nmix(iorb)
  max_nmix=max(max_nmix,nmix)
  do i=1,nmix
   max_nmixcoeff=max(max_nmixcoeff,orbpool%nmixcoeff(i,iorb))
  enddo ! i
 enddo ! iorb

 ! Allocate new orbital pool.
 allocate(new_orbpool,stat=ialloc)
 call check_alloc(ialloc,'PURGE_ORBPOOL','container')
 new_orbpool%norb=norb
 new_orbpool%norb_alloc=norb
 new_orbpool%max_nmix=max_nmix
 new_orbpool%nmix_alloc=max_nmix
 new_orbpool%max_nmixcoeff=max_nmixcoeff
 new_orbpool%nmixcoeff_alloc=max_nmixcoeff
 allocate(new_orbpool%nmix(norb),new_orbpool%imix(max_nmix,norb),&
  &new_orbpool%nmixcoeff(max_nmix,norb),&
  &new_orbpool%imixcoeff_num(max_nmixcoeff,max_nmix,norb),&
  &new_orbpool%imixcoeff_den(max_nmixcoeff,max_nmix,norb),&
  &stat=ialloc)
 call check_alloc(ialloc,'PURGE_ORBPOOL','contents')
 ! Explicitly initialize arrays with potentially unused elements.
 new_orbpool%imix=0
 new_orbpool%nmixcoeff=0
 new_orbpool%imixcoeff_num=0
 new_orbpool%imixcoeff_den=0

 ! Copy orbitals.
 jorb=0
 do iorb=1,orbpool%norb
  if(.not.is_used(iorb))cycle
  jorb=jorb+1
  nmix=orbpool%nmix(iorb)
  new_orbpool%nmix(jorb)=nmix
  do i=1,nmix
   new_orbpool%imix(i,jorb)=orbpool%imix(i,iorb)
   nmixcoeff=orbpool%nmixcoeff(i,iorb)
   new_orbpool%nmixcoeff(i,jorb)=nmixcoeff
   new_orbpool%imixcoeff_num(1:nmixcoeff,i,jorb)=&
    &orbpool%imixcoeff_num(1:nmixcoeff,i,iorb)
   new_orbpool%imixcoeff_den(1:nmixcoeff,i,jorb)=&
    &orbpool%imixcoeff_den(1:nmixcoeff,i,iorb)
  enddo ! i
 enddo ! iorb

 ! Delete orbpool and point at new.
 call deallocate_orbpool(orbpool)
 orbpool=>new_orbpool
 nullify(new_orbpool)

 ! Update PAREXP to refer to new orbpool.
 do idet=1,parexp%ndet
  do ispin=1,nspin
   do iorb=1,nele(ispin)
    parexp%orbmap(iorb,ispin,idet)=remap(parexp%orbmap(iorb,ispin,idet))
   enddo
  enddo ! ispin
 enddo ! idet

 END SUBROUTINE purge_orbpool


! Subroutines for dealing with analytical quotients of products


 SUBROUTINE divide_quot(n1,inum1,iden1,n2,inum2,iden2,n,inum,iden,flip_sign)
 !------------------------------------------------------------------------!
 ! Perform the division of quotients:                                     !
 !                                                                        !
 ! PROD_I=1,N C_INUM(I)   PROD_I=1,N1 C_INUM1(I)   PROD_I=1,N2 C_IDEN2(I) !
 ! -------------------- = ---------------------- * ---------------------- !
 ! PROD_I=1,N C_IDEN(I)   PROD_I=1,N1 C_IDEN1(I)   PROD_I=1,N2 C_INUM2(I) !
 !                                                                        !
 ! simplifying the result by cancelling repeated factors.                 !
 !------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n1,inum1(n1),iden1(n1),n2,inum2(n2),iden2(n2)
 INTEGER,INTENT(inout) :: n,inum(n1+n2),iden(n1+n2)
 LOGICAL,INTENT(inout) :: flip_sign
 if(n1>0.and.n2>0)then
  n=n1+n2
  inum(1:n1)=inum1(1:n1)
  iden(1:n1)=iden1(1:n1)
  inum(n1+1:n)=iden2(1:n2)
  iden(n1+1:n)=inum2(1:n2)
  call simplify_quot(n,inum,iden,flip_sign)
 elseif(n1>0)then
  n=n1
  inum(1:n1)=inum1(1:n1)
  iden(1:n1)=iden1(1:n1)
  call simplify_quot(n,inum,iden,flip_sign)
 elseif(n2>0)then
  n=n2
  inum(1:n2)=iden2(1:n2)
  iden(1:n2)=inum2(1:n2)
  call simplify_quot(n,inum,iden,flip_sign)
 else
  n=0
  flip_sign=.false.
 endif
 END SUBROUTINE divide_quot


 SUBROUTINE simplify_quot(n,inum,iden,flip_sign)
 !-----------------------------------------------------------!
 ! Simplify a quotient of coefficient products analytically, !
 !      C_INUM(1) * ... * C_INUM(N)                          !
 !   ---------------------------------                       !
 !    C_IDENOM(1) * ... * C_IDENOM(N)                        !
 ! cancelling out any repeated coefficients.  On output,     !
 ! FLIP_SIGN flags the need to flip the sign of the result.  !
 !-----------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(inout) :: n,inum(n),iden(n)
 LOGICAL,INTENT(inout) :: flip_sign
 INTEGER j,k,n0,j1,j2

 ! Initialize
 flip_sign=.false.
 if(n==0)return
 n0=n

 ! Simplify numerator and denominator
 do j=1,n0
  if(inum(j)==0)cycle
  do k=1,n0
   if(abs(inum(j))==abs(iden(k)))then
    if(inum(j)==-iden(k))flip_sign=.not.flip_sign
    inum(j)=0 ; iden(k)=0
    exit
   endif
  enddo ! k
 enddo ! j

 ! Sort numerator and denominator by coefficient index
 do j=1,n0
  j1=j ; j2=j
  do k=j+1,n0
   if(inum(k)/=0.and.(inum(j1)==0.or.abs(inum(j1))>abs(inum(k))))j1=k
   if(iden(k)/=0.and.(iden(j2)==0.or.abs(iden(j2))>abs(iden(k))))j2=k
  enddo ! k
  if(j1>j)call iswap1(inum(j1),inum(j))
  if(j2>j)call iswap1(iden(j2),iden(j))
 enddo ! j

 ! Update n
 do n=1,n0
  if(inum(n)==0.and.iden(n)==0)exit
 enddo
 n=n-1

 END SUBROUTINE simplify_quot


! OPSET manipulation


 SUBROUTINE update_opset(op_idet_in,opset,detactive)
 !-----------------------------------------------!
 ! Update operation set OPSET to account for the !
 ! application of operation IOP.                 !
 !-----------------------------------------------!
 IMPLICIT NONE
 TYPE(operation_set),POINTER :: opset
 INTEGER,POINTER :: op_idet_in(:)
 LOGICAL,INTENT(inout) :: detactive(:)
 INTEGER jop,i,op_ndet,op_idet(size(op_idet_in))
 LOGICAL recalc_max_ndet
 ! Sections of components of OPSET for faster access.
 INTEGER,POINTER :: opset_ndet_jop=>null(),opset_idet_jop(:)=>null(),&
  &opset_iorb_jop(:)=>null()

 ! Copy input vector since the pointer target gets overwritten
 ! in the loop below.
 op_ndet=size(op_idet_in)
 op_idet(1:op_ndet)=op_idet_in(1:op_ndet)

 ! Flag need for recalculation of max_ndet
 recalc_max_ndet=op_ndet==opset%max_ndet

 ! Flag determinants as inactive.
 do i=1,op_ndet
  detactive(op_idet(i))=.false.
 enddo ! i

 ! Remove operations referencing the involved terms
 do jop=1,opset%nop
  opset_ndet_jop=>opset%ndet(jop)
  if(opset_ndet_jop<1)cycle
  ! Point at sections of components of opset.
  opset_idet_jop=>opset%idet(:,jop)
  opset_iorb_jop=>opset%iorb(:,jop)
  ! Check each determinant in JOP to see if it needs to be removed.
  do i=opset_ndet_jop,1,-1
   if(any(op_idet(1:op_ndet)==opset_idet_jop(i)))then
    recalc_max_ndet=recalc_max_ndet.or.opset_ndet_jop==opset%max_ndet
    opset_idet_jop(i)=opset_idet_jop(opset_ndet_jop)
    opset_iorb_jop(i)=opset_iorb_jop(opset_ndet_jop)
    opset_ndet_jop=opset_ndet_jop-1
    if(opset_ndet_jop<1)then
     opset_idet_jop(1)=0
     opset_ndet_jop=0
    endif
    if(opset_ndet_jop==0)exit
   endif
  enddo ! i
 enddo ! jop

 ! Update max_ndet if needed.
 if(recalc_max_ndet)opset%max_ndet=maxval(opset%ndet(1:opset%nop))

 ! Clean up.
 nullify(opset_ndet_jop,opset_idet_jop,opset_iorb_jop)

 END SUBROUTINE update_opset


 LOGICAL FUNCTION compare_numbers(t1,t2)
 !-------------------------------------------------------!
 ! Compare T1 and T2 and return .true. if they are equal !
 ! within the appropriate tolerance, .false. otherwise.  !
 !-------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: t1,t2
 REAL(dp) tbig,tsmall
 if(abs(t1)<abs(t2))then
  tbig=t2 ; tsmall=t1
 else
  tbig=t1 ; tsmall=t2
 endif
 if(abs(tbig)>tol_zero)then
  compare_numbers=abs(1.d0-tsmall/tbig)<tol_relative
 else
  compare_numbers=abs(tbig-tsmall)<tol_absolute
 endif
 END FUNCTION compare_numbers


 INTEGER FUNCTION test_compression(orig,dedup,comp,orbpool,quiet)
 !-------------------------------------------------!
 ! Reconstruct the deduplicated expansion from the !
 ! compressed version and compare the two.         !
 ! The return value of the function is an error    !
 ! code with the following values:                 !
 ! - Regarded as success:                          !
 !   * 0: test passed                              !
 !   * 1: NaN encountered (due to zero coeffs in   !
 !        deduplicated expansion)                  !
 ! - Regarded as soft/hard failure depending on    !
 !   purpose of test:                              !
 !   * 2: coefficient value mismatch               !
 !   * 3: coefficient sign mismatch                !
 !   * 4: coefficients diverge                     !
 ! - Regarded as failure:                          !
 !   * 5: reconstruction larger than deduplicated  !
 !   * 6: reconstruction shorter than deduplicated !
 !   * 7: reconstruction and deduplicated have     !
 !        different determinants                   !
 !-------------------------------------------------!
 USE casl
 IMPLICIT NONE
 TYPE(original),POINTER :: orig
 TYPE(deduplicated),POINTER :: dedup
 TYPE(parametrized_expansion),POINTER :: comp
 TYPE(orbital_pool),POINTER :: orbpool
 LOGICAL,OPTIONAL :: quiet
 ! Local variables
 INTEGER temp_orbmap(nemax,nspin,dedup%ndet),list_NaN(dedup%ndet)
 LOGICAL temp_NaN(dedup%ndet),temp_Inf(dedup%ndet)
 REAL(dp) temp_detcoef(dedup%ndet)
 CHARACTER(80) tmpr1,tmpr2
 INTEGER ii,iorb,jorb,ispin,idet,jdet,indices(netot),order(netot),im,nm,i,&
  &iswp,n_NaN
 LOGICAL be_quiet,is_first,found(dedup%ndet),same,any_bad
 REAL(dp) t1,t2,xcoef_num,xcoef_den,dcoef_num,dcoef_den,dcoef

 test_compression=0

 be_quiet=.false.
 if(present(quiet))be_quiet=quiet

 ! Initialize
 temp_orbmap=0
 temp_detcoef=0.d0
 temp_NaN=.false.
 temp_Inf=.false.

 ! Reconstruct expansion
 jdet=0
 do idet=1,comp%ndet
  ! Initialize determinant coefficient
  dcoef_num=1.d0
  dcoef_den=1.d0
  do i=1,comp%ndetcoef(idet)
   dcoef_num=dcoef_num*eval_dedup_detcoef(orig,dedup,comp%idetcoef_num(i,idet))
   dcoef_den=dcoef_den*eval_dedup_detcoef(orig,dedup,comp%idetcoef_den(i,idet))
  enddo ! i
  ! Initialize iterator limits for this determinant product.
  order=0
  ii=0
  do ispin=1,nspin
   ! This determinant has nele(ispin) orbitals, each of which has
   ! orbpool%nmix(comp%orbmap(...)) contributions.
   do iorb=1,nele(ispin)
    ii=ii+1
    order(ii)=orbpool%nmix(comp%orbmap(iorb,ispin,idet))
   enddo ! ispin
  enddo ! iorb
  ! Iterate over all possible combinations of contributions to orbitals
  ! from (1,1,...,1) to (order(1),order(2),...,order(n)).
  is_first=.true.
  indices=0
  do while(iterate_indices_multi(netot,order,indices,is_first))
   jdet=jdet+1
   if(jdet>dedup%ndet)then
    test_compression=5
    return
   endif
   ii=0
   xcoef_num=dcoef_num
   xcoef_den=dcoef_den
   do ispin=1,nspin
    do iorb=1,nele(ispin)
     ii=ii+1
     jorb=comp%orbmap(iorb,ispin,idet)
     temp_orbmap(iorb,ispin,jdet)=orbpool%imix(indices(ii),jorb)
     ! Compute contribution to determinant coefficient
     nm=orbpool%nmixcoeff(indices(ii),jorb)
     do im=1,nm
      xcoef_num=xcoef_num*eval_dedup_detcoef(orig,dedup,&
       &orbpool%imixcoeff_num(im,indices(ii),jorb))
      xcoef_den=xcoef_den*eval_dedup_detcoef(orig,dedup,&
       &orbpool%imixcoeff_den(im,indices(ii),jorb))
     enddo ! im
    enddo ! iorb
   enddo ! ispin
   if(abs(xcoef_den)<tol_zero)then
    if(abs(xcoef_num)<tol_zero)then
     temp_NaN(jdet)=.true.
    else
     temp_Inf(jdet)=.true.
    endif
    temp_detcoef(jdet)=0.d0
   else
    temp_detcoef(jdet)=xcoef_num/xcoef_den
   endif
   ! Sort orbmap for determinant jdet and make it positive
   do ispin=1,nspin
    do iorb=1,nele(ispin)
     iswp=iorb
     do jorb=iorb+1,nele(ispin)
      if(abs(temp_orbmap(iswp,ispin,jdet))>abs(temp_orbmap(jorb,ispin,jdet)))&
       &iswp=jorb
     enddo ! jorb
     if(iswp>iorb)then
      call iswap1(temp_orbmap(iswp,ispin,jdet),temp_orbmap(iorb,ispin,jdet))
      temp_detcoef(jdet)=-temp_detcoef(jdet)
     endif
     if(temp_orbmap(iorb,ispin,jdet)<0)then
      temp_orbmap(iorb,ispin,jdet)=-temp_orbmap(iorb,ispin,jdet)
      temp_detcoef(jdet)=-temp_detcoef(jdet)
     endif
    enddo ! iorb
   enddo ! ispin
  enddo ! indices(:)
 enddo ! idet

 ! Error out if the number of determinants does not match.
 if(jdet<dedup%ndet)then
  test_compression=6
  if(.not.be_quiet)call wout('Reconstruction has '//trim(i2s(jdet))//' terms, &
   &deduplicated has '//trim(i2s(dedup%ndet))//'terms.')
  return
 endif

 ! Compare expansions
 any_bad=.false.
 found=.false.
 do idet=1,dedup%ndet
  dcoef=eval_dedup_detcoef(orig,dedup,idet)
  do jdet=1,dedup%ndet
   if(found(jdet))cycle
   if(all(dedup%orbmap(:,:,idet)==temp_orbmap(:,:,jdet)))then
    found(jdet)=.true.
    t1=dcoef
    t2=temp_detcoef(jdet)
    same=compare_numbers(t1,t2)
    if(temp_NaN(jdet))then
     test_compression=max(test_compression,1)
    elseif(temp_Inf(jdet))then
     test_compression=max(test_compression,4)
     if(.not.be_quiet)call wout(' Warning: reconstruction of coeff #'//&
      &trim(i2s(idet))//' diverges.')
     any_bad=.true.
    elseif(.not.same)then
     write(tmpr1,*)t1
     if(compare_numbers(t1,-t2))then
      test_compression=max(test_compression,3)
      if(.not.be_quiet)call wout(' Coeff #'//trim(i2s(idet))//' is '//&
       &trim(adjustl(tmpr1))//' but we got the sign wrong')
     else
      test_compression=max(test_compression,2)
      write(tmpr2,*)t2
      if(.not.be_quiet)call wout(' Coeff #'//trim(i2s(idet))//' is '//&
       &trim(adjustl(tmpr1))//' but we get '//trim(adjustl(tmpr2)))
     endif
     any_bad=.true.
    endif
    exit
   endif
  enddo ! jdet
  if(jdet<=dedup%ndet)cycle
  ! Flag that the deduplicated term was not found in the reconstruction and
  ! report.
  test_compression=max(test_compression,7)
  if(.not.be_quiet)then
   call wout(' Det #'//trim(i2s(idet))//' in deduplicated expansion not &
    &found in reconstruction.')
   do ispin=1,nspin
    call wout('   Spin '//trim(i2s(ispin))//': ',&
     &dedup%orbmap(1:nele(ispin),ispin,idet))
   enddo ! ispin
  endif ! .not.be_quiet
 enddo ! idet

 ! Reverse check: check that all determinants in reconstruction are in the
 ! deduplicated expansion.
 do idet=1,dedup%ndet
  if(found(idet))cycle
  call wout(' Det #'//trim(i2s(idet))//' in reconstruction not present in &
   &deduplicated expansion.')
  do ispin=1,nspin
   call wout('   Spin '//trim(i2s(ispin))//': ',&
    &temp_orbmap(1:nele(ispin),ispin,idet))
  enddo ! ispin
 enddo ! idet
 if(count(found)/=dedup%ndet)then
  call wout(trim(i2s(count(found)))//' terms in the reconstruction match, '//&
   &trim(i2s(dedup%ndet-count(found)))//' do not.')
  test_compression=7
  return
 endif ! mismatch

 ! Flag NaNs and report
 if(any(temp_NaN))then
  if(.not.be_quiet)then
   call wout(' Note: could not check coefficient values for the following &
    &determinants')
   call wout(' in the de-duplicated expansion because re-expanded coefficients &
    &were 0/0:')
   n_NaN=0
   do idet=1,dedup%ndet
    if(temp_NaN(idet))then
     n_NaN=n_NaN+1
     list_NaN(n_NaN)=idet
    endif
   enddo
   call wout('',list_NaN(1:n_NaN),rfmt='(i8)',rsep='')
  endif ! .not.be_quiet
  return
 endif ! any NaNs

 END FUNCTION test_compression


 SUBROUTINE test_coeff_proportionality(orig,dedup,comp,orbpool)
 !-----------------------------------------------------------!
 ! Test that the compressed expansion remains valid when the !
 ! original coefficients are modified in a manner consistent !
 ! with the proportionality constraints associated with the  !
 ! coefficient labels.                                       !
 !-----------------------------------------------------------!
 IMPLICIT NONE
 TYPE(original),POINTER :: orig
 TYPE(deduplicated),POINTER :: dedup
 TYPE(parametrized_expansion),POINTER :: comp
 TYPE(orbital_pool),POINTER :: orbpool
 TYPE(original),POINTER :: orig_prime=>null()
 INTEGER ialloc,itest,nlabel,ilabel,idet,ierr,nfail
 REAL(dp),ALLOCATABLE :: xfactor(:)
 REAL(dp) x
 ! Variables for built-in random number generator.
 INTEGER csize,i,dum1,dum2
 INTEGER,ALLOCATABLE :: cseed(:)
 ! Number of random tests to perform.
 INTEGER,PARAMETER :: ntest=20
 ! Proportion of soft-failed tests allowed to pass the test.
 REAL(dp),PARAMETER :: allow_rel_fail=.11d0
 ! Maximum relative amount to change coefficients by.
 REAL(dp),PARAMETER :: xchange=.5d0

 ! Nothing to check if there is only one label.
 nlabel=maxval(orig%detcoef_label)
 if(nlabel==1)return
 allocate(xfactor(nlabel),stat=ialloc)
 call check_alloc(ialloc,'TEST_COEFF_PROPORTIONALITY','xfactor')

 ! Initialize random number generator from timer.
 call random_seed(size=csize)
 allocate(cseed(csize),stat=ialloc)
 call check_alloc(ialloc,'TEST_COEFF_PROPORTIONALITY','cseed')
 do i=1,csize
  call system_clock(cseed(i),dum1,dum2)
 enddo
 call random_seed(put=cseed)

 ! Instantiate ORIG_PRIME and copy from/point at ORIG.
 allocate(orig_prime)
 orig_prime%ndet=orig%ndet
 orig_prime%norb=orig%norb
 orig_prime%orbmap=>orig%orbmap
 orig_prime%detcoef_label=>orig%detcoef_label

 ! Allocate a separate set of DETCOEFs.
 allocate(orig_prime%detcoef(orig_prime%ndet),stat=ialloc)
 call check_alloc(ialloc,'TEST_COEFF_PROPORTIONALITY','detcoef')

 ! Loop over random tests.
 nfail=0
 do itest=1,ntest
  ! Generate random modification factors.
  do ilabel=1,nlabel
   call random_number(x)
   xfactor(ilabel)=1.d0+(2.d0*x-1.d0)*xchange
  enddo ! ilabel
  ! Generate set of modified coefficients.
  do idet=1,orig_prime%ndet
   orig_prime%detcoef(idet)=orig%detcoef(idet)*&
    &xfactor(orig_prime%detcoef_label(idet))
  enddo ! ilabel
  ! Test the modified expansion.
  ierr=test_compression(orig_prime,dedup,comp,orbpool,quiet=.true.)
  if(ierr>4)then
   call errstop('TEST_COEFF_PROPORTIONALITY','Failed proportionality test: &
    &wrong expansion produced.')
  elseif(ierr>1)then
   nfail=nfail+1
  endif
 enddo ! itest

 ! Error out if number of soft-failures is greater than tolerance.
 if(dble(nfail)>dble(ntest)*allow_rel_fail)call errstop(&
  &'TEST_COEFF_PROPORTIONALITY','Failed proportionality test: &
  &too many soft failures.')

 ! Clean up.
 deallocate(orig_prime%detcoef)
 nullify(orig_prime%orbmap,orig_prime%detcoef,orig_prime%detcoef_label)
 deallocate(orig_prime)
 deallocate(xfactor)

 END SUBROUTINE test_coeff_proportionality


 REAL(dp) FUNCTION eval_dedup_detcoef(orig,dedup,icoef)
 !------------------------------------------------------!
 ! Return value of ICOEF-th coefficient in deduplicated !
 ! expansion.                                           !
 !------------------------------------------------------!
 TYPE(original),POINTER :: orig
 TYPE(deduplicated),POINTER :: dedup
 INTEGER,INTENT(in) :: icoef
 INTEGER id,nd,ic,idc
 REAL(dp) t1
 if(icoef==0)then
  eval_dedup_detcoef=1.d0
 else
  ! FIXME: verify that we indeed allow ICOEF<0, else simplify code.
  ic=abs(icoef)
  eval_dedup_detcoef=0.d0
  nd=dedup%ndetcoef(ic)
  do id=1,nd
   idc=dedup%idetcoef(id,ic)
   t1=orig%detcoef(abs(idc))
   if(idc<0)t1=-t1
   eval_dedup_detcoef=eval_dedup_detcoef+t1
  enddo ! id
  if(icoef<0)eval_dedup_detcoef=-eval_dedup_detcoef
 endif
 END FUNCTION eval_dedup_detcoef


 SUBROUTINE write_cmdet_casl(orig,dedup,comp,orbpool)
 !-----------------------------------------------------------!
 ! Write an cmdet.casl file containing all the information   !
 ! required to use a compressed multi-determinant expansion. !
 !-----------------------------------------------------------!
 USE casl
 IMPLICIT NONE
 TYPE(original),POINTER :: orig
 TYPE(deduplicated),POINTER :: dedup
 TYPE(parametrized_expansion),POINTER :: comp
 TYPE(orbital_pool),POINTER :: orbpool
 CHARACTER(512) errmsg
 INTEGER ispin,idet,iorb,i,j
 LOGICAL ignore_pop

 ! Create CASL structure.
 call wout('Writing cmdet.casl.')
 call wout()
 do while(pop_casl_context())
  continue
 enddo
 call set_casl_block('cmdet.casl:CMDET',errmsg,push=.true.)
 if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))

 ! Data for book-keeping purposes
 call set_casl_item('Title',trim(adjustl(title)),errmsg)
 if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))

 ! Write original expansion
 call set_casl_block('Original expansion',errmsg,push=.true.)
 if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
 do idet=1,orig%ndet
  call set_casl_block('Term '//trim(i2s(idet))//':Coefficient',errmsg,&
   &prefer_inline=.true.)
  if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
  call set_casl_item('Term '//trim(i2s(idet))//':Coefficient:%u1',&
   &orig%detcoef(idet),errmsg)
  if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
  call set_casl_item('Term '//trim(i2s(idet))//':Coefficient:Group',&
   &orig%detcoef_label(idet),errmsg)
  if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
  do ispin=1,nspin
   call set_casl_block('Term '//trim(i2s(idet))//':Spin '//trim(i2s(ispin)),&
    &errmsg,prefer_inline=.true.)
   if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
   do iorb=1,nele(ispin)
    call set_casl_item('Term '//trim(i2s(idet))//':Spin '//trim(i2s(ispin))//&
     &':%u'//trim(i2s(iorb)),orig%orbmap(iorb,ispin,idet),errmsg)
    if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
   enddo ! iorb
  enddo ! ispin
 enddo ! idet
 ignore_pop=pop_casl_context() ! Original expansion

 ! Write deduplicated expansion
 call set_casl_block('Deduplicated coefficients',errmsg,push=.true.)
 if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
 do idet=1,dedup%ndet
  call set_casl_block('c_'//trim(i2s(idet)),errmsg,prefer_inline=.true.)
  if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
  do i=1,dedup%ndetcoef(idet)
   call set_casl_item('c_'//trim(i2s(idet))//':%u'//trim(i2s(i)),&
    &dedup%idetcoef(i,idet),errmsg)
   if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
  enddo ! i
 enddo ! idet
 ignore_pop=pop_casl_context() ! Deduplicated coefficients

 ! Write compressed expansion
 call set_casl_block('Compressed expansion',errmsg,push=.true.)
 if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))

 ! Compressed orbital pool.
 call set_casl_block('Orbital pool',errmsg,push=.true.)
 if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
 do iorb=1,orbpool%norb
  call set_casl_block('Orbital '//trim(i2s(iorb)),errmsg,push=.true.)
  if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
  do i=1,orbpool%nmix(iorb)
   call set_casl_block('Component '//trim(i2s(i)),errmsg,prefer_inline=.true.,&
    &push=.true.)
   if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
   call set_casl_item('%u1',orbpool%imix(i,iorb),errmsg)
   if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
   do j=1,orbpool%nmixcoeff(i,iorb)
    call set_casl_item('Num:%u'//trim(i2s(j)),&
     &orbpool%imixcoeff_num(j,i,iorb),errmsg)
    if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
    call set_casl_item('Den:%u'//trim(i2s(j)),&
     &orbpool%imixcoeff_den(j,i,iorb),errmsg)
    if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
   enddo ! j
   ignore_pop=pop_casl_context() ! Component <i>
  enddo ! i
  ignore_pop=pop_casl_context() ! Orbital <iorb>
 enddo ! iorb
 ignore_pop=pop_casl_context() ! Orbital pool

 ! Compressed determinants.
 call set_casl_block('Expansion',errmsg,push=.true.)
 if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
 do idet=1,comp%ndet
  call set_casl_block('Term '//trim(i2s(idet)),errmsg)
  if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
  if(comp%ndetcoef(idet)>0)then
   call set_casl_block('Term '//trim(i2s(idet))//':Coefficient',errmsg,&
    &prefer_inline=.true.)
   if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
   do i=1,comp%ndetcoef(idet)
    call set_casl_item('Term '//trim(i2s(idet))//':Coefficient:Num:%u'//&
     &trim(i2s(i)),comp%idetcoef_num(i,idet),errmsg)
    if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
    if(i<comp%ndetcoef(idet))then
     call set_casl_item('Term '//trim(i2s(idet))//':Coefficient:Den:%u'//&
      &trim(i2s(i)),comp%idetcoef_den(i,idet),errmsg)
     if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
    endif
   enddo ! i
  endif
  do ispin=1,nspin
   call set_casl_block('Term '//trim(i2s(idet))//':Spin '//trim(i2s(ispin)),&
    &errmsg,prefer_inline=.true.)
   if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
   do iorb=1,nele(ispin)
    call set_casl_item('Term '//trim(i2s(idet))//':Spin '//trim(i2s(ispin))//&
     &':%u'//trim(i2s(iorb)),comp%orbmap(iorb,ispin,idet),errmsg)
    if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))
   enddo ! iorb
  enddo ! ispin
 enddo ! idet
 ignore_pop=pop_casl_context() ! Expansion
 ignore_pop=pop_casl_context() ! Compressed expansion

 ! Finish
 ignore_pop=pop_casl_context() ! cmdet.casl:CMDET

 ! Write the file to disk.
 call write_casl('cmdet.casl','cmdet.casl',errmsg)
 if(len_trim(errmsg)>0)call errstop('WRITE_CMDET_CASL',trim(errmsg))

 END SUBROUTINE write_cmdet_casl


! Allocators and deallocators


 SUBROUTINE parexp_make_room(parexp)
 !-----------------------------------------------------!
 ! Extend pointers in PAREXP to accomodate a change in !
 ! ndet and/or ndetcoeff.                              !
 !-----------------------------------------------------!
 IMPLICIT NONE
 TYPE(parametrized_expansion),POINTER :: parexp
 INTEGER ndet,max_ndetcoef
 LOGICAL realloc_ndet,realloc_ndetcoef
 ndet=parexp%ndet
 max_ndetcoef=parexp%max_ndetcoef
 ! Find new dimensions.
 realloc_ndet=.false.
 do while(ndet>parexp%ndet_alloc)
  realloc_ndet=.true.
  parexp%ndet_alloc=parexp%ndet_alloc+parexp_realloc_ndet
 enddo
 realloc_ndetcoef=.false.
 do while(max_ndetcoef>parexp%ndetcoef_alloc)
  realloc_ndetcoef=.true.
  parexp%ndetcoef_alloc=parexp%ndetcoef_alloc+parexp_realloc_ndetcoef
 enddo
 ! Reallocate
 if(realloc_ndet)then
  call resize_pointer((/nemax,nspin,parexp%ndet_alloc/),parexp%orbmap)
  call resize_pointer((/parexp%ndet_alloc/),parexp%ndetcoef)
 endif
 if(realloc_ndet.or.realloc_ndetcoef)then
  call resize_pointer((/parexp%ndetcoef_alloc,parexp%ndet_alloc/),&
   &parexp%idetcoef_num)
  call resize_pointer((/parexp%ndetcoef_alloc,parexp%ndet_alloc/),&
   &parexp%idetcoef_den)
 endif
 END SUBROUTINE parexp_make_room


 SUBROUTINE deallocate_parexp(parexp)
 !---------------------------------------------------!
 ! Deallocate PAREXP of type PARAMETRIZED_EXPANSION. !
 !---------------------------------------------------!
 IMPLICIT NONE
 TYPE(parametrized_expansion),POINTER :: parexp
 if(.not.associated(parexp))return
 if(associated(parexp%orbmap))deallocate(parexp%orbmap)
 if(associated(parexp%ndetcoef))deallocate(parexp%ndetcoef)
 if(associated(parexp%idetcoef_num))deallocate(parexp%idetcoef_num)
 if(associated(parexp%idetcoef_den))deallocate(parexp%idetcoef_den)
 deallocate(parexp)
 nullify(parexp)
 END SUBROUTINE deallocate_parexp


 SUBROUTINE orbpool_make_room(orbpool)
 !------------------------------------------------------!
 ! Extend pointers in ORBPOOL to accomodate a change in !
 ! norb, nmix and/or nmixcoeff.                         !
 !------------------------------------------------------!
 IMPLICIT NONE
 TYPE(orbital_pool),POINTER :: orbpool
 INTEGER norb,max_nmix,max_nmixcoeff
 LOGICAL realloc_norb,realloc_nmix,realloc_nmixcoeff
 norb=orbpool%norb
 max_nmix=orbpool%max_nmix
 max_nmixcoeff=orbpool%max_nmixcoeff
 ! Find new dimensions.
 realloc_norb=.false.
 do while(norb>orbpool%norb_alloc)
  realloc_norb=.true.
  orbpool%norb_alloc=orbpool%norb_alloc+orbpool_realloc_norb
 enddo
 realloc_nmix=.false.
 do while(max_nmix>orbpool%nmix_alloc)
  realloc_nmix=.true.
  orbpool%nmix_alloc=orbpool%nmix_alloc+orbpool_realloc_nmix
 enddo
 realloc_nmixcoeff=.false.
 do while(max_nmixcoeff>orbpool%nmixcoeff_alloc)
  realloc_nmixcoeff=.true.
  orbpool%nmixcoeff_alloc=orbpool%nmixcoeff_alloc+orbpool_realloc_nmixcoeff
 enddo
 ! Reallocate
 if(realloc_norb)then
  call resize_pointer((/orbpool%norb_alloc/),orbpool%nmix)
 endif
 if(realloc_norb.or.realloc_nmix)then
  call resize_pointer((/orbpool%nmix_alloc,orbpool%norb_alloc/),&
   &orbpool%imix)
  call resize_pointer((/orbpool%nmix_alloc,orbpool%norb_alloc/),&
   &orbpool%nmixcoeff)
 endif
 if(realloc_norb.or.realloc_nmix.or.realloc_nmixcoeff)then
  call resize_pointer((/orbpool%nmixcoeff_alloc,orbpool%nmix_alloc,&
   &orbpool%norb_alloc/),orbpool%imixcoeff_num)
  call resize_pointer((/orbpool%nmixcoeff_alloc,orbpool%nmix_alloc,&
   &orbpool%norb_alloc/),orbpool%imixcoeff_den)
 endif
 END SUBROUTINE orbpool_make_room


 SUBROUTINE deallocate_orbpool(orbpool)
 !------------------------------------------!
 ! Deallocate ORBPOOL of type ORBITAL_POOL. !
 !------------------------------------------!
 IMPLICIT NONE
 TYPE(orbital_pool),POINTER :: orbpool
 if(.not.associated(orbpool))return
 if(associated(orbpool%nmix))deallocate(orbpool%nmix)
 if(associated(orbpool%imix))deallocate(orbpool%imix)
 if(associated(orbpool%nmixcoeff))deallocate(orbpool%nmixcoeff)
 if(associated(orbpool%imixcoeff_num))deallocate(orbpool%imixcoeff_num)
 if(associated(orbpool%imixcoeff_den))deallocate(orbpool%imixcoeff_den)
 deallocate(orbpool)
 nullify(orbpool)
 END SUBROUTINE deallocate_orbpool


 SUBROUTINE opset_make_room(opset,detail)
 !-------------------------------------------------!
 ! Make sure OPSET can fit at least NOP operations !
 ! of up to NDET determinants each.                !
 !-------------------------------------------------!
 IMPLICIT NONE
 TYPE(operation_set),POINTER :: opset
 LOGICAL,INTENT(in),OPTIONAL :: detail
 INTEGER ndet,nop
 LOGICAL realloc_ndet,realloc_nop,store_detail
 store_detail=.true.
 if(present(detail))store_detail=detail
 ndet=opset%max_ndet
 nop=opset%nop
 ! Find new dimensions
 realloc_ndet=.false.
 do while(ndet>opset%ndet_alloc)
  realloc_ndet=.true.
  opset%ndet_alloc=opset%ndet_alloc+opset_realloc_ndet
 enddo
 realloc_nop=.false.
 do while(nop>opset%nop_alloc)
  realloc_nop=.true.
  opset%nop_alloc=opset%nop_alloc+opset_realloc_nop
 enddo
 ! Reallocate
 if(realloc_nop)then
  call resize_pointer((/opset%nop_alloc/),opset%ndet)
  if(store_detail)then
   call resize_pointer((/opset%nop_alloc/),opset%ispin)
  endif
 endif
 if(realloc_ndet.or.realloc_nop)then
  call resize_pointer((/opset%ndet_alloc,opset%nop_alloc/),opset%idet)
  if(store_detail)then
   call resize_pointer((/opset%ndet_alloc,opset%nop_alloc/),opset%iorb)
  endif
 endif
 END SUBROUTINE opset_make_room


 SUBROUTINE deallocate_opset(opset)
 !-----------------------------------------!
 ! Deallocate OPSET of type OPERATION_SET. !
 !-----------------------------------------!
 IMPLICIT NONE
 TYPE(operation_set),POINTER :: opset
 if(.not.associated(opset))return
 if(associated(opset%idet))deallocate(opset%idet)
 if(associated(opset%iorb))deallocate(opset%iorb)
 if(associated(opset%ispin))deallocate(opset%ispin)
 if(associated(opset%ndet))deallocate(opset%ndet)
 deallocate(opset)
 nullify(opset)
 END SUBROUTINE deallocate_opset


END PROGRAM det_compress
