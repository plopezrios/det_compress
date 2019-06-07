MODULE casl
!-----------------------------------------------------!
! CASL module                                         !
! ===========                                         !
! Module to read, store, modify and write data in the !
! the CASINO serialization language (CASL).           !
!                                                     !
! PLR 09.2008                                         !
!-----------------------------------------------------!
 IMPLICIT NONE
 PRIVATE
! Public procedures
 PUBLIC read_casl,write_casl,push_casl_context,pop_casl_context,&
  &query_casl_item,get_casl_item,check_unread_casl,delete_casl_item,&
  &set_casl_item,set_casl_block,first_unread_child,unique_casl_string
! Public constants
 PUBLIC CASL_KEYSIZE,CASL_VALSIZE,CASL_FULLKEYSIZE

! Metachar type
 TYPE metachar
  PRIVATE
  CHARACTER,POINTER :: chars(:)=>null()
 END TYPE metachar

! Type for items in a CASL file
 TYPE casl_item
  PRIVATE
  TYPE(casl_item),POINTER :: next=>null(),prev=>null(),parent=>null(),&
   &first_child=>null(),last_child=>null()
  LOGICAL :: is_block=.false.,is_implicit=.false.,is_inline=.false.,&
   &been_read=.false.
  INTEGER :: ilevel=0,indent=-1,nimplicit=0
  TYPE(metachar) name,value,unique_name,full_unique_name
! AVL binary search tree data
  TYPE(casl_item),POINTER :: avl_left=>null(),avl_right=>null(),&
   &avl_parent=>null()
  INTEGER :: avl_depth=0,avl_imbalance=0
 END TYPE casl_item

! The master data item
 TYPE(casl_item),POINTER :: casl_master=>null()

! AVL tree root item
 TYPE(casl_item),POINTER :: casl_avltree_root=>null()

! Type for linear lists of CASL items (e.g. context stack, list of
! nodes flagged for deletion, etc).
 TYPE casl_list
  TYPE(casl_item),POINTER :: item=>null()
  TYPE(casl_list),POINTER :: next=>null(),prev=>null()
 END TYPE casl_list

! Context stack
 TYPE(casl_list),POINTER :: context_stack=>null()

! Parameters
! ==========
! Default character sizes:
! - Size for key
 INTEGER,PARAMETER :: CASL_KEYSIZE=32
! - Max size for full path to key (used to refer to a few levels at most)
 INTEGER,PARAMETER :: CASL_FULLKEYSIZE=256
! - Size for scalar values
 INTEGER,PARAMETER :: CASL_VALSIZE=2048
! - Number of spaces used for indentation per nesting level when writing files
 INTEGER,PARAMETER :: CASL_TABSTOP=2
! - Number of spaces used for continuation-line indenting when writing files
 INTEGER,PARAMETER :: CASL_TABCONT=5

! Interfaces
 INTERFACE get_casl_item
  MODULE PROCEDURE get_casl_item_D,get_casl_item_R,get_casl_item_Z,&
   &get_casl_item_X,get_casl_item_I,get_casl_item_C,get_casl_item_L
 END INTERFACE
 INTERFACE set_casl_item
  MODULE PROCEDURE set_casl_item_D,set_casl_item_R,set_casl_item_Z,&
   &set_casl_item_X,set_casl_item_I,set_casl_item_C,set_casl_item_L
 END INTERFACE


CONTAINS


 SUBROUTINE initialize_casl
!--------------------------------------------------------------!
! Allocate the pointer corresponding to the container of all   !
! documents, and point the contextual reference pointer at it. !
!--------------------------------------------------------------!
 IMPLICIT NONE
 if(.not.associated(casl_master))then
  allocate(casl_master)
  casl_master%is_block=.true.
 endif
 if(.not.associated(context_stack))then
  allocate(context_stack)
  context_stack%item=>casl_master
 endif
 if(.not.associated(casl_avltree_root))casl_avltree_root=>casl_master
 END SUBROUTINE initialize_casl


 SUBROUTINE read_casl(filename,errmsg)
!-------------------------------------------------!
! Read CASL file named filaname and store it in   !
! the CASL tree. To be called once per CASL file. !
!-------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: filename
 CHARACTER(512),INTENT(inout) :: errmsg
 INTEGER ierr,ipos,ipos0,indent,iline,io,nkill
 LOGICAL file_present,in_use,item_present
 CHARACTER check_eof
 TYPE(casl_item),POINTER :: document,current_item,item
 TYPE(metachar) string,label
 TYPE virtual_line
  INTEGER :: indent=-1,iline=-1,cont_indent=-1
  TYPE(metachar) line
  TYPE(virtual_line),POINTER :: next=>null()
 END TYPE virtual_line
 TYPE(virtual_line),POINTER :: vline1,vline,current_vline
 TYPE(casl_list),POINTER :: kill_list

! Initialize.
 errmsg=''
 call initialize_casl

! Open file.
 inquire(file=trim(filename),exist=file_present)
 if(.not.file_present)return
 io=9
 in_use=.true.
 do while(in_use)
  io=io+1
  inquire(io,opened=in_use)
  if(io>=99)then
   errmsg='Could not find free i/o unit.'
   return
  endif
 enddo
 open(unit=io,file=trim(filename),status='old',action='read',iostat=ierr)
 if(ierr/=0)then
  errmsg='Problem opening '//trim(filename)//'.'
  return
 endif

! Pass 1: load entire file in linked list of strings. Undo line folding,
! remove comments, trailing blanks and empty lines while at it.
 nullify(vline1,current_vline,vline)
 iline=0
 do
! Read line.
  iline=iline+1
  read(io,*,iostat=ierr)check_eof
  if(ierr/=0)exit
  backspace(io)
  call read_metachar(io,string)
! Replace tabs with spaces.
  do ipos=1,len_metachar(string)
   if(string%chars(ipos)==char(9))string%chars(ipos)=' '
  enddo ! ipos
! Remove comments, and if the line is empty, skip it.
  ipos=scan_metachar(string,'#')
  if(ipos>0)call cut_metachar(string,1,ipos-1)
  if(len_trim_metachar(string)==0)cycle
! Get indentation level and remove indentation.
  do indent=0,len_metachar(string)-1
   if(string%chars(indent+1)/=' ')exit
  enddo
  call unpad_metachar(string)
! Determine whether this is a continuation line of the previous line.
  if(associated(vline1))then
   if(indent>current_vline%indent)then
    ipos=syntactical_parse(char_metachar(current_vline%line,&
     &len_metachar(current_vline%line)),0,look_for=':')
    if(ipos<-1)then ! syntax error in line
     errmsg='Parse pass 1 problem: syntax error at line '//trim(i2s(iline))//&
      &' of file '//trim(filename)//'.'
     call unset_metachar(string)
     call clean_for_abort
     close(io)
     return
    elseif(ipos/=len_trim_metachar(current_vline%line))then
! Candidate for continuation line
     if(current_vline%cont_indent==-1.or.current_vline%cont_indent==indent)then
! Append to previous line, turning the newline into a space.
      current_vline%cont_indent=indent
      call append_char_to_metachar(current_vline%line,' ')
      call append_metachar_to_metachar(current_vline%line,string)
      call unset_metachar(string)
      cycle
     else ! error
      errmsg='Parse pass 1 problem: continuation line at line '//&
       &trim(i2s(iline))//' of file '//trim(filename)//' has different &
       &indentation from previous continuation line(s) of the same line.'
      call unset_metachar(string)
      call clean_for_abort
      close(io)
      return
     endif ! correct continuation indentation or not
    endif ! previous line does not end in a (syntactically correct) ':'
   endif ! greater indentation than previous line
  endif ! not first line
! Store the contents of the new line.
  allocate(vline)
  call copy_metachar(string,vline%line)
  call unset_metachar(string)
  vline%iline=iline
  vline%indent=indent
! Link into list.
  if(.not.associated(vline1))then ! first line
   vline1=>vline
  else
   current_vline%next=>vline
  endif
  vline%next=>vline1
  current_vline=>vline
  nullify(vline)
 enddo
 nullify(current_vline)

! Close file.
 close(io)

! Return if file is empty (or contains only comments).
 if(.not.associated(vline1))return

! Prepare document.
 allocate(document)
 document%is_block=.true.
 call copy_char_to_metachar(trim(filename),document%name)
 call copy_char_to_metachar(trim(filename),document%unique_name)
 call copy_char_to_metachar(trim(filename),document%full_unique_name)

! Insert document in master.
 document%parent=>casl_master
 if(.not.associated(casl_master%first_child))then
! Insert as first child.
  casl_master%first_child=>document
  casl_master%last_child=>document
 else
! Insert as last child.
  casl_master%last_child%next=>document
  document%prev=>casl_master%last_child
  casl_master%last_child=>document
 endif
 call insert_avltree_item(document,item_present)
 if(item_present)then
  errmsg='File '//trim(filename)//' read twice.'
  call clean_for_abort
  return
 endif

! Prepare kill_list
 nullify(kill_list)

! Set context for labeling items.
 call push_casl_context(':'//trim(filename))

! Pass 2: parse data.
 nkill=0
 vline=>vline1
 current_item=>document
 do

! Create a new item.
  allocate(item)
  item%indent=vline%indent

! Check indentation and link into list.
  if(item%indent>current_item%indent)then

! ITEM is the first child of CURRENT_ITEM.
   current_item%first_child=>item
   current_item%last_child=>item
   item%parent=>current_item
   item%ilevel=current_item%ilevel+1

  else

! ITEM is a sibling of an item at the level of CURRENT_ITEM or higher.
! Move CURRENT_ITEM up the tree until it has no more indentation than ITEM.
   do while(current_item%indent>item%indent)
    if(.not.associated(current_item%parent))exit
    current_item=>current_item%parent
   enddo
! Check that indentations of CURRENT_ITEM and ITEM are the same.
   if(current_item%indent/=item%indent)then
    errmsg='Bad indentation at line '//trim(i2s(vline%iline))//' of file '//&
     &trim(filename)//': no matching sibling.'
    deallocate(item)
    call clean_for_abort
    return
   endif

! ITEM is now the next sibling of CURRENT_ITEM.
   current_item%next=>item
   item%prev=>current_item
   item%parent=>current_item%parent
   item%parent%last_child=>item
   item%ilevel=current_item%ilevel

  endif ! child or sibling

! Update CURRENT_ITEM to point at ITEM.
  current_item=>item

! Determine if there is a colon on this line. Must be followed by a
! space or newline to have a syntactical meaning.
  ipos0=0
  do
   ipos=syntactical_parse(char_metachar(vline%line,len_metachar(vline%line)),&
    &ipos0,look_for=':')
   select case(ipos)
    case(-3)
     errmsg='Internal error.'
     call clean_for_abort
     return
    case(-2)
     errmsg='Mismatched syntactical delimiters ] ) } found when looking for &
      &indicator character in line '//trim(i2s(vline%iline))//' of file '//&
      &trim(filename)//'.'
     call clean_for_abort
     return
    case(-1)
     errmsg='Mismatched syntactical delimiters [ ( { " found when looking for &
      &indicator character in line '//trim(i2s(vline%iline))//' of file '//&
      &trim(filename)//'.'
     call clean_for_abort
     return
    case(0)
     exit
   end select
   if(ipos==len_trim_metachar(vline%line))exit ! ':<newline>'
   if(vline%line%chars(ipos+1)==' ')exit ! ':<space>'
   ipos0=ipos
  enddo

  if(ipos>0)then

! There is ':' on this line, so extract key name.
   call copy_metachar(vline%line,current_item%name)
   call cut_metachar(current_item%name,1,ipos-1)
   call unpad_metachar(current_item%name)

! Check item name is valid.  The presence of ':' in an item's name is
! forbidden, since we use ':' to divide up the fully classified item name
! into crumbs.
   if(scan_metachar(current_item%name,':')>0)then
    errmsg='Missing space after ":" or misnamed item at line '//&
     &trim(i2s(vline%iline))//' of file '//trim(filename)//': an item name &
     &cannot contain ":".'
    call clean_for_abort
    return
   endif

! Deal with CASL directives.
   if(current_item%name%chars(1)=='%')then
    if(len_metachar(current_item%name)>1)then
     if(current_item%name%chars(2)=='!')then
! Syntactic comment: re-write name, flag for deletion.
      call cut_metachar(current_item%name,3)
      call unpad_metachar(current_item%name)
      call prepend_char_to_metachar(current_item%name,&
       &'%delete'//trim(i2s(nkill))//'-')
      if(.not.associated(kill_list))then
       allocate(kill_list)
      else
       allocate(kill_list%next)
       kill_list%next%prev=>kill_list
       kill_list=>kill_list%next
      endif
      kill_list%item=>current_item
      nkill=nkill+1
     else ! unrecognized type
      call label_casl_item(current_item,label)
      errmsg='Illegal item name "'//char_metachar(label,len_metachar(label))//&
       &'": character "%" at the beginning of the name is reserved for &
       &directives, and this name does not match any known directive.'
      call unset_metachar(label)
      call clean_for_abort
      return
     endif ! type of CASL directive
    else ! just '%'
     call label_casl_item(current_item,label)
     errmsg='Illegal item name "'//char_metachar(label,len_metachar(label))//&
      &'": character "%" at the beginning of the name is reserved for &
      &directives, and this name does not match any known directive.'
     call unset_metachar(label)
     call clean_for_abort
     return
    endif
   endif

! Compute unique name and full unique name
   call copy_char_to_metachar(trim(unique_casl_string(char_metachar(&
    &current_item%name,len_metachar(current_item%name)))),&
    &current_item%unique_name)
   call copy_metachar(current_item%unique_name,current_item%full_unique_name)
   call prepend_char_to_metachar(current_item%full_unique_name,':')
   call prepend_metachar_to_metachar(current_item%full_unique_name,&
    &current_item%parent%full_unique_name)

! Determine whether this is a block, a scalar or an empty scalar depending
! on the position of ':'.
   if(ipos==len_trim_metachar(vline%line))then
! ':' at end of line: a block may be starting here. Check next line.
    if(associated(vline%next,vline1))then
! Last line, so not a block but an empty-valued scalar.
     call unset_metachar(current_item%value)
    else
     if(vline%next%indent>vline%indent)then
! Next line has greater indentation, so this is a block.
      current_item%is_block=.true.
     else
! Next line has lesser indentation, so this is an empty-valued scalar.
      call unset_metachar(current_item%value)
     endif
    endif
   else
! This is a scalar keyword-value association (or an inline block, dealt with
! later), so get value.
    call copy_metachar(vline%line,current_item%value)
    call cut_metachar(current_item%value,ipos+1)
    call unpad_metachar(current_item%value)
   endif

  else ! no ':' on this line

! This is an implicitly-named item, so compute name and get value.
   current_item%is_implicit=.true.
   current_item%parent%nimplicit=current_item%parent%nimplicit+1
   call copy_char_to_metachar('%u'//trim(i2s(current_item%parent%nimplicit)),&
    &current_item%name)
   call copy_char_to_metachar('%u'//trim(i2s(current_item%parent%nimplicit)),&
    &current_item%unique_name)
   call copy_metachar(current_item%unique_name,current_item%full_unique_name)
   call prepend_char_to_metachar(current_item%full_unique_name,':')
   call prepend_metachar_to_metachar(current_item%full_unique_name,&
    &current_item%parent%full_unique_name)
   call copy_metachar(vline%line,current_item%value)
   call unpad_metachar(vline%line)

  endif

! Insert in AVL tree.
  call insert_avltree_item(current_item,item_present)
  if(item_present)then
   call label_casl_item(current_item,label)
   errmsg='Item '//char_metachar(label,len_metachar(label))//&
    &' found twice (second occurrence at line '//trim(i2s(vline%iline))//&
    &' of file '//trim(filename)//').  This may be due to bad indentation in &
    &file '//trim(filename)//'.'
   call clean_for_abort
   return
  endif

! Detect inline blocks (so far flagged as scalar).
  if(.not.(current_item%is_block.or.current_item%is_implicit))then
   if(len_metachar(current_item%value)>0)then
    if(current_item%value%chars(1)=='[')then
     current_item%is_block=.true.
     current_item%is_inline=.true.
     call parse_inline_block(current_item,errmsg)
     if(len_trim(errmsg)>0)then
      call clean_for_abort
      return
     endif
    endif
   endif
  endif

! Go to next line.
  if(associated(vline%next,vline1))exit
  vline=>vline%next
 enddo

! Apply '%!' directives, if any.
 call apply_kill_list

! Release memory.
 call clean_vline_storage

! Pop casl context stack.
 call pop_casl_context()


 CONTAINS


  RECURSIVE SUBROUTINE parse_inline_block(item,errmsg)
!----------------------------------------------------------------!
! Parse the value of ITEM, which is an inline block, and place   !
! the contents in the pointer tree. The inline block may contain !
! any number of inline blocks at any depth level.                !
!----------------------------------------------------------------!
  IMPLICIT NONE
  TYPE(casl_item),POINTER :: item
  CHARACTER(512),INTENT(inout) :: errmsg
  INTEGER ipos,ipos0
  LOGICAL item_present
  TYPE(casl_item),POINTER :: current_child
  TYPE(metachar) label

! Initialize.
  errmsg=''

! Get contents of brackets and do error checking.
  ipos=syntactical_parse(char_metachar(item%value,len_metachar(item%value)),1)
  select case(ipos)
   case(-3)
    errmsg='Internal error.'
    return
   case(-2)
    call label_casl_item(item,label)
    errmsg='Mismatched syntactical delimiters ] ) } found when looking for &
     &end of inline block "'//char_metachar(label,len_metachar(label))//'".'
    call unset_metachar(label)
    return
   case(-1)
    call label_casl_item(item,label)
    errmsg='Mismatched syntactical delimiters [ ( { " found when looking for &
     &end of inline block "'//char_metachar(label,len_metachar(label))//'".'
    call unset_metachar(label)
    return
   case(0)
    call label_casl_item(item,label)
    errmsg='Cannot find closing square bracket for inline block "'//&
     &char_metachar(label,len_metachar(label))//'".'
    call unset_metachar(label)
    return
  end select
  if(ipos<len_metachar(item%value))then
   call label_casl_item(item,label)
   errmsg='Trailing characters after inline block "'//&
    &char_metachar(label,len_metachar(label))//'".'
   call unset_metachar(label)
   return
  endif
  call cut_metachar(item%value,2,ipos-1)

! Empty block?
  if(len_trim_metachar(item%value)==0)return

! Loop over items in item%value.
  nullify(current_child)
  do

! Create a new item.
   if(.not.associated(current_child))then
    allocate(current_child)
    item%first_child=>current_child
    item%last_child=>current_child
   else
    allocate(current_child%next)
    current_child%next%prev=>current_child
    current_child=>current_child%next
    item%last_child=>current_child
   endif
   current_child%parent=>item
   current_child%ilevel=item%ilevel+1
   current_child%indent=item%indent+1

! Get next item in line, up to a comma followed by a space.
   ipos0=0
   do
    ipos=syntactical_parse(char_metachar(item%value,len_metachar(item%value)),&
     &ipos0,look_for=',')
    select case(ipos)
     case(-3)
      errmsg='Internal error.'
      return
     case(-2)
      call label_casl_item(item,label)
      errmsg='Mismatched syntactical delimiters ] ) } found when looking for &
       &separator in inline block "'//&
       &char_metachar(label,len_metachar(label))//'".'
      call unset_metachar(label)
      return
     case(-1)
      call label_casl_item(item,label)
      errmsg='Mismatched syntactical delimiters [ ( { " found when looking &
       &for separator in inline block "'//&
       &char_metachar(label,len_metachar(label))//'".'
      call unset_metachar(label)
      return
     case(0)
      exit
    end select
    if(ipos==0)exit
    if(ipos==len_trim_metachar(item%value))exit
    if(item%value%chars(ipos+1)==' ')exit
    ipos0=ipos
   enddo

   call copy_metachar(item%value,current_child%value)
   if(ipos<1)then
! Only one item remaining. Empty the string to flag exit.
    call unset_metachar(item%value)
   else
! More items remain. Remove the one we've just read.
    call cut_metachar(current_child%value,1,ipos-1)
    call cut_metachar(item%value,ipos+1)
   endif

! Determine if there is a colon on this line.
   ipos=scan_metachar(current_child%value,':')
   if(ipos>0)then
    call copy_metachar(current_child%value,current_child%name)
    call cut_metachar(current_child%name,1,ipos-1)
    call unpad_metachar(current_child%name)
! Check item name is valid:
! - The presence of ':' in an item's name is forbidden, since we use ':' to
!   divide up the fully classified item name into crumbs.
    if(any(current_child%name%chars(:)==':'))then
     errmsg='Missing space or misnamed item at line '//trim(i2s(vline%iline))&
      &//' of file '//trim(filename)//': an item name cannot contain ":".'
     return
    endif
    if(current_child%name%chars(1)=='%')then
! CASL directive
     if(len_metachar(current_child%name)>1)then
      if(current_child%name%chars(2)=='!')then
! Syntactic comment: re-write name, flag for deletion.
       call cut_metachar(current_child%name,3)
       call unpad_metachar(current_child%name)
       call prepend_char_to_metachar(current_child%name,&
        &'%delete'//trim(i2s(nkill))//'-')
       if(.not.associated(kill_list))then
        allocate(kill_list)
       else
        allocate(kill_list%next)
        kill_list%next%prev=>kill_list
        kill_list=>kill_list%next
       endif
       kill_list%item=>current_child
       nkill=nkill+1
      else
       call label_casl_item(current_child,label)
       errmsg='Illegal item name "'//char_metachar(label,len_metachar(label))//&
        &'": character "%" at the beginning of the name is reserved for &
        &directives, and this name does not match any known directive.'
       call unset_metachar(label)
       return
      endif ! type of CASL directive
     else ! just '%'
      call label_casl_item(current_child,label)
      errmsg='Illegal item name "'//char_metachar(label,len_metachar(label))//&
       &'": character "%" at the beginning of the name is reserved for &
       &directives, and this name does not match any known directive.'
      call unset_metachar(label)
      return
     endif
    endif
    call copy_char_to_metachar(trim(unique_casl_string(char_metachar&
     &(current_child%name,len_metachar(current_child%name)))),&
     &current_child%unique_name)
    call copy_metachar(current_child%unique_name,&
     &current_child%full_unique_name)
    call prepend_char_to_metachar(current_child%full_unique_name,':')
    call prepend_metachar_to_metachar(current_child%full_unique_name,&
     &current_child%parent%full_unique_name)
    call cut_metachar(current_child%value,ipos+1)
    call unpad_metachar(current_child%value)
   else
! This is an implicitly-named item.
    current_child%is_implicit=.true.
    item%nimplicit=item%nimplicit+1
    call copy_char_to_metachar('%u'//trim(i2s(item%nimplicit)),&
     &current_child%name)
    call copy_char_to_metachar('%u'//trim(i2s(item%nimplicit)),&
     &current_child%unique_name)
    call copy_metachar(current_child%unique_name,&
     &current_child%full_unique_name)
    call prepend_char_to_metachar(current_child%full_unique_name,':')
    call prepend_metachar_to_metachar(current_child%full_unique_name,&
     &current_child%parent%full_unique_name)
    call unpad_metachar(current_child%value)
    if(len_trim_metachar(current_child%value)==0)then
     call label_casl_item(item,label)
     errmsg='Syntax error in contents of "'//&
      &char_metachar(label,len_metachar(label))//'": unnamed scalar with &
      &zero-sized value encountered (i.e., empty item #'//&
      &trim(i2s(item%nimplicit))//' in comma-separated list).'
     call unset_metachar(label)
     return
    endif
   endif

! Insert in AVL tree.
   call insert_avltree_item(current_child,item_present)
   if(item_present)then
    call label_casl_item(current_child,label)
    if(.not.current_child%is_implicit)then
     errmsg='Item '//char_metachar(label,len_metachar(label))//' found twice.'
    else
     errmsg='Problem with generation of internal names for implicitly named &
      &items: label "'//char_metachar(label,len_metachar(label))//'" generated &
      &twice. This is a bug.'
    endif
    call unset_metachar(label)
    return
   endif

! Deal with deeper levels of inline blocks.
   if(.not.current_child%is_implicit)then
    if(len_metachar(current_child%value)>0)then
     if(current_child%value%chars(1)=='[')then
      current_child%is_block=.true.
      current_child%is_inline=.true.
      call parse_inline_block(current_child,errmsg)
      if(len_trim(errmsg)>0)return
     endif
    endif
   endif

   if(len_metachar(item%value)==0)exit

  enddo ! Loop over items in line.

! Free memory used by in-line block string once it has been parsed into
! structured data.
  call unset_metachar(item%value)

  END SUBROUTINE parse_inline_block


  RECURSIVE FUNCTION syntactical_parse(line,i0,look_for,ignore_other) &
   &RESULT(syn_parse)
!---------------------------------------------------------------------!
! Find the position of the first instance at the current syntactical  !
! level of character LOOK_FOR in LINE after position I0, or, if       !
! LOOK_FOR is the empty string, find the closing bracket correspoding !
! to the opening bracket at position I0. If IGNORE_OTHER = T, other   !
! syntax elements are ignored and not checked for correctness.        !
!                                                                     !
! Error codes:                                                        !
!  0 : syntax correct (if checked), character not found               !
! -1 : string contains an unmatched left delimiter                    !
! -2 : string contains an unmatched right delimiter                   !
! -3 : don't know what to look for (missing LOOK_FOR/wrong I0)        !
!---------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: i0
  LOGICAL,INTENT(in),OPTIONAL :: ignore_other
  CHARACTER(*),INTENT(in) :: line
  CHARACTER,INTENT(in),OPTIONAL :: look_for
  INTEGER syn_parse,i,n
  LOGICAL no_syntax
  CHARACTER left_delim,right_delim,c1
  syn_parse=-1
  no_syntax=.false. ; if(present(ignore_other))no_syntax=ignore_other
  if(.not.present(look_for))then
! If not provided, the I0-th character of LINE should be a left delimiter, and
! we are looking for the right delimiter, which we define here:
   if(i0<1)then
    syn_parse=-3 ; return
   endif
   left_delim=line(i0:i0)
   select case(left_delim)
    case('[')    ; right_delim=']'
    case('(')    ; right_delim=')'
    case('{')    ; right_delim='}'
    case('"')    ; right_delim='"'
   case default ; syn_parse=-3 ; return
   end select
  else
   right_delim=look_for
  endif
  i=i0 ; n=len(line)
  do
   i=i+1 ; if(i>n)exit
   c1=line(i:i)
   if(c1==right_delim)then
! Character found.
    syn_parse=i ; return
   elseif(.not.no_syntax)then
! Check for changes in the syntactical level depth.
    select case(c1)
     case('[','(','{') ; i=syntactical_parse(line,i)
     case('"')         ; i=syntactical_parse(line,i,ignore_other=.true.)
     case(']',')','}') ; syn_parse=-2 ; return
     case default      ; cycle ! normal character
    end select
    if(i<0)then ! Propagate whatever error was returned above.
     syn_parse=i ; return
    elseif(i==0)then ! Generate "unmatched left delimiter" error.
     syn_parse=-1 ; return
    endif
   endif
  enddo
! Not found.
  syn_parse=0
  END FUNCTION syntactical_parse


  SUBROUTINE clean_vline_storage
!------------------------------------!
! Release memory from VLINE storage. !
!------------------------------------!
  IMPLICIT NONE
  TYPE(virtual_line),POINTER :: current_vline,vline
  if(associated(vline1))then
   current_vline=>vline1%next
   do while(.not.associated(current_vline,vline1))
    vline=>current_vline%next
    call unset_metachar(current_vline%line)
    deallocate(current_vline)
    current_vline=>vline
   enddo
   call unset_metachar(vline1%line)
   deallocate(vline1)
   nullify(vline1)
  endif
  END SUBROUTINE clean_vline_storage


  SUBROUTINE apply_kill_list(clean_only)
!--------------------------------!
! Apply '%!' directives, if any. !
!--------------------------------!
  IMPLICIT NONE
  LOGICAL,OPTIONAL :: clean_only
  LOGICAL just_clean
  if(associated(kill_list))then
   just_clean=.false. ; if(present(clean_only))just_clean=clean_only
   do
    if(.not.just_clean)call kill_casl_item(kill_list%item)
    if(associated(kill_list%prev))then
     kill_list=>kill_list%prev
     deallocate(kill_list%next)
    else
     deallocate(kill_list)
     nullify(kill_list)
     exit
    endif
   enddo
  endif
  END SUBROUTINE apply_kill_list


  SUBROUTINE clean_for_abort
!--------------------------------------------------!
! Release memory before returning due to an error. !
!--------------------------------------------------!
  IMPLICIT NONE
  call clean_vline_storage
  call apply_kill_list(clean_only=.true.)
  if(associated(document))call kill_casl_item(document)
  END SUBROUTINE clean_for_abort


 END SUBROUTINE read_casl


 SUBROUTINE write_casl(label,file,errmsg)
!---------------------------------------------------------------!
! Write the contents of the CASL tree under item labelled LABEL !
! to file FILE, overwriting any existing file named FILE.       !
!---------------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label,file
 CHARACTER(512),INTENT(inout) :: errmsg
 INTEGER io,ierr,delta_indent
 LOGICAL in_use
 TYPE(casl_item),POINTER :: item
 TYPE(metachar) line

! Initialize.
 errmsg=''

! Locate root and make sure it's non-empty.
 call initialize_casl
 call find_avltree_item(label,item)
 if(.not.associated(item))return
 if(.not.associated(item%first_child))return
 if(.not.associated(item%parent,casl_master))return

! Open file.
 io=9
 in_use=.true.
 do while(in_use)
  io=io+1
  inquire(io,opened=in_use)
  if(io>=99)then
   errmsg='Could not find free i/o unit.'
   return
  endif
 enddo
 open(io,file=trim(file),status='replace',action='write',iostat=ierr)
 if(ierr/=0)then
  errmsg='Problem opening '//trim(file)//'.'
  return
 endif

! Detect base indentation.
 delta_indent=item%ilevel+1

! Write CASL tree.
 item=>item%first_child
 do
  call write_casl_item(item,io,delta_indent,line,.false.,0,.false.)
  if(.not.associated(item%next))exit
  write(io,'(a)')''
  item=>item%next
 enddo

! Close file.
 close(io)


 CONTAINS


  RECURSIVE SUBROUTINE write_casl_item(item,io,delta_indent,line,inline,&
   &base_ilevel,need_comma)
!------------------------------------------------------------------------!
! Write item ITEM to unit IO, with indentation ITEM%ILEVEL-DELTA_INDENT. !
! LINE, INLINE, BASE_ILEVEL and NEED_COMMA are passed in order to write  !
! inline-mode blocks correctly.                                          !
!------------------------------------------------------------------------!
  IMPLICIT NONE
  TYPE(casl_item),POINTER :: item
  INTEGER,INTENT(in) :: io,delta_indent,base_ilevel
  LOGICAL,INTENT(in) :: inline,need_comma
  TYPE(metachar),INTENT(inout) :: line
  TYPE(casl_item),POINTER :: scan_item

  if(.not.inline)then
! We are not in inline mode.
   call unset_metachar(line)
   if(item%is_block)then
! ITEM is a block.
    if(associated(item%first_child))then
! ITEM is a non-empty block.
     if(.not.item%is_inline)then
! We are not entering inline mode.
      write(io,'(a)')repeat(' ',CASL_TABSTOP*(item%ilevel-delta_indent))//&
       &char_metachar(item%name,len_metachar(item%name))//':'
      scan_item=>item%first_child
      do
       call write_casl_item(scan_item,io,delta_indent,line,.false.,0,.false.)
       if(.not.associated(scan_item%next))exit
       scan_item=>scan_item%next
      enddo
     else
! We are entering inline mode.
      call copy_char_to_metachar(repeat(' ',&
       &CASL_TABSTOP*(item%ilevel-delta_indent))//&
       &char_metachar(item%name,len_metachar(item%name))//': [ ',line)
      scan_item=>item%first_child
      do
       call write_casl_item(scan_item,io,delta_indent,line,.true.,item%ilevel,&
        &associated(scan_item%next))
       if(.not.associated(scan_item%next))exit
       scan_item=>scan_item%next
      enddo
      call append_char_to_metachar(line,' ]')
      write(io,'(a)')char_metachar(line,len_metachar(line))
      call unset_metachar(line)
     endif ! entering inline or not
    else
! ITEM is an empty block.
     write(io,'(a)')repeat(' ',CASL_TABSTOP*(item%ilevel-delta_indent))//&
      &char_metachar(item%name,len_metachar(item%name))//': [ ]'
    endif
   else
! ITEM is a regular item.
    if(.not.item%is_implicit)then
     write(io,'(a)')repeat(' ',CASL_TABSTOP*(item%ilevel-delta_indent))//&
      &char_metachar(item%name,len_metachar(item%name))//': '//&
      &char_metachar(item%value,len_metachar(item%value))
    else
     write(io,'(a)')repeat(' ',CASL_TABSTOP*(item%ilevel-delta_indent))//&
      &char_metachar(item%value,len_metachar(item%value))
    endif
   endif ! block or regular item
  else
! We are already in inline mode.
   if(item%is_block)then
! ITEM is a block.
    call append_or_flush(line,&
     &char_metachar(item%name,len_metachar(item%name))//': [ ',&
     &base_ilevel-delta_indent)
    if(associated(item%first_child))then
! ITEM is a non-empty block.
     scan_item=>item%first_child
     do
      call write_casl_item(scan_item,io,delta_indent,line,.true.,base_ilevel,&
       &associated(scan_item%next))
      if(.not.associated(scan_item%next))exit
      scan_item=>scan_item%next
     enddo
    endif
    if(need_comma)then
     call append_or_flush(line,' ], ',base_ilevel-delta_indent)
    else
     call append_or_flush(line,' ]',base_ilevel-delta_indent)
    endif
   else
! ITEM is a regular item.
    if(.not.item%is_implicit)then
     if(need_comma)then
      call append_or_flush(line,&
       &char_metachar(item%name,len_metachar(item%name))//': '//&
       &char_metachar(item%value,len_metachar(item%value))//', ',&
       &base_ilevel-delta_indent)
     else
      call append_or_flush(line,&
       &char_metachar(item%name,len_metachar(item%name))//': '//&
       &char_metachar(item%value,len_metachar(item%value)),&
       &base_ilevel-delta_indent)
     endif
    else
     if(need_comma)then
      call append_or_flush(line,&
       &char_metachar(item%value,len_metachar(item%value))//', ',&
       &base_ilevel-delta_indent)
     else
      call append_or_flush(line,&
       &char_metachar(item%value,len_metachar(item%value)),&
       &base_ilevel-delta_indent)
     endif
    endif
   endif
  endif
  END SUBROUTINE write_casl_item


  SUBROUTINE append_or_flush(line,string,indent)
!---------------------------------------------------------------------!
! Append STRING to LINE. If the concatenation is longer than a fixed  !
! internal parameter, LINE is flushed and reset to INDENT indentation !
! spaces before appending.                                            !
!---------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: indent
  CHARACTER(*),INTENT(in) :: string
  TYPE(metachar),INTENT(inout) :: line
  INTEGER i1,i2,len_line,len_string,len_all
  INTEGER,PARAMETER :: max_line_length=79
  CHARACTER(len(string)) adjustl_string
  len_line=len_metachar(line)
  len_string=len(string)
  len_all=len_line+len_string
  if(len_all>max_line_length)then
   write(io,'(a)')char_metachar(line,len_line)
   adjustl_string=adjustl(string)
   i1=len_trim(string)
   i2=len_trim(adjustl_string)
   if(i1>i2)then
    call copy_char_to_metachar(repeat(' ',CASL_TABSTOP*indent+CASL_TABCONT)//&
     &adjustl_string(1:len_string-i1+i2),line)
   else
    call copy_char_to_metachar(repeat(' ',CASL_TABSTOP*indent+CASL_TABCONT)//&
     &string,line)
   endif
  else
   call append_char_to_metachar(line,string)
  endif
  END SUBROUTINE append_or_flush


 END SUBROUTINE write_casl


 FUNCTION unique_casl_string(string) RESULT(unique_string)
!---------------------------------------------!
! Turn string to lowercase and remove spaces. !
!---------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: string
 INTEGER i,ic,ich
 CHARACTER(len(string)) unique_string
! Convert to lowercase.
 do i=1,len_trim(string)
  ich=ichar(string(i:i))
  if(ich>=ichar('A').and.ich<=ichar('Z'))ich=ich+ichar('a')-ichar('A')
  unique_string(i:i)=char(ich)
 enddo ! i
! Remove spaces.
 ic=scan(trim(unique_string),' ')
 do while(ic>0)
  unique_string(ic:)=trim(unique_string(ic+1:))
  ic=scan(trim(unique_string),' ')
 enddo ! ic>0
 END FUNCTION unique_casl_string


 SUBROUTINE push_casl_context(label,respect_spelling)
!---------------------------------------------------------------!
! Add LABEL at the top of the context stack and move the stack  !
! pointer to it. The function returns .false. if LABEL does not !
! exist or if it's not a block item, in which case no action is !
! taken.                                                        !
!---------------------------------------------------------------!
 IMPLICIT NONE
 LOGICAL,INTENT(in),OPTIONAL :: respect_spelling
 CHARACTER(*),INTENT(in) :: label
 LOGICAL keep_name,valid
 CHARACTER(len(label)) crumb
 TYPE(casl_item),POINTER :: item

! Get optional parameters or set defaults.
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling

! Make sure we have a CASL structure.
 call initialize_casl

! Locate item to push, must be a block.
 call find_avltree_item(label,item)
 valid=associated(item)
 if(valid)valid=item%is_block

 if(valid)then
! Found and is a block: create a new stack level and switch to it.
  allocate(context_stack%next)
  context_stack%next%prev=>context_stack
  context_stack%next%item=>item
  context_stack=>context_stack%next
! Re-spell name if requested.
  if(.not.keep_name.and.len_trim(label)>0)then
   call extract_last_crumb(trim(label),crumb)
   call copy_char_to_metachar(trim(crumb),item%name)
  endif
 endif

 END SUBROUTINE push_casl_context


 SUBROUTINE pop_casl_context()
!-----------------------------------------------------------------!
! Remove the top item in the context stack and move the stack     !
! pointer one item down. The function returns .false. if there is !
! only one item in the stack, in which case no action is taken.   !
!-----------------------------------------------------------------!
 IMPLICIT NONE
! Make sure we have a CASL structure.
 call initialize_casl
! Is this the bottom of the stack?
 if(associated(context_stack%prev))then
  context_stack=>context_stack%prev
  deallocate(context_stack%next)
  nullify(context_stack%next)
 endif
 END SUBROUTINE pop_casl_context


 SUBROUTINE split_first_crumb_metachar(label,crumb)
!--------------------------------------------------!
! Split LABEL = 'crumb1:crumb2:...:lastcrumb' into !
! LABEL='crumb2:...:lastcrumb' and crumb='crumb1'. !
!--------------------------------------------------!
 IMPLICIT NONE
 TYPE(metachar),INTENT(inout) :: label,crumb
 INTEGER ipos
 ipos=scan_metachar(label,':')
 call copy_metachar(label,crumb)
 if(ipos>0)then
  call cut_metachar(crumb,1,ipos-1)
  call cut_metachar(label,ipos+1)
 else
  call unset_metachar(label)
 endif
 END SUBROUTINE split_first_crumb_metachar


 INTEGER FUNCTION count_crumbs(label)
!-----------------------------------------!
! Count the number of crumbs n in LABEL = !
! 'crumb1:crumb2:...:crumbn'.             !
!-----------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 TYPE(metachar) part_label
 INTEGER ipos
 count_crumbs=0
 if(len_trim(label)==0)return
 count_crumbs=1
 call copy_char_to_metachar(trim(label),part_label)
 do
  ipos=scan_metachar(part_label,':')
  if(ipos<1)exit
  count_crumbs=count_crumbs+1
  call cut_metachar(part_label,ipos+1)
 enddo
 call unset_metachar(part_label)
 END FUNCTION count_crumbs


 SUBROUTINE remove_last_crumb_metachar(label)
!---------------------------------------------!
! Given LABEL = 'crumb1:crumb2:...:lastcrumb' !
! return LABEL='crumb1:crumb2:...'.           !
!---------------------------------------------!
 IMPLICIT NONE
 TYPE(metachar),INTENT(inout) :: label
 INTEGER ipos
 ipos=scan_metachar(label,':',back=.true.)
 if(ipos>0)then
  call cut_metachar(label,1,ipos-1)
 else
  call unset_metachar(label)
 endif
 END SUBROUTINE remove_last_crumb_metachar


 SUBROUTINE extract_last_crumb(label,crumb)
!-----------------------------------------------------!
! Extract lastcrumb in 'crumb1:crumb2:...:lastcrumb'. !
!-----------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 CHARACTER(len(label)),INTENT(out) :: crumb
 INTEGER ipos
 ipos=scan(label,':',back=.true.)
 if(ipos>0)then
  crumb=label(ipos+1:)
 else
  crumb=label
 endif
 END SUBROUTINE extract_last_crumb


 SUBROUTINE label_casl_item(item,label,full,unique)
!------------------------------------------------------------!
! Return the full label of ITEM with respect to CONTEXT_BASE !
! in string LABEL.                                           !
!------------------------------------------------------------!
 IMPLICIT NONE
 TYPE(casl_item),POINTER :: item
 TYPE(metachar),INTENT(inout) :: label
 LOGICAL,INTENT(in),OPTIONAL :: full,unique
 LOGICAL within_context,reduce
 TYPE(casl_item),POINTER :: scan_item
 within_context=.true. ; if(present(full))within_context=.not.full
 reduce=.false. ; if(present(unique))reduce=unique
 call unset_metachar(label)
 call initialize_casl
 if(.not.associated(item))return
 scan_item=>item
 do
  if((within_context.and.associated(context_stack%item,scan_item)).or.&
   &associated(casl_master,scan_item))exit
  if(len_metachar(label)==0)then
   if(.not.reduce)then
    call copy_metachar(scan_item%name,label)
   else
    call copy_metachar(scan_item%unique_name,label)
   endif
  else
   if(.not.reduce)then
    call prepend_char_to_metachar(label,&
     &trim(char_metachar(scan_item%name,len_metachar(scan_item%name))//':'))
   else
    call prepend_char_to_metachar(label,&
     &trim(char_metachar(scan_item%unique_name,&
     &len_metachar(scan_item%unique_name))//':'))
   endif
  endif
  scan_item=>scan_item%parent
 enddo
 END SUBROUTINE label_casl_item


 SUBROUTINE get_casl_item_D(label,value,ierr,respect_spelling)
!------------------------------------------------------------!
! Put value of CASL item labelled LABEL into VALUE (dp real) !
! and "project" the string value to real.                    !
!------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(out) :: ierr
 REAL(kind(1.d0)),INTENT(out) :: value
 LOGICAL,INTENT(in),OPTIONAL :: respect_spelling
 CHARACTER(*),INTENT(in) :: label
 LOGICAL keep_name
 CHARACTER(len(label)) crumb
 CHARACTER(80) tempr
 TYPE(casl_item),POINTER :: item
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 value=0.d0 ; ierr=-1
 call find_avltree_item(label,item)
 if(.not.associated(item))return
 if(.not.keep_name.and.len_trim(label)>0)then
  call extract_last_crumb(trim(label),crumb)
  call copy_char_to_metachar(trim(crumb),item%name)
 endif
 call read_from_char_D(char_metachar(item%value,len_metachar(item%value)),&
  &value,ierr)
 call flag_casl_read(item)
 write(tempr,*)value
 call copy_char_to_metachar(trim(adjustl(tempr)),item%value)
 END SUBROUTINE get_casl_item_D


 SUBROUTINE read_from_char_D(string,value,ierr)
!---------------------------------------!
! Read dp real from a character string. !
!---------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(out) :: ierr
 REAL(kind(1.d0)),INTENT(out) :: value
 CHARACTER(*),INTENT(in) :: string
 REAL(kind(1.d0)) readtemp
 value=0.d0 ; ierr=0
 read(string,*,iostat=ierr)readtemp
 if(ierr==0)value=readtemp
 END SUBROUTINE read_from_char_D


 SUBROUTINE get_casl_item_R(label,value,ierr,respect_spelling)
!------------------------------------------------------------!
! Put value of CASL item labelled LABEL into VALUE (sp real) !
! and "project" the string value to real.                    !
!------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(out) :: ierr
 REAL(kind(1.0)),INTENT(out) :: value
 LOGICAL,INTENT(in),OPTIONAL :: respect_spelling
 CHARACTER(*),INTENT(in) :: label
 LOGICAL keep_name
 CHARACTER(len(label)) crumb
 CHARACTER(80) tempr
 TYPE(casl_item),POINTER :: item
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 value=0. ; ierr=-1
 call find_avltree_item(label,item)
 if(.not.associated(item))return
 if(.not.keep_name.and.len_trim(label)>0)then
  call extract_last_crumb(trim(label),crumb)
  call copy_char_to_metachar(trim(crumb),item%name)
 endif
 call read_from_char_R(char_metachar(item%value,len_metachar(item%value)),&
  &value,ierr)
 call flag_casl_read(item)
 write(tempr,*)value
 call copy_char_to_metachar(trim(adjustl(tempr)),item%value)
 END SUBROUTINE get_casl_item_R


 SUBROUTINE read_from_char_R(string,value,ierr)
!---------------------------------------!
! Read sp real from a character string. !
!---------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(out) :: ierr
 REAL(kind(1.0)),INTENT(out) :: value
 CHARACTER(*),INTENT(in) :: string
 REAL(kind(1.0)) readtemp
 value=0.0 ; ierr=0
 read(string,*,iostat=ierr)readtemp
 if(ierr==0)value=readtemp
 END SUBROUTINE read_from_char_R


 SUBROUTINE get_casl_item_Z(label,value,ierr,respect_spelling)
!---------------------------------------------------------------!
! Put value of CASL item labelled LABEL into VALUE (dp complex) !
! and "project" the string value to complex.                    !
!---------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(out) :: ierr
 COMPLEX(kind(1.d0)),INTENT(out) :: value
 LOGICAL,INTENT(in),OPTIONAL :: respect_spelling
 CHARACTER(*),INTENT(in) :: label
 LOGICAL keep_name
 CHARACTER(len(label)) crumb
 CHARACTER(80) tempr
 TYPE(casl_item),POINTER :: item
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 value=cmplx(0.d0,0.d0,kind(1.d0)) ; ierr=-1
 call find_avltree_item(label,item)
 if(.not.associated(item))return
 if(.not.keep_name.and.len_trim(label)>0)then
  call extract_last_crumb(trim(label),crumb)
  call copy_char_to_metachar(trim(crumb),item%name)
 endif
 call read_from_char_Z(char_metachar(item%value,len_metachar(item%value)),&
  &value,ierr)
 call flag_casl_read(item)
 write(tempr,*)value
 call copy_char_to_metachar(trim(adjustl(tempr)),item%value)
 END SUBROUTINE get_casl_item_Z


 SUBROUTINE read_from_char_Z(string,value,ierr)
!------------------------------------------!
! Read dp complex from a character string. !
!------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(out) :: ierr
 COMPLEX(kind(1.d0)),INTENT(out) :: value
 CHARACTER(*),INTENT(in) :: string
 COMPLEX(kind(1.d0)) readtemp
 value=cmplx(0.d0,0.d0,kind(1.d0)) ; ierr=0
 read(string,*,iostat=ierr)readtemp
 if(ierr==0)value=readtemp
 END SUBROUTINE read_from_char_Z


 SUBROUTINE get_casl_item_X(label,value,ierr,respect_spelling)
!---------------------------------------------------------------!
! Put value of CASL item labelled LABEL into VALUE (sp complex) !
! and "project" the string value to complex.                    !
!---------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(out) :: ierr
 COMPLEX(kind(1.0)),INTENT(out) :: value
 LOGICAL,INTENT(in),OPTIONAL :: respect_spelling
 CHARACTER(*),INTENT(in) :: label
 LOGICAL keep_name
 CHARACTER(len(label)) crumb
 CHARACTER(80) tempr
 TYPE(casl_item),POINTER :: item
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 value=cmplx(0.,0.,kind(1.0)) ; ierr=-1
 call find_avltree_item(label,item)
 if(.not.associated(item))return
 if(.not.keep_name.and.len_trim(label)>0)then
  call extract_last_crumb(trim(label),crumb)
  call copy_char_to_metachar(trim(crumb),item%name)
 endif
 call read_from_char_X(char_metachar(item%value,len_metachar(item%value)),&
  &value,ierr)
 call flag_casl_read(item)
 write(tempr,*)value
 call copy_char_to_metachar(trim(adjustl(tempr)),item%value)
 END SUBROUTINE get_casl_item_X


 SUBROUTINE read_from_char_X(string,value,ierr)
!------------------------------------------!
! Read sp complex from a character string. !
!------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(out) :: ierr
 COMPLEX(kind(1.0)),INTENT(out) :: value
 CHARACTER(*),INTENT(in) :: string
 COMPLEX(kind(1.0)) readtemp
 value=cmplx(0.0,0.0,kind(1.0)) ; ierr=0
 read(string,*,iostat=ierr)readtemp
 if(ierr==0)value=readtemp
 END SUBROUTINE read_from_char_X


 SUBROUTINE get_casl_item_I(label,value,ierr,respect_spelling)
!------------------------------------------------------------!
! Put value of CASL item labelled LABEL into VALUE (integer) !
! and "project" the string value to integer.                 !
!------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(out) :: ierr,value
 LOGICAL,INTENT(in),OPTIONAL :: respect_spelling
 CHARACTER(*),INTENT(in) :: label
 LOGICAL keep_name
 CHARACTER(len(label)) crumb
 TYPE(casl_item),POINTER :: item
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 value=0 ; ierr=-1
 call find_avltree_item(label,item)
 if(.not.associated(item))return
 if(.not.keep_name.and.len_trim(label)>0)then
  call extract_last_crumb(trim(label),crumb)
  call copy_char_to_metachar(trim(crumb),item%name)
 endif
 call read_from_char_I(char_metachar(item%value,len_metachar(item%value)),&
  &value,ierr)
 if(ierr/=0)return
 call flag_casl_read(item)
 call copy_char_to_metachar(trim(i2s(value)),item%value)
 END SUBROUTINE get_casl_item_I


 SUBROUTINE read_from_char_I(string,value,ierr)
!---------------------------------------!
! Read integer from a character string. !
!---------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(out) :: ierr,value
 CHARACTER(*),INTENT(in) :: string
 INTEGER readtemp
 value=0 ; ierr=0
 read(string,*,iostat=ierr)readtemp
 if(ierr==0)value=readtemp
 END SUBROUTINE read_from_char_I


 SUBROUTINE get_casl_item_C(label,value,ierr,respect_spelling)
!---------------------------------------------------------------!
! Put value of CASL item labelled LABEL into VALUE (character). !
!---------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(out) :: ierr
 LOGICAL,INTENT(in),OPTIONAL :: respect_spelling
 CHARACTER(*),INTENT(in) :: label
 CHARACTER(*),INTENT(out) :: value
 LOGICAL keep_name
 CHARACTER(len(label)) crumb
 TYPE(casl_item),POINTER :: item
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 value="" ; ierr=-1
 call find_avltree_item(label,item)
 if(.not.associated(item))return
 if(.not.keep_name.and.len_trim(label)>0)then
  call extract_last_crumb(trim(label),crumb)
  call copy_char_to_metachar(trim(crumb),item%name)
 endif
 ierr=0
 call copy_metachar_to_char(item%value,value)
 call flag_casl_read(item)
 END SUBROUTINE get_casl_item_C


 SUBROUTINE get_casl_item_L(label,value,ierr,respect_spelling)
!------------------------------------------------------------!
! Put value of CASL item labelled LABEL into VALUE (logical) !
! and "project" the string value to logical.                 !
!------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(out) :: ierr
 LOGICAL,INTENT(in),OPTIONAL :: respect_spelling
 LOGICAL,INTENT(out) :: value
 CHARACTER(*),INTENT(in) :: label
 LOGICAL keep_name
 CHARACTER(len(label)) crumb
 TYPE(casl_item),POINTER :: item
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 value=.false. ; ierr=-1
 call find_avltree_item(label,item)
 if(.not.associated(item))return
 if(.not.keep_name.and.len_trim(label)>0)then
  call extract_last_crumb(trim(label),crumb)
  call copy_char_to_metachar(trim(crumb),item%name)
 endif
 select case(unique_casl_string(trim(char_metachar(item%value,&
  &len_metachar(item%value)))))
 case('.true.','true','t','yes','y','1') ; value=.true.
 case('.false.','false','f','no','n','0') ; value=.false.
 case default
  return
 end select
 ierr=0
 call flag_casl_read(item)
 if(value)then
  call copy_char_to_metachar('T',item%value)
 else
  call copy_char_to_metachar('F',item%value)
 endif
 END SUBROUTINE get_casl_item_L


 SUBROUTINE flag_casl_read(item)
!--------------------------------------------------------!
! Flag an item as read, along with all of its ancestors. !
!--------------------------------------------------------!
 IMPLICIT NONE
 TYPE(casl_item),POINTER :: item
 TYPE(casl_item),POINTER :: scan_item
 scan_item=>item
 do
  if(scan_item%been_read)exit
  scan_item%been_read=.true.
  if(.not.associated(scan_item%parent))exit
  scan_item=>scan_item%parent
 enddo
 END SUBROUTINE flag_casl_read


 SUBROUTINE force_casl_noninline(item)
!---------------------------------------------------------------!
! Make a CASL item non-inline, along with all of its ancestors. !
!---------------------------------------------------------------!
 IMPLICIT NONE
 TYPE(casl_item),POINTER :: item
 TYPE(casl_item),POINTER :: scan_item
 item%is_inline=.true.
 scan_item=>item
 do
  if(.not.scan_item%is_inline)exit
  scan_item%is_inline=.false.
  if(.not.associated(scan_item%parent))exit
  scan_item=>scan_item%parent
 enddo
 END SUBROUTINE force_casl_noninline


 SUBROUTINE query_casl_item(label,exists,is_block,nchildren,respect_spelling)
!------------------------------------!
! Returns properties of a CASL item. !
!------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(out),OPTIONAL :: nchildren
 LOGICAL,INTENT(in),OPTIONAL :: respect_spelling
 LOGICAL,INTENT(out),OPTIONAL :: exists,is_block
 CHARACTER(*),INTENT(in) :: label
 LOGICAL keep_name
 CHARACTER(len(label)) crumb
 TYPE(casl_item),POINTER :: item
 call find_avltree_item(label,item)
 if(.not.associated(item))then
  if(present(exists))exists=.false.
  if(present(is_block))is_block=.false.
  if(present(nchildren))nchildren=0
 else
  if(present(exists))exists=.true.
  if(present(is_block))is_block=item%is_block
! Use preferred name.
  keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
  if(.not.keep_name.and.len_trim(label)>0)then
   call extract_last_crumb(label,crumb)
   call copy_char_to_metachar(trim(crumb),item%name)
  endif
  if(present(nchildren))then
   nchildren=0
   if(item%is_block.and.associated(item%first_child))then
    item=>item%first_child
    do
     nchildren=nchildren+1
     if(.not.associated(item%next))exit
     item=>item%next
    enddo
   endif
  endif
 endif
 END SUBROUTINE query_casl_item


 SUBROUTINE first_unread_child(label,child_unique_name,ierr,flag_as_read)
!-----------------------------------------------------!
! Return the name of the first unread child of LABEL. !
!-----------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(out) :: ierr
 LOGICAL,INTENT(in),OPTIONAL :: flag_as_read
 CHARACTER(*),INTENT(in) :: label
 CHARACTER(CASL_KEYSIZE),INTENT(out) :: child_unique_name
 TYPE(casl_item),POINTER :: item
 child_unique_name='' ; ierr=-1
 call find_avltree_item(label,item)
 if(.not.associated(item))return
 ierr=-2
 if(.not.(item%is_block.and.associated(item%first_child)))return
 ierr=-3
 item=>item%first_child
 do
  if(.not.item%been_read)then
   ierr=0 ; call copy_metachar_to_char(item%unique_name,child_unique_name)
   if(present(flag_as_read))then
    if(flag_as_read)item%been_read=.true.
   endif
   return
  endif
  if(.not.associated(item%next))exit
  item=>item%next
 enddo
 END SUBROUTINE first_unread_child


 SUBROUTINE check_unread_casl(label,errmsg)
!------------------------------------------------------------------!
! Check if there are any unread CASL items under LABEL. All unread !
! items are printed as a warning, and if IS_ERROR is supplied and  !
! is T, the routine raises an error instead of returning.          !
! NB, LABEL must correspond to a block, else no check is done.     !
! NB2, empty sub-blocks will not be checked.                       !
!------------------------------------------------------------------!
! FIXME: this is mostly untested
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 CHARACTER(*),INTENT(inout) :: errmsg
 INTEGER count_unread
 TYPE(casl_item),POINTER :: item,scan_item

 count_unread=0
 errmsg=''

 call find_avltree_item(label,item)
 if(.not.associated(item))return ! technically there are no unread items
 scan_item=>item
 do
  if(scan_item%is_block.and.associated(scan_item%first_child))then
   scan_item=>scan_item%first_child
   cycle
  endif
  if(.not.scan_item%been_read)count_unread=count_unread+1
  if(associated(scan_item,item))exit
  do while(.not.associated(scan_item%next))
   scan_item=>scan_item%parent
   if(associated(scan_item,item))exit
  enddo
  if(associated(scan_item,item))exit
  scan_item=>scan_item%next
 enddo

 if(count_unread>0)errmsg='CASL item "'//label//'" contains '//&
  &trim(i2s(count_unread))//' unrecognized entries.'

 END SUBROUTINE check_unread_casl


 SUBROUTINE delete_casl_item(label)
!-------------------------!
! Remove CASL item LABEL. !
!-------------------------!
 CHARACTER(*),INTENT(in) :: label
 TYPE(casl_item),POINTER :: item
 call find_avltree_item(label,item)
 if(.not.associated(item))return
 call kill_casl_item(item)
 END SUBROUTINE delete_casl_item


 SUBROUTINE kill_casl_item(item)
!---------------------------------------------!
! Remove CASL item ITEM, including all of its !
! contents if it is a block.                  !
!---------------------------------------------!
 IMPLICIT NONE
 TYPE(casl_item),POINTER :: item
 TYPE(casl_item),POINTER :: parent,scan_item,scan_item_next
 INTEGER iimplicit,jimplicit,idum
 LOGICAL ldum
 CHARACTER(32) impl_name

! Kill contents.
 if(item%is_block)then
  if(associated(item%first_child))then
   scan_item=>item%first_child
   do
    if(scan_item%is_block.and.associated(scan_item%first_child))then
     scan_item=>scan_item%first_child
     nullify(scan_item%parent%first_child)
     cycle
    endif
    call delete_avltree_item(scan_item)
    call unset_metachar(scan_item%value)
    call unset_metachar(scan_item%name)
    call unset_metachar(scan_item%unique_name)
    call unset_metachar(scan_item%full_unique_name)
    if(associated(scan_item%next))then
     scan_item_next=>scan_item%next
    else
     scan_item_next=>scan_item%parent
    endif
    deallocate(scan_item)
    scan_item=>scan_item_next
    if(associated(scan_item,item))exit
   enddo
  endif
 else ! not a block
  call unset_metachar(item%value)
 endif ! block or not

! Bridge over item in list.
 if(associated(item%parent))then
  parent=>item%parent
  if(associated(parent%first_child,item))then ! first child
   if(.not.associated(item%next))then ! only child
    nullify(parent%first_child,parent%last_child)
   else ! not only child
    parent%first_child=>item%next
    nullify(item%next%prev)
   endif
  elseif(associated(parent%last_child,item))then ! last, not only child
   parent%last_child=>item%prev
   nullify(item%prev%next)
  else ! middle child
   item%prev%next=>item%next
   item%next%prev=>item%prev
  endif
 else
  nullify(parent)
 endif ! item has parent

! Remove item from AVL tree.
 call delete_avltree_item(item)

! Adjust implicit-item count in parent and implicit names in siblings.
 if(associated(parent))then
  call copy_metachar_to_char(item%name,impl_name)
  if(impl_name(1:2)=='%u')then
   parent%nimplicit=parent%nimplicit-1
   if(associated(parent%first_child))then
    call read_from_char_I(trim(impl_name(3:)),iimplicit,idum)
    if(iimplicit<=parent%nimplicit)then
     scan_item=>parent%first_child
     do
      call copy_metachar_to_char(scan_item%name,impl_name)
      if(impl_name(1:2)=='%u')then
       call read_from_char_I(trim(impl_name(3:)),jimplicit,idum)
       if(jimplicit>iimplicit)then
! Remove from AVL tree.
        call delete_avltree_item(scan_item)
! Rename item.
        call copy_char_to_metachar('%u'//trim(i2s(jimplicit-1)),scan_item%name)
        call copy_char_to_metachar('%u'//trim(i2s(jimplicit-1)),&
         &scan_item%unique_name)
        call copy_metachar(scan_item%unique_name,scan_item%full_unique_name)
        call prepend_char_to_metachar(scan_item%full_unique_name,':')
        call prepend_metachar_to_metachar(scan_item%full_unique_name,&
         &scan_item%parent%full_unique_name)
! Reinsert in AVL tree.
        call insert_avltree_item(scan_item,ldum)
       endif
       if(jimplicit==parent%nimplicit)exit
      endif
      if(.not.associated(scan_item%next))exit
      scan_item=>scan_item%next
     enddo
    endif
   endif
  endif
 endif

! Delete metachars.
 call unset_metachar(item%name)
 call unset_metachar(item%unique_name)
 call unset_metachar(item%full_unique_name)
 call unset_metachar(item%value)

! Purge item.
 deallocate(item)
 nullify(item)

 END SUBROUTINE kill_casl_item


 SUBROUTINE create_casl_item(label,item,errmsg,value,only_if_unset,as_block,&
  &respect_spelling)
!-------------------------------------------------------!
! Create a CASL item if it doesn't exist, including all !
! ancestors if necessary, and return a pointer ITEM.    !
! If AS_BLOCK=T, create a block, else create a scalar   !
! item and set ITEM%VALUE to VALUE, unless ITEM already !
! exists and ONLY_IF_UNSET=T.                           !
!-------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 TYPE(casl_item),POINTER :: item
 CHARACTER(512),INTENT(inout) :: errmsg
 CHARACTER(*),INTENT(in),OPTIONAL :: value
 LOGICAL,INTENT(in),OPTIONAL :: only_if_unset,as_block,respect_spelling
 INTEGER iseq,ierr,icrumb,ncrumb_mlabel,ncrumb_label
 LOGICAL item_present,create_block,if_unset,keep_name
 CHARACTER(len(label)) crumb
 TYPE(casl_item),POINTER :: parent
 TYPE(metachar) mlabel,mlabel2,ulabel,ulabel2

! Initialize.
 call initialize_casl
 errmsg=''
 create_block=.false. ; if(present(as_block))create_block=as_block
 if_unset=.false. ; if(present(only_if_unset))if_unset=only_if_unset
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling

! Get metachar version of full unique label.
 call full_path_of_label(label,mlabel)

! Just return if it exists.
 call find_avltree_item(':'//char_metachar(mlabel,len_metachar(mlabel)),item)
 if(associated(item))then
  if(create_block.neqv.item%is_block)then
   errmsg='Block status mismatch for '//label//'.  Check if you have defined &
    &this item as a scalar/block in the relevant .casl file when it should &
    &have been a block/scalar (e.g., you may have tried to define an empty &
    &block but forgot the square brackets "[ ]" after the colon, or you may &
    &have accidentally indented the line after an empty scalar).'
   call unset_metachar(mlabel)
   return
  endif
! Replace name if not told otherwise.
  if(.not.keep_name.and.len_trim(label)>0)then
   call extract_last_crumb(trim(label),crumb)
   call copy_char_to_metachar(trim(crumb),item%name)
  endif
! Replace value if not told otherwise.
  if(.not.(if_unset.or.create_block))call copy_char_to_metachar(value,&
   &item%value)
! Flag the item as read and return.
  call flag_casl_read(item)
  call unset_metachar(mlabel)
  return
 endif

! Find closest ancestor.
 call copy_metachar(mlabel,mlabel2)
 do while(len_metachar(mlabel2)>0)
  call remove_last_crumb_metachar(mlabel2)
  call find_avltree_item(':'//char_metachar(mlabel2,len_metachar(mlabel2)),&
   &parent)
  if(associated(parent))exit
 enddo
 if(len_metachar(mlabel2)>0)call cut_metachar(mlabel,len_metachar(mlabel2)+2)

! Construct un-reduced partial label.
 ncrumb_mlabel=count_crumbs(char_metachar(mlabel,len_metachar(mlabel)))
 ncrumb_label=count_crumbs(trim(label))
 call copy_char_to_metachar(trim(label),ulabel)
 do icrumb=ncrumb_mlabel,ncrumb_label-1
  call split_first_crumb_metachar(ulabel,ulabel2)
 enddo ! icrumb

! Create missing ancestors.
 call split_first_crumb_metachar(mlabel,mlabel2)
 call split_first_crumb_metachar(ulabel,ulabel2)
 if(len_metachar(mlabel)>0)then
  do while(len_metachar(mlabel)>0)
   allocate(item)
   call copy_metachar(ulabel2,item%name)
   if(associated(parent,casl_master))then
! This is a filename, so do not reduce name.
    call copy_metachar(ulabel2,item%unique_name)
    call copy_metachar(ulabel2,item%full_unique_name)
   else
    call copy_metachar(mlabel2,item%unique_name)
    call copy_metachar(item%unique_name,item%full_unique_name)
    call prepend_char_to_metachar(item%full_unique_name,':')
    call prepend_metachar_to_metachar(item%full_unique_name,&
     &parent%full_unique_name)
   endif
   item%is_block=.true.
   item%ilevel=parent%ilevel+1
   item%is_inline=parent%is_inline
   item%parent=>parent
   nullify(item%next)
   if(associated(parent%first_child))then
    parent%last_child%next=>item
    item%prev=>parent%last_child
    parent%last_child=>item
   else ! first child
    parent%first_child=>item
    parent%last_child=>item
    nullify(item%prev)
   endif
   call insert_avltree_item(item,item_present)
   if(item_present)then
    errmsg='Failed to insert new item "'//char_metachar(item%full_unique_name,&
     &len_metachar(item%full_unique_name))//'" in AVL tree since it appears &
     &to be already present (ancestor).'
    call unset_metachar(mlabel)
    call unset_metachar(mlabel2)
    call unset_metachar(ulabel)
    call unset_metachar(ulabel2)
    return
   endif
   parent=>item
   nullify(item)
   call split_first_crumb_metachar(mlabel,mlabel2)
   call split_first_crumb_metachar(ulabel,ulabel2)
  enddo ! len_metachar(mlabel)>0
 endif ! len_metachar(mlabel)>0

! Clean metachars.
 call unset_metachar(mlabel)
 call unset_metachar(ulabel)

! Create the item and set its name.
 allocate(item)
 call copy_metachar(ulabel2,item%name)
 if(associated(parent,casl_master))then
! This is a filename, so do not reduce name.
  call copy_metachar(ulabel2,item%unique_name)
 else
  call copy_metachar(mlabel2,item%unique_name)
 endif
 item%is_block=create_block
 item%ilevel=parent%ilevel+1
 item%been_read=.true.
 item%is_inline=parent%is_inline
 item%parent=>parent
 if(.not.create_block)then
! Regular item: create value, check name is valid, deal with implicit items.
  call copy_metachar_to_char(mlabel2,crumb)
  if(crumb(1:1)=='%')then
   if(crumb(2:2)=='u')then
    item%is_implicit=.true.
    parent%nimplicit=parent%nimplicit+1
    if(trim(crumb)=='%u')then
     call append_char_to_metachar(mlabel2,trim(i2s(parent%nimplicit)))
     call copy_metachar(mlabel2,item%name)
     call copy_metachar(mlabel2,item%unique_name)
    else
     call read_from_char_I(crumb(3:),iseq,ierr)
     if(ierr/=0)then
      errmsg='Problem reading implicit name index in label.'
      call unset_metachar(mlabel2)
      call unset_metachar(ulabel2)
      return
     endif
     if(parent%nimplicit/=iseq)then
      errmsg='Problem counting implicitly-named items.'
      call unset_metachar(mlabel2)
      call unset_metachar(ulabel2)
      return
     endif
    endif
   else
    errmsg='Cannot create non-implicit item whose name starts with "%".'
    call unset_metachar(mlabel2)
    call unset_metachar(ulabel2)
    return
   endif
  endif
  if(present(value))then
   call copy_char_to_metachar(trim(value),item%value)
  else
   call unset_metachar(item%value)
  endif
 else
! Block item: check name is valid.
  if(mlabel2%chars(1)=='%')then
   errmsg='Cannot create block item whose name starts with "%".'
   call unset_metachar(mlabel2)
   call unset_metachar(ulabel2)
   return
  endif
 endif

! Clean metachars.
 call unset_metachar(mlabel2)
 call unset_metachar(ulabel2)

! Full name
 call copy_metachar(item%unique_name,item%full_unique_name)
 if(.not.associated(parent,casl_master))then
  call prepend_char_to_metachar(item%full_unique_name,':')
  call prepend_metachar_to_metachar(item%full_unique_name,&
   &parent%full_unique_name)
 endif

! Insert in tree.
 if(associated(parent%first_child))then
  parent%last_child%next=>item
  item%prev=>parent%last_child
  parent%last_child=>item
 else ! first child
  parent%first_child=>item
  parent%last_child=>item
  nullify(item%prev)
 endif

! Insert in AVL tree.
 call insert_avltree_item(item,item_present)
 if(item_present)then
  errmsg='Failed to insert new item "'//char_metachar(item%full_unique_name,&
   &len_metachar(item%full_unique_name))//'" in AVL tree since it appears &
   &to be already present.'
  return
 endif

 END SUBROUTINE create_casl_item


 SUBROUTINE full_path_of_label(label,path)
!---------------------------------------------------------------!
! Compute full path to LABEL (which is potentially relative to  !
! context, unless LABEL starts with ':') and return it in PATH. !
!---------------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 TYPE(metachar),INTENT(inout) :: path
 TYPE(metachar) tstring,tstring2

! Concatenate path-to-context and LABEL into PATH.
 call copy_char_to_metachar(label,path)
 if(len_metachar(path)>0)then
  if(path%chars(1)/=':')then
   call prepend_char_to_metachar(path,':')
   call prepend_metachar_to_metachar(path,context_stack%item%full_unique_name)
  endif
 endif

! Remove any ':' at beginning
 if(len_metachar(path)>0)then
  do while(path%chars(1)==':')
   call cut_metachar(path,2)
   if(len_metachar(path)==0)exit
  enddo
 endif

 if(len_metachar(path)>0)then

! Extract first crumb.
  call split_first_crumb_metachar(path,tstring)

! Reduce everything but first crumb and reconstruct path.
  if(len_metachar(tstring)>0.and.len_metachar(path)>0)then
   call copy_char_to_metachar(&
    &trim(unique_casl_string(char_metachar(path,len_metachar(path)))),tstring2)
   call prepend_char_to_metachar(tstring2,':')
   call prepend_metachar_to_metachar(tstring2,tstring)
   call copy_metachar(tstring2,path)
   call unset_metachar(tstring2)
  elseif(len_metachar(tstring)>0)then
   call copy_metachar(tstring,path)
  endif
  call unset_metachar(tstring)

 endif

 END SUBROUTINE full_path_of_label


 SUBROUTINE set_casl_item_D(label,value,errmsg,only_if_unset,respect_spelling)
!------------------------------------------------------!
! Set the value of CASL item LABEL to VALUE (real dp), !
! creating item LABEL and its parents if necessary.    !
!------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 REAL(kind(1.d0)),INTENT(in) :: value
 CHARACTER(512),INTENT(inout) :: errmsg
 LOGICAL,INTENT(in),OPTIONAL :: only_if_unset,respect_spelling
 LOGICAL if_unset,keep_name
 CHARACTER(80) tmpr
 TYPE(casl_item),POINTER :: item
 errmsg=''
 if_unset=.false. ; if(present(only_if_unset))if_unset=only_if_unset
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 write(tmpr,*)value
 call create_casl_item(label,item,errmsg,value=trim(adjustl(tmpr)),&
  &only_if_unset=if_unset,respect_spelling=keep_name)
 END SUBROUTINE set_casl_item_D


 SUBROUTINE set_casl_item_R(label,value,errmsg,only_if_unset,respect_spelling)
!------------------------------------------------------!
! Set the value of CASL item LABEL to VALUE (real sp), !
! creating item LABEL and its parents if necessary.    !
!------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 REAL(kind(1.0)),INTENT(in) :: value
 CHARACTER(512),INTENT(inout) :: errmsg
 LOGICAL,INTENT(in),OPTIONAL :: only_if_unset,respect_spelling
 LOGICAL if_unset,keep_name
 CHARACTER(80) tmpr
 TYPE(casl_item),POINTER :: item
 errmsg=''
 if_unset=.false. ; if(present(only_if_unset))if_unset=only_if_unset
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 write(tmpr,*)value
 call create_casl_item(label,item,errmsg,value=trim(adjustl(tmpr)),&
  &only_if_unset=if_unset,respect_spelling=keep_name)
 END SUBROUTINE set_casl_item_R


 SUBROUTINE set_casl_item_Z(label,value,errmsg,only_if_unset,respect_spelling)
!---------------------------------------------------------!
! Set the value of CASL item LABEL to VALUE (complex dp), !
! creating item LABEL and its parents if necessary.       !
!---------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 COMPLEX(kind(1.d0)),INTENT(in) :: value
 CHARACTER(512),INTENT(inout) :: errmsg
 LOGICAL,INTENT(in),OPTIONAL :: only_if_unset,respect_spelling
 LOGICAL if_unset,keep_name
 CHARACTER(80) tmpr
 TYPE(casl_item),POINTER :: item
 errmsg=''
 if_unset=.false. ; if(present(only_if_unset))if_unset=only_if_unset
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 write(tmpr,*)value
 call create_casl_item(label,item,errmsg,value=trim(adjustl(tmpr)),&
  &only_if_unset=if_unset,respect_spelling=keep_name)
 END SUBROUTINE set_casl_item_Z


 SUBROUTINE set_casl_item_X(label,value,errmsg,only_if_unset,respect_spelling)
!---------------------------------------------------------!
! Set the value of CASL item LABEL to VALUE (complex sp), !
! creating item LABEL and its parents if necessary.       !
!---------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 COMPLEX(kind(1.0)),INTENT(in) :: value
 CHARACTER(512),INTENT(inout) :: errmsg
 LOGICAL,INTENT(in),OPTIONAL :: only_if_unset,respect_spelling
 LOGICAL if_unset,keep_name
 CHARACTER(80) tmpr
 TYPE(casl_item),POINTER :: item
 errmsg=''
 if_unset=.false. ; if(present(only_if_unset))if_unset=only_if_unset
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 write(tmpr,*)value
 call create_casl_item(label,item,errmsg,value=trim(adjustl(tmpr)),&
  &only_if_unset=if_unset,respect_spelling=keep_name)
 END SUBROUTINE set_casl_item_X


 SUBROUTINE set_casl_item_I(label,value,errmsg,only_if_unset,respect_spelling)
!------------------------------------------------------!
! Set the value of CASL item LABEL to VALUE (integer), !
! creating item LABEL and its parents if necessary.    !
!------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 INTEGER,INTENT(in) :: value
 CHARACTER(512),INTENT(inout) :: errmsg
 LOGICAL,INTENT(in),OPTIONAL :: only_if_unset,respect_spelling
 LOGICAL if_unset,keep_name
 CHARACTER(80) tmpr
 TYPE(casl_item),POINTER :: item
 errmsg=''
 if_unset=.false. ; if(present(only_if_unset))if_unset=only_if_unset
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 write(tmpr,*)value
 call create_casl_item(label,item,errmsg,value=trim(adjustl(tmpr)),&
  &only_if_unset=if_unset,respect_spelling=keep_name)
 END SUBROUTINE set_casl_item_I


 SUBROUTINE set_casl_item_C(label,value,errmsg,only_if_unset,respect_spelling)
!--------------------------------------------------------!
! Set the value of CASL item LABEL to VALUE (character), !
! creating item LABEL and its parents if necessary.      !
!--------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label,value
 CHARACTER(512),INTENT(inout) :: errmsg
 LOGICAL,INTENT(in),OPTIONAL :: only_if_unset,respect_spelling
 LOGICAL if_unset,keep_name
 TYPE(casl_item),POINTER :: item
 errmsg=''
 if_unset=.false. ; if(present(only_if_unset))if_unset=only_if_unset
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 call create_casl_item(label,item,errmsg,value=value,only_if_unset=if_unset,&
  &respect_spelling=keep_name)
 END SUBROUTINE set_casl_item_C


 SUBROUTINE set_casl_item_L(label,value,errmsg,only_if_unset,respect_spelling)
!------------------------------------------------------!
! Set the value of CASL item LABEL to VALUE (logical), !
! creating item LABEL and its parents if necessary.    !
!------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 LOGICAL,INTENT(in) :: value
 CHARACTER(512),INTENT(inout) :: errmsg
 LOGICAL,INTENT(in),OPTIONAL :: only_if_unset,respect_spelling
 LOGICAL if_unset,keep_name
 CHARACTER(1) char1
 TYPE(casl_item),POINTER :: item
 errmsg=''
 if_unset=.false. ; if(present(only_if_unset))if_unset=only_if_unset
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 if(value)then
  char1='T'
 else
  char1='F'
 endif
 call create_casl_item(label,item,errmsg,value=char1,only_if_unset=if_unset,&
  &respect_spelling=keep_name)
 END SUBROUTINE set_casl_item_L


 SUBROUTINE set_casl_block(label,errmsg,prefer_inline,force_noninline,&
  &respect_spelling,push)
!----------------------------------------------------!
! Create CASL item LABEL of type block if it doesn't !
! exist, also creating its parents if necessary.     !
! Optional flags have the following effects:         !
! - PREFER_INLINE: if .true., flag the contents of   !
!   this block to be written in inline mode.  A      !
!   value of .false. has no effect if any ancestor   !
!   of this block is in inline mode.                 !
! - FORCE_NONINLINE: if .true. flag this block and   !
!   its ancestors to be written in normal            !
!   (non-inline) mode.                               !
! - RESPECT_SPELLING: if .true. and the block        !
!   already exists, rewrite the block's unreduced    !
!   name with the spelling given in LABEL.           !
! - PUSH: if .true., push this block onto the CASL   !
!   context stack.                                   !
!----------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 CHARACTER(512),INTENT(inout) :: errmsg
 LOGICAL,INTENT(in),OPTIONAL :: prefer_inline,force_noninline,&
  &respect_spelling,push
 LOGICAL keep_name
 TYPE(casl_item),POINTER :: item
 errmsg=''
 keep_name=.false. ; if(present(respect_spelling))keep_name=respect_spelling
 call create_casl_item(label,item,errmsg,as_block=.true.,&
  &respect_spelling=keep_name)
 if(present(prefer_inline))then
  item%is_inline=prefer_inline
  if(associated(item%parent))then ! can't allow normal mode inside inline mode
   if(item%parent%is_inline)item%is_inline=.true.
  endif
 endif
 if(present(force_noninline))then
  if(force_noninline)call force_casl_noninline(item)
 endif
 if(present(push))then
  if(push)then
   allocate(context_stack%next)
   context_stack%next%prev=>context_stack
   context_stack%next%item=>item
   context_stack=>context_stack%next
  endif
 endif
 END SUBROUTINE set_casl_block


! Tools to handle AVL self-balancing binary search trees.


 SUBROUTINE insert_avltree_item(item,already_present)
!----------------------------------------------------!
! Insert ITEM into AVL tree, and ALREADY_PRESENT = T !
! if another item with the same unique_full_label is !
! already there.  ITEM%FULL_UNIQUE_NAME must be set  !
! before calling this routine, other components of   !
! ITEM are not used here.  ITEM%AVL_* are set by     !
! this routine, and the tree is rebalanced.          !
!----------------------------------------------------!
 IMPLICIT NONE
 TYPE(casl_item),POINTER :: item
 LOGICAL,INTENT(out) :: already_present
 TYPE(casl_item),POINTER :: scan_item
 already_present=.false.
! FIXME - could call find_avltree_item, which in turn would allow for
! inserts to be cached too.
 scan_item=>casl_avltree_root
 do
  select case(compare_metachar(item%full_unique_name,&
   &scan_item%full_unique_name))
  case('=')
   already_present=.true.
   return
  case('>')
   if(associated(scan_item%avl_right))then
    scan_item=>scan_item%avl_right
   else
    scan_item%avl_right=>item
    item%avl_parent=>scan_item
    item%avl_depth=1
    call rebalance_avltree_up_from_item(item)
    return
   endif
  case('<')
   if(associated(scan_item%avl_left))then
    scan_item=>scan_item%avl_left
   else
    scan_item%avl_left=>item
    item%avl_parent=>scan_item
    item%avl_depth=1
    call rebalance_avltree_up_from_item(item)
    return
   endif
  end select
 enddo
 END SUBROUTINE insert_avltree_item


 SUBROUTINE rebalance_avltree_up_from_item(item)
!--------------------------------------------------------!
! Ensure the AVL tree is balanced at ITEM and all of its !
! ancestors.  This routine should be called following    !
! every insertion (with ITEM => inserted item) or        !
! deletion (with ITEM => lowest-level item whose         !
! children have been re-attached).                       !
!--------------------------------------------------------!
 IMPLICIT NONE
 TYPE(casl_item),POINTER :: item
 TYPE(casl_item),POINTER :: scan_item
 scan_item=>item
 do
  call update_avltree_depth(scan_item)
  call balance_avltree_item(scan_item)
  if(.not.associated(scan_item%avl_parent))then
   casl_avltree_root=>scan_item
   exit
  endif
  scan_item=>scan_item%avl_parent
 enddo
 END SUBROUTINE rebalance_avltree_up_from_item


 SUBROUTINE balance_avltree_item(item)
!--------------------------------------------------------!
! Carry out an individual balance operation so that the  !
! child trees of ITEM have a depth difference of at most !
! +1/-1.  On output, ITEM points at the item that has    !
! replaced the original ITEM (if any), so that recursive !
! balance operation can always continue at               !
! ITEM%AVL_PARENT after calling this routine.            !
!--------------------------------------------------------!
 IMPLICIT NONE
 TYPE(casl_item),POINTER :: item
 TYPE(casl_item),POINTER :: child,grandchild
 INTEGER imbal1,imbal2

 imbal1=item%avl_imbalance

 if(imbal1<-1)then

  child=>item%avl_left
  imbal2=child%avl_imbalance

  if(imbal2<0)then
! Left-left imbalance
   if(associated(child%avl_right))then
    item%avl_left=>child%avl_right
    item%avl_left%avl_parent=>item
   else
    nullify(item%avl_left)
   endif
   child%avl_right=>item
   if(associated(item%avl_parent))then
    child%avl_parent=>item%avl_parent
    if(associated(item%avl_parent%avl_left,item))then
     item%avl_parent%avl_left=>child
    else
     item%avl_parent%avl_right=>child
    endif
   else
    nullify(child%avl_parent)
   endif
   item%avl_parent=>child
   call update_avltree_depth(item)
   call update_avltree_depth(child)
   item=>child

  elseif(imbal2>0)then
! Left-right imbalance
   grandchild=>child%avl_right
   if(associated(grandchild%avl_left))then
    child%avl_right=>grandchild%avl_left
    child%avl_right%avl_parent=>child
   else
    nullify(child%avl_right)
   endif
   if(associated(grandchild%avl_right))then
    item%avl_left=>grandchild%avl_right
    item%avl_left%avl_parent=>item
   else
    nullify(item%avl_left)
   endif
   grandchild%avl_left=>child
   grandchild%avl_right=>item
   if(associated(item%avl_parent))then
    grandchild%avl_parent=>item%avl_parent
   else
    nullify(grandchild%avl_parent)
   endif
   if(associated(item%avl_parent))then
    grandchild%avl_parent=>item%avl_parent
    if(associated(item%avl_parent%avl_left,item))then
     grandchild%avl_parent%avl_left=>grandchild
    else
     grandchild%avl_parent%avl_right=>grandchild
    endif
   else
    nullify(child%avl_parent)
   endif
   grandchild%avl_left=>child
   child%avl_parent=>grandchild
   grandchild%avl_right=>item
   item%avl_parent=>grandchild
   call update_avltree_depth(child)
   call update_avltree_depth(item)
   call update_avltree_depth(grandchild)
   item=>grandchild
  endif

 elseif(imbal1>1)then

  child=>item%avl_right
  imbal2=child%avl_imbalance

  if(imbal2<0)then
! Right-left imbalance
   grandchild=>child%avl_left
   if(associated(grandchild%avl_right))then
    child%avl_left=>grandchild%avl_right
    child%avl_left%avl_parent=>child
   else
    nullify(child%avl_left)
   endif
   if(associated(grandchild%avl_left))then
    item%avl_right=>grandchild%avl_left
    item%avl_right%avl_parent=>item
   else
    nullify(item%avl_right)
   endif
   grandchild%avl_right=>child
   grandchild%avl_left=>item
   if(associated(item%avl_parent))then
    grandchild%avl_parent=>item%avl_parent
   else
    nullify(grandchild%avl_parent)
   endif
   if(associated(item%avl_parent))then
    grandchild%avl_parent=>item%avl_parent
    if(associated(item%avl_parent%avl_right,item))then
     grandchild%avl_parent%avl_right=>grandchild
    else
     grandchild%avl_parent%avl_left=>grandchild
    endif
   else
    nullify(child%avl_parent)
   endif
   grandchild%avl_right=>child
   child%avl_parent=>grandchild
   grandchild%avl_left=>item
   item%avl_parent=>grandchild
   call update_avltree_depth(child)
   call update_avltree_depth(item)
   call update_avltree_depth(grandchild)
   item=>grandchild

  elseif(imbal2>0)then
! Right-right imbalance
   if(associated(child%avl_left))then
    item%avl_right=>child%avl_left
    item%avl_right%avl_parent=>item
   else
    nullify(item%avl_right)
   endif
   child%avl_left=>item
   if(associated(item%avl_parent))then
    child%avl_parent=>item%avl_parent
    if(associated(item%avl_parent%avl_right,item))then
     item%avl_parent%avl_right=>child
    else
     item%avl_parent%avl_left=>child
    endif
   else
    nullify(child%avl_parent)
   endif
   item%avl_parent=>child
   call update_avltree_depth(item)
   call update_avltree_depth(child)
   item=>child
  endif

 endif

 END SUBROUTINE balance_avltree_item


 SUBROUTINE update_avltree_depth(item)
!--------------------------------------------------------!
! Update ITEM%AVL_*DEPTH and ITEM%AVL_IMBALANCE using    !
! information from ITEM's children, which are assumed to !
! be up-to-date.  This routine should be called for all  !
! ancestors of ITEM afterwards.                          !
!--------------------------------------------------------!
 IMPLICIT NONE
 TYPE(casl_item),POINTER :: item
 INTEGER left_depth,right_depth
 if(associated(item%avl_left))then
  left_depth=item%avl_left%avl_depth
 else
  left_depth=0
 endif
 if(associated(item%avl_right))then
  right_depth=item%avl_right%avl_depth
 else
  right_depth=0
 endif
 item%avl_depth=max(left_depth,right_depth)+1
 item%avl_imbalance=right_depth-left_depth
 END SUBROUTINE update_avltree_depth


 SUBROUTINE delete_avltree_item(item)
!------------------------------------------------------!
! Delete ITEM from the AVL tree.  This routine can be  !
! called before or after deleting the CASL information !
! from ITEM (but it obviously must be called before    !
! destroying ITEM).  The AVL tree is rebalanced after  !
! deletion.                                            !
!------------------------------------------------------!
 IMPLICIT NONE
 TYPE(casl_item),POINTER :: item
 TYPE(casl_item),POINTER :: new_top,new_bottom

 nullify(new_top,new_bottom)

 if(.not.associated(item%avl_left).and..not.associated(item%avl_right))then
! Item has no children, so it can be safely deleted.
  if(associated(item%avl_parent))then
   if(associated(item%avl_parent%avl_left,item))then
    nullify(item%avl_parent%avl_left)
   else
    nullify(item%avl_parent%avl_right)
   endif
   new_bottom=>item%avl_parent
  else ! no parent either, so we don't have an AVL tree
   nullify(casl_avltree_root)
  endif

 elseif(.not.associated(item%avl_left))then
! Item has only a right child, which replaces it.
  if(associated(item%avl_parent))then
   if(associated(item%avl_parent%avl_left,item))then
    item%avl_parent%avl_left=>item%avl_right
   else
    item%avl_parent%avl_right=>item%avl_right
   endif
   item%avl_right%avl_parent=>item%avl_parent
   new_bottom=>item%avl_right
  else
   casl_avltree_root=>item%avl_right
   nullify(item%avl_right%avl_parent)
  endif

 elseif(.not.associated(item%avl_right))then
! Item has only a left child, which replaces it.
  if(associated(item%avl_parent))then
   if(associated(item%avl_parent%avl_left,item))then
    item%avl_parent%avl_left=>item%avl_left
   else
    item%avl_parent%avl_right=>item%avl_left
   endif
   item%avl_left%avl_parent=>item%avl_parent
   new_bottom=>item%avl_left
  else
   casl_avltree_root=>item%avl_left
   nullify(item%avl_left%avl_parent)
  endif

 else
! Item has two children, so search for replacement following them.
  if(item%avl_imbalance>0)then
! Item has a right imbalance, so search on the right bearing left.
   new_top=>item%avl_right
   do while(associated(new_top%avl_left))
    new_top=>new_top%avl_left
   enddo
   new_bottom=>new_top%avl_parent
   if(associated(new_bottom,item))then
    nullify(new_bottom)
   else
    if(associated(new_top%avl_right))then
     new_bottom%avl_left=>new_top%avl_right
     new_bottom%avl_left%avl_parent=>new_bottom
    else
     nullify(new_bottom%avl_left)
    endif
   endif
   if(associated(item%avl_left))then
    new_top%avl_left=>item%avl_left
    new_top%avl_left%avl_parent=>new_top
   else
    nullify(new_top%avl_left)
   endif
   if(associated(new_bottom))then ! else new_top would become its own child!
    new_top%avl_right=>item%avl_right
    new_top%avl_right%avl_parent=>new_top
   endif
   if(associated(item%avl_parent))then
    new_top%avl_parent=>item%avl_parent
    if(associated(item%avl_parent%avl_left,item))then
     new_top%avl_parent%avl_left=>new_top
    else
     new_top%avl_parent%avl_right=>new_top
    endif
   else
    nullify(new_top%avl_parent)
    casl_avltree_root=>new_top
   endif
   if(.not.associated(new_bottom))new_bottom=>new_top

  else
! Item has a left imbalance (or none), so search on the left bearing right.
   new_top=>item%avl_left
   do while(associated(new_top%avl_right))
    new_top=>new_top%avl_right
   enddo
   new_bottom=>new_top%avl_parent
   if(associated(new_bottom,item))then
    nullify(new_bottom)
   else
    if(associated(new_top%avl_left))then
     new_bottom%avl_right=>new_top%avl_left
     new_bottom%avl_right%avl_parent=>new_bottom
    else
     nullify(new_bottom%avl_right)
    endif
   endif
   if(associated(item%avl_right))then
    new_top%avl_right=>item%avl_right
    new_top%avl_right%avl_parent=>new_top
   else
    nullify(new_top%avl_right)
   endif
   if(associated(new_bottom))then ! else new_top would become its own child!
    new_top%avl_left=>item%avl_left
    new_top%avl_left%avl_parent=>new_top
   endif
   if(associated(item%avl_parent))then
    new_top%avl_parent=>item%avl_parent
    if(associated(item%avl_parent%avl_right,item))then
     new_top%avl_parent%avl_right=>new_top
    else
     new_top%avl_parent%avl_left=>new_top
    endif
   else
    nullify(new_top%avl_parent)
    casl_avltree_root=>new_top
   endif
   if(.not.associated(new_bottom))new_bottom=>new_top

  endif

 endif

! Rebalance tree.
 if(associated(new_bottom))call rebalance_avltree_up_from_item(new_bottom)

! Fully disconnect item.
 nullify(item%avl_parent,item%avl_left,item%avl_right)
 item%avl_depth=0
 item%avl_imbalance=0

 END SUBROUTINE delete_avltree_item


 SUBROUTINE find_avltree_item(label,item)
!--------------------------------------------------------------!
! Locate item corresponding to LABEL in AVL tree and return it !
! via pointer ITEM, or point at NULL() if not found.           !
!--------------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 TYPE(casl_item),POINTER :: item
 TYPE(metachar) mlabel

! Prepare metachar version of full unique label.
 call full_path_of_label(label,mlabel)
 if(len_metachar(mlabel)==0)then
  item=>casl_master ; return
 endif

! Initialize search at tree root.
 item=>casl_avltree_root

! Scan tree to find item.
 do
  select case(compare_metachar(mlabel,item%full_unique_name))
  case('=')
   call unset_metachar(mlabel)
   return
  case('>')
   if(.not.associated(item%avl_right))then
    nullify(item)
    call unset_metachar(mlabel)
    return
   endif
   item=>item%avl_right
  case('<')
   if(.not.associated(item%avl_left))then
    nullify(item)
    call unset_metachar(mlabel)
    return
   endif
   item=>item%avl_left
  end select
 enddo

 END SUBROUTINE find_avltree_item


! Tools replacing the functionality previously provided by
! the problematic iso_varying_string module.


 SUBROUTINE unset_metachar(string)
!---------------!
! Unset STRING. '
!---------------!
 IMPLICIT NONE
 TYPE(metachar),INTENT(inout) :: string
 if(associated(string%chars))then
  deallocate(string%chars)
  nullify(string%chars)
 endif
 END SUBROUTINE unset_metachar


 SUBROUTINE copy_metachar(string1,string2)
!--------------------------!
! Copy STRING1 to STRING2. !
!--------------------------!
 IMPLICIT NONE
 TYPE(metachar),INTENT(in) :: string1
 TYPE(metachar),INTENT(inout) :: string2
 INTEGER len_string,ipos
 call unset_metachar(string2)
 len_string=len_metachar(string1)
 if(len_string==0)return
 allocate(string2%chars(len_string))
 do ipos=1,len_string
  string2%chars(ipos)=string1%chars(ipos)
 enddo ! ipos
 END SUBROUTINE copy_metachar


 SUBROUTINE copy_char_to_metachar(chars,string)
!-----------------------!
! Copy CHARS to STRING. !
!-----------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: chars
 TYPE(metachar),INTENT(inout) :: string
 INTEGER len_chars,ipos
 call unset_metachar(string)
 len_chars=len(chars)
 if(len_chars==0)return
 allocate(string%chars(len_chars))
 do ipos=1,len_chars
  string%chars(ipos)=chars(ipos:ipos)
 enddo ! ipos
 END SUBROUTINE copy_char_to_metachar


 SUBROUTINE copy_metachar_to_char(string,chars)
!-----------------------!
! Copy STRING to CHARS. !
!-----------------------!
 IMPLICIT NONE
 TYPE(metachar),INTENT(in) :: string
 CHARACTER(*),INTENT(out) :: chars
 INTEGER ipos,len_chars,len_string
 len_chars=len(chars)
 if(len_chars<1)return
 chars=''
 if(.not.associated(string%chars))return
 len_string=size(string%chars)
 do ipos=1,min(len_chars,len_string)
  chars(ipos:ipos)=string%chars(ipos)
 enddo ! ipos
 END SUBROUTINE copy_metachar_to_char


 INTEGER FUNCTION len_metachar(string)
!--------------------------!
! Return length of STRING. !
!--------------------------!
 IMPLICIT NONE
 TYPE(metachar),INTENT(in) :: string
 len_metachar=0
 if(.not.associated(string%chars))return
 len_metachar=size(string%chars)
 END FUNCTION len_metachar


 INTEGER FUNCTION len_trim_metachar(string)
!--------------------------------------------------!
! Return length of STRING without trailing spaces. !
!--------------------------------------------------!
 IMPLICIT NONE
 TYPE(metachar),INTENT(in) :: string
 INTEGER len_string,i0,i1
 len_trim_metachar=0
 len_string=len_metachar(string)
 if(len_string==0)return
 i0=1
 do i1=len_string,1,-1
  if(string%chars(i1)/=' ')exit
 enddo
 len_trim_metachar=max(i1-i0+1,0)
 END FUNCTION len_trim_metachar


!!$ LOGICAL FUNCTION metachar_equals_char(string,chars)
!!$!------------------------------------------------------------!
!!$! Return whether STRING represents the same string as CHARS. !
!!$!------------------------------------------------------------!
!!$ IMPLICIT NONE
!!$ TYPE(metachar),INTENT(in) :: string
!!$ CHARACTER(*),INTENT(in) :: chars
!!$ INTEGER ipos,len_string,len_chars
!!$ metachar_equals_char=.false.
!!$ len_chars=len(chars)
!!$ len_string=len_metachar(string)
!!$ if(len_chars/=len_string)return
!!$ do ipos=1,len_chars
!!$  if(string%chars(ipos)/=chars(ipos:ipos))return
!!$ enddo ! ipos
!!$ metachar_equals_char=.true.
!!$ END FUNCTION metachar_equals_char


 CHARACTER FUNCTION compare_metachar(string1,string2)
!-----------------------------------------------------------!
! Return '<', '=', '>' if STRING1 comes before, at the same !
! position as, or after STRING2 in lexicographical order,   !
! respectively.                                             !
!-----------------------------------------------------------!
 IMPLICIT NONE
 TYPE(metachar),INTENT(in) :: string1,string2
 INTEGER len_string1,len_string2,ipos
 CHARACTER c1,c2
 len_string1=len_metachar(string1)
 len_string2=len_metachar(string2)
 if(len_string1==0.or.len_string2==0)then
  compare_metachar='='
  if(len_string1/=0)then
   compare_metachar='>'
  elseif(len_string2/=0)then
   compare_metachar='<'
  endif
  return
 endif
 do ipos=1,min(len_string1,len_string2)
  c1=string1%chars(ipos)
  c2=string2%chars(ipos)
  if(c1/=c2)then
   compare_metachar='<'
   if(c1>c2)compare_metachar='>'
   return
  endif
 enddo ! ipos
 compare_metachar='='
 if(len_string1>len_string2)then
  compare_metachar='>'
 elseif(len_string1<len_string2)then
  compare_metachar='<'
 endif
 END FUNCTION compare_metachar


 SUBROUTINE read_metachar(io,string)
!-------------------------------------------------------!
! Read a full line from unit IO and place it in STRING. !
!-------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: io
 TYPE(metachar),INTENT(inout) :: string
 INTEGER,PARAMETER :: buffer_size=80
 CHARACTER(buffer_size) buffer
 INTEGER nchar_read,istatus
 call unset_metachar(string)
 do
  read(io,fmt='(a)',advance='no',eor=9999,size=nchar_read,iostat=istatus)&
   &buffer(1:buffer_size)
  if(istatus/=0)return
  call append_char_to_metachar(string,buffer)
 enddo
 9999 continue
 call append_char_to_metachar(string,buffer(1:nchar_read))
 END SUBROUTINE read_metachar


 SUBROUTINE append_char_to_metachar(string,chars)
!-------------------------!
! Append CHARS to STRING. !
!-------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: chars
 TYPE(metachar),INTENT(inout) :: string
 INTEGER len_chars,len_string,ipos
 CHARACTER,ALLOCATABLE :: temp_chars(:)
 len_chars=len(chars)
 if(len_chars==0)return
 len_string=0
 if(associated(string%chars))len_string=size(string%chars)
 if(len_string>0)then
  allocate(temp_chars(len_string))
  temp_chars(1:len_string)=string%chars(1:len_string)
 endif
 call unset_metachar(string)
 allocate(string%chars(len_string+len_chars))
 if(len_string>0)then
  string%chars(1:len_string)=temp_chars(1:len_string)
  deallocate(temp_chars)
 endif
 do ipos=1,len_chars
  string%chars(len_string+ipos)=chars(ipos:ipos)
 enddo ! ipos
 END SUBROUTINE append_char_to_metachar


 SUBROUTINE append_metachar_to_metachar(string,string1)
!---------------------------!
! Append STRING1 to STRING. !
!---------------------------!
 IMPLICIT NONE
 TYPE(metachar),INTENT(in) :: string1
 TYPE(metachar),INTENT(inout) :: string
 INTEGER len_string1,len_string
 CHARACTER,ALLOCATABLE :: temp_chars(:)
 len_string1=len_metachar(string1)
 if(len_string1==0)return
 len_string=0
 if(associated(string%chars))len_string=size(string%chars)
 if(len_string>0)then
  allocate(temp_chars(len_string))
  temp_chars(1:len_string)=string%chars(1:len_string)
 endif
 call unset_metachar(string)
 allocate(string%chars(len_string+len_string1))
 if(len_string>0)then
  string%chars(1:len_string)=temp_chars(1:len_string)
  deallocate(temp_chars)
 endif
 string%chars(len_string+1:len_string+len_string1)=string1%chars(1:len_string1)
 END SUBROUTINE append_metachar_to_metachar


 SUBROUTINE prepend_char_to_metachar(string,chars)
!--------------------------!
! Prepend CHARS to STRING. !
!--------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: chars
 TYPE(metachar),INTENT(inout) :: string
 INTEGER len_chars,len_string,ipos
 CHARACTER,ALLOCATABLE :: temp_chars(:)
 len_chars=len(chars)
 if(len_chars==0)return
 len_string=0
 if(associated(string%chars))len_string=size(string%chars)
 if(len_string>0)then
  allocate(temp_chars(len_string))
  temp_chars(1:len_string)=string%chars(1:len_string)
 endif
 call unset_metachar(string)
 allocate(string%chars(len_string+len_chars))
 do ipos=1,len_chars
  string%chars(ipos)=chars(ipos:ipos)
 enddo ! ipos
 if(len_string>0)then
  string%chars(len_chars+1:len_string+len_chars)=temp_chars(1:len_string)
  deallocate(temp_chars)
 endif
 END SUBROUTINE prepend_char_to_metachar


 SUBROUTINE prepend_metachar_to_metachar(string,string1)
!----------------------------!
! Prepend STRING1 to STRING. !
!----------------------------!
 IMPLICIT NONE
 TYPE(metachar),INTENT(in) :: string1
 TYPE(metachar),INTENT(inout) :: string
 INTEGER len_string1,len_string
 CHARACTER,ALLOCATABLE :: temp_chars(:)
 len_string1=len_metachar(string1)
 if(len_string1==0)return
 len_string=0
 if(associated(string%chars))len_string=size(string%chars)
 if(len_string>0)then
  allocate(temp_chars(len_string))
  temp_chars(1:len_string)=string%chars(1:len_string)
 endif
 call unset_metachar(string)
 allocate(string%chars(len_string+len_string1))
 string%chars(1:len_string1)=string1%chars(1:len_string1)
 if(len_string>0)then
  string%chars(len_string1+1:len_string1+len_string)=temp_chars(1:len_string)
  deallocate(temp_chars)
 endif
 END SUBROUTINE prepend_metachar_to_metachar


 FUNCTION char_metachar(string,len_chars) RESULT(chars)
!--------------------------------------------------------!
! Return the the Fortran string of length LEN_CHARS that !
! corresponds to STRING.                                 !
!--------------------------------------------------------!
 IMPLICIT NONE
 TYPE(metachar),INTENT(in) :: string
 INTEGER,INTENT(in) :: len_chars
 CHARACTER(len_chars) chars
 INTEGER len_string,ipos
 len_string=0
 chars=''
 if(.not.associated(string%chars))return
 len_string=size(string%chars)
 do ipos=1,min(len_chars,len_string)
  chars(ipos:ipos)=string%chars(ipos)
 enddo ! ipos
 END FUNCTION char_metachar


 SUBROUTINE cut_metachar(string,i0,j1)
!---------------------------------------------!
! Replace STRING with its I0:J1 substring, or !
! I0: if J1 is not present.                   !
!---------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: i0
 INTEGER,INTENT(in),OPTIONAL :: j1
 TYPE(metachar),INTENT(inout) :: string
 INTEGER len_string,len_final,i1
 TYPE(metachar) temp_string
 len_string=len_metachar(string)
 if(len_string==0)return
 i1=len_string ; if(present(j1))i1=min(j1,len_string)
 if(i0>i1)then
  call unset_metachar(string)
  return
 endif
 len_final=i1-i0+1
 allocate(temp_string%chars(len_final))
 temp_string%chars(1:len_final)=string%chars(i0:i1)
 call unset_metachar(string)
 string%chars=>temp_string%chars
 nullify(temp_string%chars)
 END SUBROUTINE cut_metachar


 SUBROUTINE unpad_metachar(string)
!-------------------------------------------------!
! Remove leading and trailing blanks from STRING. !
!-------------------------------------------------!
 IMPLICIT NONE
 TYPE(metachar),INTENT(inout) :: string
 INTEGER len_string,i0,i1
 len_string=len_metachar(string)
 if(len_string==0)return
 do i0=1,len_string
  if(string%chars(i0)/=' ')exit
 enddo
 do i1=len_string,1,-1
  if(string%chars(i1)/=' ')exit
 enddo
 call cut_metachar(string,i0,i1)
 END SUBROUTINE unpad_metachar


 INTEGER FUNCTION scan_metachar(string,c,back)
!-------------------------------------------------------------!
! Return the position of the first (or last if BACK) instance !
! of character C in STRING, or zero if not found              !
!-------------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER,INTENT(in) :: c
 TYPE(metachar),INTENT(in) :: string
 LOGICAL,INTENT(in),OPTIONAL :: back
 INTEGER len_string
 LOGICAL backwards
 scan_metachar=0
 len_string=len_metachar(string)
 if(len_string==0)return
 backwards=.false. ; if(present(back))backwards=back
 if(.not.backwards)then
  do scan_metachar=1,len_string
   if(string%chars(scan_metachar)==c)return
  enddo
 else
  do scan_metachar=len_string,1,-1
   if(string%chars(scan_metachar)==c)return
  enddo
 endif
 scan_metachar=0
 END FUNCTION scan_metachar


! Integer-to-string function


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


END MODULE casl
