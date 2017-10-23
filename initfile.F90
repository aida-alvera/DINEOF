!
!  OAK, Ocean Assimilation Kit
!  Copyright(c) 2002-2011 Alexander Barth and Luc Vandenblucke
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!

! include the fortran preprocessor definitions
#include "ppdef.h"

module InitFile

  integer, parameter :: maxLength = 1024*4

  interface getInitValue
     module procedure getInitValue_integer,   &
          getInitValue_real,      &
          getInitValue_string,    &
          getInitValue_integer1D, &
          getInitValue_real1D,    &
          getInitValue_string1D
  end interface

  interface presentInitValue
     module procedure presentInitValue
  end interface

  character, parameter :: tab = char(9)

!  integer :: stderr=6
contains


  !--------------------------------------------------------------------

integer function  freeunit()
integer unit
integer, parameter :: iu = 1000
logical used

  used = .true.
  unit = iu

  do while (used)
    unit = unit+1
    inquire(unit,opened=used)
!    if (used) write(stdout,*) unit,' used'
  end do

!  write(stdout,*) unit
  freeunit = unit
end function

  !--------------------------------------------------------------------

  logical function iswhitespace(c)
   character, intent(in) :: c
   iswhitespace = c == ' ' .or. c == tab
  end function iswhitespace

  !--------------------------------------------------------------------

  subroutine parseItem(buffer,first_index,last_index,error)
    character(len=*), intent(in) :: buffer
    integer, intent(inout) :: first_index
    integer, intent(out) :: last_index
    integer, intent(out) :: error

    error = 0

    ! eat blanks

    do while (iswhitespace(buffer(first_index:first_index)))
       first_index = first_index + 1

       if (first_index.gt.len(buffer)) exit
    enddo

    if (first_index.le.len(buffer)) then
       select case (buffer(first_index:first_index))
       ! comments
       case ('!','#')
          first_index = len(buffer)+1
          last_index = len(buffer)

       case ('=',',',';','*','/')
          last_index = first_index

       case (']')
          error = -1

       case ('[')
          last_index = first_index+1
          do while ((last_index.le.len(buffer)).and.  &
               (buffer(last_index:last_index).ne.']'))
             last_index = last_index + 1  
          end do
          if (last_index.gt.len(buffer)) then
             error = -1
          end if

       case ('''')
          ! for strings
          last_index = first_index+1
          do while ((last_index.le.len(buffer)).and.  &
               .not.((last_index.eq.len(buffer)).and.(buffer(last_index:last_index).eq.'''')).and. &
               .not.((buffer(last_index:last_index).eq.'''').and.(buffer(last_index+1:last_index+1).ne.'''')))
             if ((buffer(last_index:last_index).eq.'''').and.(buffer(last_index+1:last_index+1).eq.'''')) then
                last_index = last_index + 2
             else
                last_index = last_index + 1  
             end if
          end do
          if (last_index.gt.len(buffer)) then
             error = -1
          end if

       case default
          last_index = first_index

          do while ((last_index.lt.len(buffer)).and.  &
               (buffer(last_index+1:last_index+1).ne.'=').and. &
               (buffer(last_index+1:last_index+1).ne.'!').and. &
               (buffer(last_index+1:last_index+1).ne.'''').and. &
               (buffer(last_index+1:last_index+1).ne.'[').and. &
               (buffer(last_index+1:last_index+1).ne.',').and. &
               (buffer(last_index+1:last_index+1).ne.';').and. &
               (buffer(last_index+1:last_index+1).ne.']').and. &
               (.not.iswhitespace(buffer(last_index+1:last_index+1))))
             last_index = last_index + 1  
          enddo
       end select
    else
       last_index = len(buffer)
    endif

  end subroutine parseItem

  !--------------------------------------------------------------------

  subroutine parseString(source,dest)
    implicit none
    character(len=*), intent(in) :: source
    character(len=*), intent(out) :: dest

    integer :: position
    integer :: length

    position = 1
    length = 0
    dest = ''

    do while ((position.le.len(source)).and.(length.lt.len(dest)))
       dest = dest(1:length)//source(position:position)
       length = length+1

       if (source(position:position).eq.'''') then
          position = position +2
       else
          position = position + 1
       endif
    end do
  end subroutine parseString

  !--------------------------------------------------------------------

  function countElements(string,error)
    integer countElements
    character (len=*), intent(in) :: string
    integer, intent(out) :: error

    integer :: length, first_index, last_index, ParseError, count


    length = len_trim(string)
    error = 0

    if ((string(1:1).eq.'[').and.(string(length:length).eq.']')) then
       length = len_trim(string(1:length-1))
       first_index = 2

       call parseItem(string,first_index,last_index,ParseError)    
       if ((ParseError.ne.0).or.(first_index.gt.last_index)) then
          error = -1
          return
       endif
       count = 1

       do while (last_index.ne.length)
          first_index = last_index+1
          call parseItem(string,first_index,last_index,ParseError)    

          if ((ParseError.ne.0).or.(first_index.gt.last_index).or. &
               (string(first_index:last_index).ne.',')) then
             error = -1
             return
          else
             count = count+1
          end if

          first_index = last_index+1
          call parseItem(string,first_index,last_index,ParseError)    

          if ((ParseError.ne.0).or.(first_index.gt.last_index)) then
             error = -1
             return
          endif
       end do
    else
       error = -1
    end if

    countElements = count
  end function countElements

  !________________________________________________________________________________
  !

  subroutine hashVec(string,array,error,errormsg)
    implicit none
    character(len=*), intent(in) :: string
    character(len=maxLength), pointer :: array(:)
    integer, intent(out) :: error
    character(len=*), intent(out) :: errormsg

    integer :: dim, ParseError,fi,li,i,istat

    dim = countElements(string,ParseError)
    if (ParseError.ne.0) then
       write(errormsg,*) 'hashArray: parsing error in string "',trim(string),'"'
       error = -2
       return
    end if

    fi = 1; 
    call parseItem(string,fi,li,ParseError)    

    if (ParseError.ne.0) then
       write(errormsg,*) 'hashArray: parsing error in string "',trim(string),'"'
       error = -2
       return
    end if

    allocate(array(dim))

    fi = fi+1
    call parseItem(string,fi,li,ParseError)    
    array(1) = string(fi:li)

    do i=2,dim
       fi=li+1
       call parseItem(string,fi,li,ParseError)    

       fi=li+1
       call parseItem(string,fi,li,ParseError)    
       array(i) = string(fi:li)
    end do
  end subroutine hashVec

  !________________________________________________________________________________
  !
  logical function domatch(str,pattern)
  implicit none
  character(*), intent(in) :: str,pattern

# ifdef GMATCH
  integer ::  match
   external match
  domatch = match(trim(str)//char(0),trim(pattern)//char(0)).eq.1    
# else
  domatch = trim(str).eq.trim(pattern)
#endif
  end function
  !________________________________________________________________________________
  !


  subroutine getInitString(filename,variable,string,lnNumber,error,errormsg)
    implicit none
    character(len=*), intent(in) :: filename, variable
    character(len=*), intent(out) :: string
    integer, intent(out) :: lnNumber, error
    character(len=*), intent(out) :: errormsg

    ! error = 0   every thing is fine
    ! error = -1  variable not found
    ! error = -2  syntax error in file


    character(len=maxLength) :: line
    integer iostat, lineNumber, ParseError
    integer fi,li, iu,istat
#ifdef _OPENMP
    integer omp_get_thread_num
#endif
    error = -1
    lineNumber = 1

!$omp critical (reserveFreeUnit)
    iu = freeunit()

!    iu=20234+omp_get_thread_num()
!    write(6,*) 'open iu',iu,omp_get_thread_num()
    open(iu,file=filename,status='old')
!    write(6,*) 'opened iu',iu,omp_get_thread_num()

    read(iu,'(A)',iostat=iostat) line

    do while (iostat.eq.0) 
       fi = 1
       call parseItem(line,fi,li,ParseError)
       if (ParseError.ne.0) then
          write(errormsg,*) 'getInitString: parsing error at line number ', &
               lineNumber,' in file "',trim(filename),'"'
          error = -2
          exit
       end if

!       if ((fi.le.len(line)).and.(line(fi:li).eq.trim(variable))) then
       if ((fi.le.len(line)).and.(domatch(trim(variable),line(fi:li)))) then
          error = 0
          lnNumber = lineNumber

          fi = li+1
          call parseItem(line,fi,li,ParseError)
          if (ParseError.ne.0) then
             write(errormsg,*) 'getInitString: parsing error at line number ', &
                  lineNumber,' in file "',trim(filename),'"'
             error = -2
             exit
          end if

          if ((fi.gt.len(line)).or.(line(fi:li).ne.'=')) then
             write(errormsg,*) 'getInitString: Missing "=" at line number ', &
                  lineNumber,' in file "',trim(filename),'"'
             error = -2
             exit
          end if

          fi = li+1
          call parseItem(line,fi,li,ParseError)
          if (ParseError.ne.0) then
             write(errormsg,*) 'getInitString: parsing error at line number ', &
                  lineNumber,' in file "',trim(filename),'"'
             error = -2
             exit
          end if

          string = line(fi:li)
       end if

       read(iu,'(A)',iostat=iostat) line
       lineNumber = lineNumber + 1
    end do

!!$omp critical (reserveFreeUnit)
    close(iu)
!$omp end critical (reserveFreeUnit)

!!$omp end critical

    if (error.eq.-1) then
       write(errormsg,*) 'getInitString: variable "',trim(variable),'" not found ', &
            ' in file "',trim(filename),'"'
    end if

  end subroutine getInitString

  !________________________________________________________________________________
  !

  subroutine getInitValue_integer(filename,variable,value,err,default)
    implicit none
    character(len=*), intent(in) :: filename, variable
    integer, intent(out) :: value
    integer, intent(out), optional :: err
    integer, intent(in), optional :: default

    ! error = 0   every thing is fine
    ! error = -1  variable not found
    ! error = -2  syntax error in file


    character(len=maxLength) :: string
    character(len=maxLength) :: errormsg
    integer lineNumber, ConversionError, fi,li,iu, error, istat

    call getInitString(filename,variable,string,lineNumber,error,errormsg)

    if (error.eq.0) then
       read(string,*,iostat=ConversionError) value

       if (ConversionError.ne.0) then
          write(errormsg,*) 'getInitValue_integer: error at line ',lineNumber, &
             ' in file "',trim(filename),'". Integer expected but "', &
             trim(string),'" found.'
          error = -2
       end if
    end if

    if (error.eq.-1.and.present(default)) then
      value=default
      error = 0
    end if

    ! error handling

    if (present(err)) then
    ! if an error occurs, it will be treated at the calling level
       err = error
    else
       if (error.ne.0) then
         write(stderr,*) trim(errormsg)
         stop 'Program stoped in initfile.F90 around line 361'
       end if
    end if
  end subroutine getInitValue_integer

  !________________________________________________________________________________
  !

  subroutine getInitValue_real(filename,variable,value,err,default)
    implicit none
    character(len=*), intent(in) :: filename, variable
    real, intent(out) :: value
    integer, intent(out), optional :: err
    real, intent(in), optional :: default

    ! error = 0   every thing is fine
    ! error = -1  variable not found
    ! error = -2  syntax error in file


    character(len=maxLength) :: string
    character(len=maxLength) :: errormsg
    integer lineNumber, ConversionError, fi,li,iu, error,istat

    call getInitString(filename,variable,string,lineNumber,error,errormsg)

    if (error.eq.0) then
       read(string,*,iostat=ConversionError) value

       if (ConversionError.ne.0) then
          write(errormsg,*) 'getInitValue_real: error at line ',lineNumber, &
             ' in file "',trim(filename),'". Real expected but "', &
             trim(string),'" found.'
          error = -2
       end if
    end if

    if (error.eq.-1.and.present(default)) then
      value=default
      error = 0
    end if

    ! error handling

    if (present(err)) then
    ! if an error occurs, it will be treated at the calling level
       err = error
    else
       if (error.ne.0) then
         write(stderr,*) trim(errormsg)
         stop 'Program stoped in initfile.F90 around line 412'
       end if
    end if

  end subroutine getInitValue_real

  !________________________________________________________________________________
  !

  subroutine getInitValue_string(filename,variable,value,err,default)
    implicit none
    character(len=*), intent(in) :: filename, variable
    character(len=*), intent(out) :: value
    integer, intent(out), optional :: err
    character(len=*), intent(in), optional :: default

    ! error = 0   every thing is fine
    ! error = -1  variable not found
    ! error = -2  syntax error in file


    character(len=maxLength) :: string
    character(len=maxLength) :: errormsg
    integer lineNumber, ConversionError, fi,li,iu, error,istat

    call getInitString(filename,variable,string,lineNumber,error,errormsg)

    if (error.eq.0) then
       read(string,*,iostat=ConversionError) value

       if (ConversionError.ne.0) then
          write(errormsg,*) 'getInitValue_string: error at line ',lineNumber, &
             ' in file "',trim(filename),'". Character string expected but "', &
             trim(string),'" found.'
          error = -2
       end if
    end if


    if (error.eq.-1.and.present(default)) then
      value=default
      error = 0
    end if

!    write(6,*) 'variable,value ',variable,value

    ! error handling

    if (present(err)) then
    ! if an error occurs, it will be treated at the calling level
       err = error
    else
       if (error.ne.0) then
         write(stderr,*) trim(errormsg)
         stop 'Program stoped in initfile.F90 around line 463'
       end if
    end if
  end subroutine getInitValue_string

  !________________________________________________________________________________
  !

  subroutine getInitValue_integer1D(filename,variable,value,err,default)
    implicit none
    character(len=*), intent(in) :: filename, variable
    integer, pointer :: value(:)
    integer, intent(out), optional :: err
    integer, intent(in), optional :: default(:)

    ! error = 0   every thing is fine
    ! error = -1  variable not found
    ! error = -2  syntax error in file


    character(len=maxLength) :: string
    character(len=maxLength) :: errormsg
    character(len=maxLength), pointer :: vec(:)
    integer lineNumber, ConversionError, error, i,istat

    call getInitString(filename,variable,string,lineNumber,error,errormsg)

    if (error.eq.0) then
       call hashVec(string,vec,error,errormsg)

       if (error.eq.0) then
          allocate(value(size(vec)))

          do i=1,size(vec)
             read(vec(i),*,iostat=ConversionError) value(i)

             if (ConversionError.ne.0) then
                write(errormsg,*) 'getInitValue_integer1D: error at line ', &
                     lineNumber,' in file "',trim(filename),'"'
                error = -2
                exit
             end if
          end do

          deallocate(vec)
       end if
    end if

    if (error.eq.-1.and.present(default)) then
      allocate(value(size(default)))
      value=default
      error = 0
    end if

    ! error handling

    if (present(err)) then
    ! if an error occurs, it will be treated at the calling level
       err = error
    else
       if (error.ne.0) then
         write(stderr,*) trim(errormsg)
         stop 'Program stoped in initfile.F90 around line 525'
       end if
    end if

  end subroutine getInitValue_integer1D

  !________________________________________________________________________________
  !

  subroutine getInitValue_real1D(filename,variable,value,err,default)
    implicit none
    character(len=*), intent(in) :: filename, variable
    real, pointer :: value(:)
    integer, intent(out), optional :: err
    real, intent(in), optional :: default(:)

    ! error = 0   every thing is fine
    ! error = -1  variable not found
    ! error = -2  syntax error in file


    character(len=maxLength) :: string
    character(len=maxLength) :: errormsg
    character(len=maxLength), pointer :: vec(:)
    integer lineNumber, ConversionError, error, i,istat

    call getInitString(filename,variable,string,lineNumber,error,errormsg)

    if (error.eq.0) then
       call hashVec(string,vec,error,errormsg)

       if (error.eq.0) then
          allocate(value(size(vec)))

          do i=1,size(vec)
             read(vec(i),*,iostat=ConversionError) value(i)

             if (ConversionError.ne.0) then
                write(errormsg,*) 'getInitValue_real1D: error at line ', &
                     lineNumber,' in file "',trim(filename),'"'
                error = -2
                exit
             end if
          end do

          deallocate(vec)
       end if
    end if

    if (error.eq.-1.and.present(default)) then
      allocate(value(size(default)))
      value=default
      error = 0
    end if

    ! error handling

    if (present(err)) then
    ! if an error occurs, it will be treated at the calling level
       err = error
    else
       if (error.ne.0) then
         write(stderr,*) trim(errormsg)
         stop 'Program stoped in initfile.F90 around line 525'
       end if
    end if

  end subroutine getInitValue_real1D

  !________________________________________________________________________________
  !

  subroutine getInitValue_string1D(filename,variable,value,err,default)
    implicit none
    character(len=*), intent(in) :: filename, variable
    character(len=*), pointer :: value(:)
    integer, intent(out), optional :: err
    character(len=*), intent(in), optional :: default(:)

    ! error = 0   every thing is fine
    ! error = -1  variable not found
    ! error = -2  syntax error in file


    character(len=maxLength) :: string
    character(len=maxLength) :: errormsg
    character(len=maxLength), pointer :: vec(:)
    integer lineNumber, ConversionError, error, i,istat

    call getInitString(filename,variable,string,lineNumber,error,errormsg)

    if (error.eq.0) then
       call hashVec(string,vec,error,errormsg)

       if (error.eq.0) then
          allocate(value(size(vec)))

          do i=1,size(vec)
             read(vec(i),*,iostat=ConversionError) value(i)

             if (ConversionError.ne.0) then
                write(errormsg,*) 'getInitValue_string1D: error at line ', &
                     lineNumber,' in file "',trim(filename),'"'
                error = -2
                exit
             end if
          end do

          deallocate(vec)
       end if
    end if

    if (error.eq.-1.and.present(default)) then
      allocate(value(size(default)))
      value=default
      error = 0
    end if

    ! error handling

    if (present(err)) then
    ! if an error occurs, it will be treated at the calling level
       err = error
    else
       if (error.ne.0) then
         write(stderr,*) trim(errormsg)
         stop 'Program stoped in initfile.F90 around line 655'
       end if
    end if

  end subroutine getInitValue_string1D

  !________________________________________________________________________________
  !

  logical function presentInitValue(filename,variable,err)
    implicit none
    character(len=*), intent(in) :: filename, variable
    integer, intent(out), optional :: err

    ! error = 0   every thing is fine
    ! error = -2  syntax error in file


    character(len=maxLength) :: string
    character(len=maxLength) :: errormsg
    integer lineNumber, ConversionError, error, i

    call getInitString(filename,variable,string,lineNumber,error,errormsg)

    if (error.eq.0) then
       presentInitValue = .true.
    elseif (error.eq.-1) then
       presentInitValue = .false.
    else
      if (present(err)) then
        err = error
      else
         write(stderr,*) trim(errormsg)
         stop 'Program stoped in initfile.F90 around line 688'
      end if
    end if
  end function

  !________________________________________________________________________________
  !




end module InitFile




