#include "ppdef.h"

!--------------------------------------------------------
!----- subroutine readMatrix ----------------------------
!--------------------------------------------------------

! reads the data and mask files and creates 'matrix', 
! containing data and missing data, no land values
! 
!--------------------------------------------------------
subroutine ReadMatrix(sstfile,matrix,m,n,valex,maskfile,norma,imax,jmax,first,fileMean,fileStd)
  use ufileformat
  use initfile
  implicit none

  character(len=200),intent(in) :: sstfile(:),maskfile(:)
  real,pointer                  :: matrix(:,:)
  integer,intent(in)            :: norma
  integer,intent(out)           :: m,n
  real,intent(out)              :: valex
  integer, pointer              :: imax(:),jmax(:),first(:)
  real,pointer                  :: fileMean(:),fileStd(:)

  integer, allocatable          :: land(:),sea(:),nlines(:)
  real,pointer                  :: file(:,:,:),mask(:,:)
  real, parameter               :: valexc = 9999.
  integer                       :: i,j,s,t,k,d,q,nbvars,testN
  real                          :: Y,maskvalex,testmean
  real,pointer                  :: valexN
  character(len=200)            :: initfilename


  interface
     subroutine stat(file,mask,mean,stddev)
       real,pointer               :: file(:,:,:),mask(:,:)
       real,intent(out)           :: mean,stddev
     end subroutine stat
  end interface

  interface
     subroutine norm(file,mean,stddev)
       real,pointer               :: file(:,:,:)
       real,intent(in)            :: mean,stddev
     end subroutine norm 
  end interface


  nbvars=size(sstfile)
  
  allocate(imax(nbvars),jmax(nbvars),land(nbvars),sea(nbvars),nlines(nbvars),first(nbvars),fileMean(nbvars),fileStd(nbvars))
 
 call getarg(1,initfilename)

write(stdout,*)'initfilename',initfilename

!if( presentInitValue(initfilename,'mask')) then
!if(maskfile = 'nomaskfile') then
! write(stdout,*)'no mask for reading '
!endif

  do q=1,nbvars


     call uload(sstfile(q),file, valex)
    
     if( presentInitValue(initfilename,'mask')) then
       call uload(maskfile(q),mask,maskvalex)
     else
       allocate(mask(size(file,1),size(file,2)))
       mask = mask+1
       
     endif

     where (mask.eq.0) mask = valexc
     

     !where (spread(mask,3,size(file,3)).eq.valexc) file = valex

     do t=1,size(file,3)
        do j=1,size(file,2)
           do i=1,size(file,1)
             if (mask(i,j).eq.valexc) file(i,j,t) = valexc
              
           enddo
        enddo
              
     enddo
     

     !    ------------------------------------
     !     check if all input matrices
     !     have the same size as their masks
     !    ------------------------------------

     if(size(file,1).ne.size(mask,1).or.size(file,2).ne.size(mask,2)) then
        write(stdout,*) 'File', trim(sstfile(q)), 'does not have the same size as its mask,', trim(maskfile(q)), 'please check.'
        STOP
     end if

     imax(q)=size(mask,1)
     jmax(q)=size(mask,2)

     !    -----------------------------------------
     !     temporal dimension, same for all files!!
     !    -----------------------------------------  
     nlines(q)=size(file,3)
    
     !    -----------------------------------------
     !     counts land points
     !    -----------------------------------------    
     land(q)=count(mask.eq.valexc)

     !    -----------------------------------------
     !     calculates first dimension of the matrix 
     !     (spatial points)
     !    -----------------------------------------    
     sea(q)=(imax(q)*jmax(q))-land(q)



     if(q.eq.1) then
        first(q)=1
     else
        first(q)=first(q-1)+sea(q-1)
     end if

     deallocate(file,mask)


  end do

  if(any(nlines.ne.nlines(1))) then
     write(stdout,*) 'The temporal sizes of both files provided are not the same'
     write(stdout,*) 'I have found the following temporal dimensions,', nlines(:)
     write(stdout,*) 'Please bear in mind that this routine is constructed for the'
     write(stdout,*) 'case when the temporal size of all files being reconstructed is the same.'


     STOP
  end if




  write(stdout,*)
  write(stdout,*) '********************************************************************'
  write(stdout,*) 'Now some statistics about your data:'
  write(stdout,*)
  
  do q=1,nbvars
  write(stdout,400)' Number of mask land points: ',land(q)
  write(stdout,300)' Dimension of file ',q,': ','',imax(q),'x',jmax(q),'x',nlines(q)
  write(stdout,*)  
  enddo

  m=sum(sea(:))
  allocate(matrix(m,nlines(1)))
      write(*,*)'allocate'
  do q=1,nbvars

    
     call uload(sstfile(q),file, valex)
     
 write(*,*)'end uload'
!     where (file.eq.valex) file = valexc 

 do t=1,size(file,3)
       do j=1,size(file,2)
         do i=1,size(file,1)
           if (file(i,j,t).eq.valex) file(i,j,t)=valexc
        enddo
     enddo
  enddo


    write(*,*)'end uload2'
     if( presentInitValue(initfilename,'mask')) then
       call uload(maskfile(q),mask,maskvalex)
     else
       allocate(mask(size(file,1),size(file,2)))
       mask = mask+1
     endif
write(*,*)'end uload mask'

!     where (mask.eq.0) mask = valexc

  do j=1,size(mask,2)
       do i=1,size(mask,1)
           if (mask(i,j).eq.0) mask(i,j)=valexc
        enddo
  enddo
  write(*,*)'end loop'


     !    -----------------------------------
     !     normalisation of the input matrices
     !     if norma = 1
     !    ------------------------------------
      
     if(norma.eq.1) then
        call stat(file,mask,fileMean(q),fileStd(q))
        call norm(file,fileMean(q),fileStd(q))
        !fileNorm=file
     else
         
        call stat(file,mask,fileMean(q),fileStd(q))
       ! fileNorm = file
       ! deallocate(file)
        
     end if
      write(*,*)'end loop2'


     !where(fileNorm.eq.9999) fileNorm=valexc 
     do t=1,size(file,3)
       do j=1,size(file,2)
         do i=1,size(file,1)
           if (file(i,j,t).eq.9999) file(i,j,t)=valexc
        enddo
     enddo
  enddo

write(*,*)'end stats'
     ! where (spread(mask,3,nlines(q)).eq.valexc) file = valexc

     ! testN = count(file.ne.valexc)
     ! testmean = sum(file,file.ne.valexc)/testN

     !    ---------------------------------------
     !     reshape input matrices into one matrix
     !     size M*N (M spatial dimension,
     !               N temporal dimension)
     !    ---------------------------------------
   
     do t=1,nlines(q)

        s = first(q)

        do i=1,imax(q)
           do j=1,jmax(q)
              if (mask(i,j).eq.valexc) file(i,j,t) = valexc

              if(mask(i,j).eq.valexc.and.file(i,j,t).ne.valexc) then
                 write(stdout,*) 'WARNING ',i,j
              end if


              
              if(mask(i,j).ne.valexc) then

                 matrix(s,t)=file(i,j,t)

                 s = s+1

              endif
           enddo
        enddo
     enddo
write(*,*)'end'
     deallocate(file,mask)

write(*,*)'end2'
call flush(stdout,istat)
  !    -----------------------------------
  !     write mean and standard deviation
  !     of input matrices
  !    ------------------------------------

  write(stdout,*)
  write(stdout,500)' Mean: ',fileMean(q)
  write(stdout,500)' Standard deviation: ',fileStd(q)
  write(stdout,*)


  enddo

  valex=valexc

  n=nlines(1)




deallocate(land,sea,nlines)
200 format(a20,f6.2)
300 format(a26,i2,a2,a15,i6,a3,i6,a3,i6)
400 format(a30,i39)
500 format(a30,f39.4)

end subroutine ReadMatrix





