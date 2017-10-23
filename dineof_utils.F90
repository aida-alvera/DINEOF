#include "ppdef.h"

module dineof_utils

contains

!---------------------------------------------------------
 !----- subroutine valsvd ---------------------------------
 !---------------------------------------------------------

 ! computes the missing data values from the EOF basis     
 ! 
 !---------------------------------------------------------

subroutine valsvd(U,S,V,X,IEX,JEX,VAL,NN,IMISST) 

implicit none
 integer,intent(in)    :: IMISST
 real, intent(in)      :: U(:,:),V(:,:),S(:,:)
 real, intent(inout)   :: X(:,:)
 integer,intent(in)    :: NN,IEX(:),JEX(:)
 real,intent(out)      :: VAL
 integer               :: K,t


do t=1,IMISST
  VAL=0 
  DO K=1,NN 

    VAL=VAL+U(IEX(t),K)*S(K,1)*V(JEX(t),K) 
    
  ENDDO
  X(IEX(t),JEX(t))=VAL
enddo
 

end subroutine




 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine dindiff(x,B,alpha,numit)
  use ufileformat
  implicit none

  real,intent(inout)         :: B(:)
  real,intent(in)            :: x(:)
  real, pointer              :: BF(:),xe(:)
  integer                    :: nsize,k,numit
  real                       :: alpha,valex


  nsize = size(B,1)

  allocate(BF(size(B,1)+1))

  allocate(xe(size(x,1)+1))
  BF=0

  ! extended x
  xe(1) = 1.5*x(1) - .5 * x(2)
  xe(2:nsize) = (x(1:nsize-1) + x(2:nsize))/2
  xe(nsize+1) = 1.5*x(nsize) - .5 * x(nsize-1)


  do k=1,numit
    !  F(2:nsize) = alpha * (f(2:nsize) - f(1:nsize-1));
    !  f = f +  (F(2:nsize+1) - F(1:nsize));  
    BF(2:nsize) = alpha * (B(2:nsize) - B(1:nsize-1))/(x(2:nsize) - x(1:nsize-1))
    B = B + (BF(2:nsize+1) - BF(1:nsize))/ (xe(2:nsize+1) - xe(1:nsize))
  end do

  deallocate(BF,xe)

 end subroutine dindiff

 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


 !---------------------------------------------------------
 !----- subroutine writeMatrix ----------------------------
 !---------------------------------------------------------

 ! reshapes matrix into 3D file, then writes results  
 ! into separate files
 ! 
 !---------------------------------------------------------

 subroutine writeMatrix(X,xmean,valex,resultfnames,maskfile,norma,imax,jmax,first,DirOutput,fileMean,fileStd)
  use ufileformat
  use initfile
  implicit none

  real,intent(in)                  :: X(:,:)
  real,intent(in),optional         :: fileMean(:),fileStd(:)
  integer,intent(in)               :: imax(:),jmax(:),first(:)
  real,intent(in)                  :: xmean
  real,intent(inout)               :: valex
  integer,intent(in)               :: norma

  character(len=200),intent(in)    :: resultfnames(:),maskfile(:)
  character(len=200),intent(in)    :: DirOutput

  real, allocatable                :: sst(:,:,:)
  real,pointer                     :: mask(:,:)
  real, parameter                  :: valexc = 9999.
  real                             :: mean,stddev,var
  integer                          :: i,j,s,t,q,nbvars,NN,N
  character(len=200)               :: initfilename


  N = size(X,2)

  nbvars=size(resultfnames)
  call getarg(1,initfilename)

  !Reshape the matrix into a 3D matrix (two spatial dimensions and one temporal dimension)
  !---------------------------------------------------------------------------------------
  do q=1,nbvars

    allocate(sst(imax(q),jmax(q),N))

    if( presentInitValue(initfilename,'mask')) then
      call uload(maskfile(q),mask,valex)
    else
      allocate(mask(size(sst,1),size(sst,2)))
      mask = mask+1 
    endif

      where (mask.eq.0) mask = valexc

      do t=1,N

      s = first(q)

      do i=1,imax(q)
        do j=1,jmax(q)

          if(mask(i,j).ne.valexc) then

            sst(i,j,t)=X(s,t);

            if(norma.eq.1) then
              sst(i,j,t)=(sst(i,j,t)*fileStd(q))+fileMean(q)
            end if

            s = s+1
          else
            sst(i,j,t)=valexc;
          endif
        enddo
      enddo
    enddo



    NN = count(sst.ne.valexc)
    mean = sum(sst,sst.ne.valexc)/NN
    var = sum((sst-mean)**2,sst.ne.valexc)/NN
    stddev = sqrt(var)


    if(norma.eq.1) then
      !   -----------------------------------------
      !   some statistics about the filled matrices
      !   -----------------------------------------

      write(stdout,*)
      write(stdout,*)'mean ',mean
      write(stdout,*)'Standard deviation',stddev,fileStd,fileMean
      write(stdout,*)
    end if

    call usave(trim(resultfnames(q)),sst(:,:,:),valexc)

    deallocate(sst,mask)

  enddo

 end subroutine writeMatrix



 !---------------------------------------------------------
 !----- subroutine writeSVD ----------------------------
 !---------------------------------------------------------

 ! writes EOFs
 ! 
 !---------------------------------------------------------


 subroutine writeSVD(DirOutput,eofvfname,eofsfname,s,v,valc,sumVarEnd)
  use ufileformat
  use netcdf
  implicit none

  character(len=200),intent(in)     :: DirOutput,eofvfname,eofsfname
  real,intent(in)                   :: s(:),v(:,:),valc(:),sumVarEnd

  integer                           :: i,j,k,ncid,status
  real, parameter                   :: valex = 9999.
  real, parameter                   :: valexc = 9999.
  real                              :: varex(size(s))


  ! call usave(trim(resultfnames(q)),sst(:,:,:),valexc)
   status = nf90_create(path = trim(DirOutput)//'/DINEOF_diagnostics.nc', cmode = nf90_clobber, ncid = ncid)
  if (status /= nf90_noerr) call handle_err(status)
  
  status = nf90_close(ncid)
  if (status /= nf90_noerr) call handle_err(status)


  !   Variance detailed
  !   _______________

  varex= 100.0*s**2/sumVarEnd


  call usave(trim(DirOutput)//'/DINEOF_diagnostics.nc#varEx',varex,valex)

  write(stdout,*) 'Sum of the squares of the singular values of the ',size(s), 'EOFs retained', sum(varex)

  !
  !   Singular values
  !   _______________


  call usave(trim(DirOutput)//'/DINEOF_diagnostics.nc#vlsng',s,valex)

 
  !
  !   Temporal modes 
  !   _______________
  
  call usave(eofvfname,v,valex)

  !---------------------------------------------------------------------------------------

  !   Valc: write expected error for each mode
  !   ----------------
  call usave(trim(DirOutput)//'/DINEOF_diagnostics.nc#valc',valc,valexc)   
  !---------------------------------------------------------------------------------------


  !stop      

500 format(a72)
501 format(e10.4)
502 format(3e10.4)
503 format(300(e10.4,1X))
504 format(a5,i3,a3,f8.4,a2)

 end subroutine writeSVD


 subroutine clear_files(filenames)
 use netcdf
 character(len=200),intent(in)    :: filenames
 integer                          :: status,ncid



  status = nf90_create(path = trim(filenames), cmode = nf90_clobber, ncid = ncid)
  if (status /= nf90_noerr) call handle_err(status)
  
  status = nf90_close(ncid)
  if (status /= nf90_noerr) call handle_err(status)

 
 end subroutine clear_files



 subroutine handle_err(status)
 use netcdf
        integer, intent ( in) :: status
     
        if(status /= nf90_noerr) then
          print *, trim(nf90_strerror(status))
          stop "Stopped"
        end if
 end subroutine handle_err




end module dineof_utils
