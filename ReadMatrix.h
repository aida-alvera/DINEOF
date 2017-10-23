    interface
      subroutine ReadMatrix(sstfile,matrix,m,n,valex,maskfile,norma,imax,jmax,first,fileMean,fileStd)
      character(len=200),intent(in)  :: sstfile(:),maskfile(:)
      real,pointer                  :: matrix(:,:),fileMean(:),fileStd(:)
      integer,intent(out)           :: m,n
      integer,intent(in)            :: norma
      real,intent(out)              :: valex
      integer, pointer              :: imax(:),jmax(:),first(:)
      end subroutine ReadMatrix
    end interface
