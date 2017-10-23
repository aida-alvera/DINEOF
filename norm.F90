#include "ppdef.h"

subroutine norm(file,mean,stddev)
  use ufileformat
  implicit none

  real,pointer               :: file(:,:,:)
  real,intent(in)            :: mean,stddev
  real,parameter             :: valexNorm = 9999.
  integer                    :: i,j,t



  do t=1,size(file,3)
     do j=1,size(file,2)
        do i=1,size(file,1)
          if (file(i,j,t).ne.valexNorm) then 
            file(i,j,t)=(file(i,j,t)-mean)/stddev
          else
            file(i,j,t)=valexNorm
          endif
        enddo
     enddo
  enddo

end subroutine norm
