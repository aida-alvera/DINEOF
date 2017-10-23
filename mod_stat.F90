module mod_stat
contains

  subroutine stat(file,mask,fileNorm,mean,stddev)
    use ufileformat
    implicit none

    real,pointer               :: file(:,:,:),mask(:,:)
    real,pointer               :: fileNorm(:,:,:)
    real,intent(out)           :: mean,stddev
    real,parameter             :: valexNorm = 9999.
    real                       :: var  
    integer                    :: N,i,j,t


    allocate(fileNorm(size(file,1),size(file,2),size(file,3)))


    ! where (spread(mask,3,size(file,3)).eq.valexNorm) file = valexNorm
    ! 
    ! probably more efficient

    do t=1,size(file,3)
       do j=1,size(file,2)
          do i=1,size(file,1)
             if (mask(i,j).eq.valexNorm) file(i,j,t) = valexNorm
          enddo
       enddo
    enddo

    N = count(file.ne.valexNorm)
    mean = sum(file,file.ne.valexNorm)/N
    var = sum((file-mean)**2,file.ne.valexNorm)/N
    stddev = sqrt(var)
    ! l = maxloc(file,file.ne.valex)
    ! m = minloc(file,file.ne.valex)


    where(file.ne.valexNorm) fileNorm=(file-mean)/stddev
    where(file.eq.valexNorm) fileNorm=valexNorm


  end subroutine stat
end module mod_stat
