!#define DEBUG
#define VERBOSE

module cross_validation

!  constants for mask

  integer, parameter :: sea  = 1
  integer, parameter :: land = 0


contains


 subroutine crossval(mask,nb,cloud_size,len,valex,X,iex,jex,mex,kex,Xex)
  implicit none

  ! argumets

  integer, intent(in)    :: mask(:,:)
  integer, intent(in)    :: nb,cloud_size

  real,    intent(in)    :: len(3)
  real,    intent(in)    :: valex
  real,    intent(inout) :: X(:,:)
  integer, intent(out)   :: iex(nb), jex(nb), mex(nb), kex(nb)
  real,    intent(out)   :: Xex(nb)

  ! local constants

  integer, parameter :: stencil = 1

  ! local variables

  real    :: prop(2*stencil+1,2*stencil+1,2*stencil+1)
  integer :: imax,jmax,kmax,N,l,nbseeds,mmax
  integer :: i,j,k,ip,jp,kp,m,nz,nzn
  logical :: already_exclued
  real    :: pmax,count,tprop,rand
  integer, allocatable :: mindex(:,:), iindex(:), jindex(:)

  imax = size(mask,1);
  jmax = size(mask,2);
  kmax = size(X,2);

  Xex = 0;

  N = imax*jmax*kmax;
  mmax = sum(mask);
  allocate(mindex(imax,jmax),iindex(mmax),jindex(mmax))

  do kp=-stencil,stencil
    do jp=-stencil,stencil
      do ip=-stencil,stencil
        prop(ip+stencil+1,jp+stencil+1,kp+stencil+1) =  &
             exp(- (ip/len(1))**2 - (jp/len(2))**2 - (kp/len(3))**2  );
      end do
    end do
  end do

  prop(2,2,2) = 0;

  nz  = 0;
  nzn = 0;

  iex = 0
  jex = 0
  kex = 0
  mex = 0

  nbseeds = nb/cloud_size

  ! maps index (i,j) to linear index m according to land-sea mask 
  m = 0
  mindex = -1

  !  change ReadMatrix
  !  do j=1,jmax
  !    do i=1,imax    

  do i=1,imax    
    do j=1,jmax
      if (mask(i,j) .eq. sea) then
        m = m+1;
        mindex(i,j) = m
        iindex(m) = i
        jindex(m) = j
      end if
    end do
  end do


  ! seeds for clouds

  nz = 0
  do while (nz .lt. nbseeds)
    m = random_integer(1,mmax)
    k = random_integer(1,kmax)

    if (X(m,k) .ne. valex) then
      nz = nz+1;
      iex(nz) = iindex(m);
      jex(nz) = jindex(m);
      kex(nz) = k;
      mex(nz) = m;
      Xex(nz) = X(m,k);
      X(m,k) = valex;       
    end if
  end do

  ! expand clouds

  do while (nz .lt. nb)
#   ifdef DEBUG
    write(6,*) 'nz ',nz
#   endif

    ! choose one clouded point

    l = random_integer(1,nz)

    i = iex(l);
    j = jex(l);
    k = kex(l);
    m = mex(l);

    ! choose one point in its neighbourhood

    kp= random_integer(max(k-stencil,1),min(k+stencil,kmax))
    jp= random_integer(max(j-stencil,1),min(j+stencil,jmax))
    ip= random_integer(max(i-stencil,1),min(i+stencil,imax))

    call random_number(rand)

    if (rand .lt. prop(ip-i+stencil+1,jp-j+stencil+1,kp-k+stencil+1) .and. mask(ip,jp) .eq. sea) then

      m = mindex(ip,jp)

      if (X(m,k) .ne. valex) then
        nz = nz+1;
        iex(nz) = ip;
        jex(nz) = jp;
        kex(nz) = k;
        mex(nz) = m;
        Xex(nz) = X(m,k);
        X(m,k) = valex;       
      end if
    end if
  end do


 end subroutine crossval

 ! returns a pseudo-random integer within [imin imax]

 integer function random_integer(imin,imax)
  implicit none
  integer, intent(in) :: imin,imax

  real :: rand
  call random_number(rand)

  random_integer = imin + (imax-imin+1)*rand

# ifdef DEBUG
  if (random_integer .gt. imax .or. random_integer .lt. imin) then
    write(6,*) 'random_integer ',random_integer,imin,imax,rand
  end if
# endif

  ! necessary due to finite precision of the product of floats

  if (random_integer.gt.imax) random_integer = imax
  if (random_integer.lt.imin) random_integer = imin

 end function random_integer




end module cross_validation


program cv
  use cross_validation
  use initfile
  use ufileformat
  implicit none

  include 'ReadMatrix.h'

  integer, pointer :: imax(:), jmax(:), first(:)
  character(len=200) :: initfilename,filename, cloudsfname,cloudmask_fname 
  character(len=200), pointer :: datafnames(:),maskfnames(:)

  integer :: mmax, nb,v, nbmots, prec, ivalex,m,n,cloud_size
  real    :: len(3) = (/ 2., 2., 0.5 /);
  integer :: i,j,k,l
  real    :: valex
  real,    pointer     :: mean(:), std(:)
  real,    pointer     :: X(:,:), Xex(:)
  !  integer, allocatable :: iex(:), jex(:), ijex(:), kex(:),  mask(:,:)
  integer, pointer :: iex(:), jex(:), ijex(:), kex(:),  mask(:,:)
  integer, allocatable :: clouds_indexes(:,:)
  logical :: isdegen
  real, allocatable  :: cloud_mask(:,:,:)

  ! external function

  integer :: iargc

  if (iargc().ne.1) then
     write(0,*) 'Usage: crossval <initfilename>'
     stop
  end if

  call getarg(1,initfilename)

  call getInitValue(initfilename,'data',datafnames)
  call getInitValue(initfilename,'mask',maskfnames)
  call getInitValue(initfilename,'clouds',cloudsfname)
  call getInitValue(initfilename,'cloud_size',cloud_size)

  v=1

  call ReadMatrix(datafnames,X,m,n,valex,maskfnames,0,imax,jmax,first,mean,std)
  write(6,*) 'readmatrix done'
  call uload(maskfnames(1),mask,ivalex)


  mmax = count(mask .eq. sea);


  ! get or compute the number of cross-validation points

  call getInitValue(initfilename,'number_cv_points',nb,default=floor(min(0.01*m*n+40,0.03*m*n)))



#ifdef VERBOSE
  write(6,*) 'Sea points according to mask: ',mmax
  write(6,*) 'Points for cross validation: ',nb
#endif

  allocate( iex(nb), jex(nb), ijex(nb), kex(nb), Xex(nb),clouds_indexes(nb,2))
  call crossval(mask,nb,cloud_size,len,valex,X,iex,jex,ijex,kex,Xex)

  clouds_indexes(:,1) = ijex;
  clouds_indexes(:,2) = kex;

  call usave(cloudsfname,real(clouds_indexes),-9999.)

  if (presentInitValue(initfilename,'cloud_mask')) then

     call getInitValue(initfilename,'cloud_mask',cloudmask_fname)

     allocate(cloud_mask(imax(1),jmax(1),n))

     cloud_mask = 1

     do l=1,nb
        cloud_mask(iex(l),jex(l),kex(l))=0
     end  do

     call usave(cloudmask_fname,cloud_mask,valex)
     deallocate(cloud_mask)

  endif

  deallocate(X)
end program cv
