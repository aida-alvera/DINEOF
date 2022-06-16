#include "ppdef.h"

program dineof
 use ufileformat
 use initfile
 use dineof_utils

 implicit none


 integer              ::   maxm, maxn, maxnev, maxncv, ldv, ldu,neini,nev,ncv,&
      &                        M,N,I,J,IMISS,IMEAN,ICROSS,IMISST,IM,k,IT,P,Popt,NN,&
      &                        NITE,ipid, ierror,R,nitemax,q,getpid,&
      &                        MAKENIT,ival,rec,nbvars,eofno,testN,norma,nitedonelast

 integer,pointer      :: imax(:),jmax(:),first(:)

! parameter               (maxnev=1000, maxncv=1010)

 real, allocatable    :: meanA(:),& 
      &                        w(:), V(:,:), U(:,:),&
      &                        workl(:), workd(:),&
      &                        S(:,:), resid(:),xcr(:,:),X2(:,:)
 real                 :: tol, t, sumSV, XMEAN, VALEX,valextime,Y,DIFF,VAL,VALCopt,valexmask,Xrec
 real                 :: t1,t2,t3,mytol,testmean,sumVarIni,SumVarEnd,varini,stdvini
 real, pointer        :: X(:,:),mask(:,:),Xp(:,:),fileMean(:),fileStd(:)

 ! indexes of cross-validation points choosen by user
 real, pointer        :: cvp_indexes(:,:)
 real                 :: valexc,toliter,cfrac

 real, allocatable    :: XEX(:,:),VALC(:),sst(:,:,:),img(:,:,:),valosc(:,:),Xlastit(:,:),nitedonef(:),valosclast(:)
 real,allocatable     :: VALC_all(:),VALC_present(:)
 integer, allocatable :: IEX(:),JEX(:),nitedone(:)

 logical,allocatable  :: select(:)
 logical              :: istrans
 character*72            comment

 character(len=200),pointer   :: sstfile(:),maskfile(:),resultfnames(:),eofufnames(:)
 character(len=200)   :: DirOutput, initfilename, cloudsfname,eofvfname,eofsfname,logfile,name1,name2,name3,name1q
      

#ifdef B_DIFF
 character(len=MaxFNameLength)   :: timefile
 real*8, pointer                 :: time(:)
 real                            :: alpha,dt
 integer                         :: numit
#endif

 integer, pointer     :: seed(:)
 integer              :: iargc

 integer              :: istat

#ifdef DIAG_CROSSVAL
 character(len=200)   :: DirDiag
 character(len=200),                                                          &
      pointer         :: diagfnames(:)

#endif

 include 'ReadMatrix.h'


! interface
!   subroutine valsvd(U,S,V,IEX,JEX,VAL,M,N,maxnev,maxncv,ldu,ldv,NN)
!    integer,intent(in) :: maxnev,maxncv,ldu,ldv
!    real, intent(in) :: U(ldu, maxnev),V(ldv,maxncv),S(maxncv,2)
!    integer,intent(in) :: M,N,NN,IEX,JEX
!    real,intent(out) :: VAL
!   end subroutine valsvd
! end interface


 

 !      ----------------------------------------
 !      check if DINEOF has been called properly
 !      ----------------------------------------  

 if(iargc().ne.1) then
   write(stderr,*) 'DINEOF usage:'
   write(stderr,*) 'dineof <init filename> '
   stop
 end if

 call getarg(1,initfilename)

 !-------------------------------------------------------------------------
 !      Read numerical parameters from file initfilename
 !-------------------------------------------------------------------------
 call getInitValue(initfilename,'neini',neini) 
 call getInitValue(initfilename,'nev',nev) 
 call getInitValue(initfilename,'ncv',ncv) 
 call getInitValue(initfilename,'tol',tol) 
 call getInitValue(initfilename,'nitemax',nitemax)
 call getInitValue(initfilename,'toliter',toliter)
 call getInitValue(initfilename,'rec',rec) 
 call getInitValue(initfilename,'eof',eofno) 
 call getInitValue(initfilename,'norm',norma) 
 call getInitValue(initfilename,'Output',DirOutput) 
#ifdef DIAG_CROSSVAL
 call getInitValue(initfilename,'DirDiag',DirDiag) 
 call getInitValue(initfilename,'diags',diagfnames) 
#endif
 call getInitValue(initfilename,'data',sstfile) 

 if( presentInitValue(initfilename,'mask')) then
  call getInitValue(initfilename,'mask',maskfile)
 else 
  write(stdout,*)'no maskfile',size(maskfile)
 endif

#ifdef B_DIFF
 call getInitValue(initfilename,'alpha',alpha) 
 call getInitValue(initfilename,'numit',numit) 
 if(alpha.ne.0)  call getInitValue(initfilename,'time',timefile)
#endif
 call getInitValue(initfilename,'results',resultfnames) 
 call getInitValue(initfilename,'EOF.U',eofufnames) 
 call getInitValue(initfilename,'EOF.V',eofvfname) 
 call getInitValue(initfilename,'EOF.Sigma',eofsfname) 

 call getInitValue(initfilename,'seed',ipid) 


 if (stdout.ne.6) then
   call getInitValue(initfilename,'logfile',logfile,default='dineof.out') 
   open(stdout,file=logfile,status='unknown')
 end if

 write(stdout,*) '********************************************************************'
 write(stdout,*) 'Numerical data read'
 write(stdout,*) 'You entered the values:'
 write(stdout,*) 'number of EOF modes you want to compute',nev
 write(stdout,*) 'maximal size for the Krylov subspace',ncv
 write(stdout,*)

 call flush(stdout,istat)

 if(norma.eq.1) then
   write(stdout,*) 'You asked for the normalisation of the input matrices'
 else
   write(stdout,*) 'You asked not to normalise the input matrices'
 end if

 write(stdout,*)

 if(eofno.eq.1) then
   write(stdout,*) 'The right and left EOFs will be written in directory ',trim(DirOutput)
 else
   write(stdout,*) 'You asked not to write the left and right EOFs used for the reconstruction'
   write(stdout,*) '(to set this option ON, change eof=0 in file dineof.init)'
 end if
 write(stdout,*)

 write(stdout,*) '********************************************************************'
 !-----------------------------------
 call flush(stdout,istat)
 maxnev=nev
 maxncv=ncv
 !nitemax=35

 allocate(select(maxncv))

 call flush(stdout,istat)

 allocate(VALC(0:nev+1),VALC_all(0:nev+1),VALC_present(0:nev+1))

 nbvars = size(sstfile)
 write(stdout,*)
 do q=1,nbvars
   write(stdout,*)'You entered filenames ', trim(sstfile(q)) 
   if(presentInitValue(initfilename,'mask')) then
     write(stdout,*)'                      ', trim(maskfile(q)) 
   else
     write(stdout,*)'                      and no mask file: no land points in your domain'
   endif
 end do
 write(stdout,*)

 !--------ReadMatrix-------------------------------------------------------------------
 !     Reads input matrices and the exclusion value VALEX 

 !     Reshapes the input matrices into one matrix with size M*N (spatial and temporal)
 !     Input(sstfile(s), maskfile(s),norma)
 !     Output(X,M,N,VALEX,mask,imax,jmax,first,fileMean,fileStd)
 !-------------------------------------------------------------------------------------

 call ReadMatrix(sstfile,X, M, N,VALEX,maskfile,norma,imax,jmax,first,fileMean,fileStd)

 !--------------------------------------------
 write(stdout,*) '***'
 write(stdout,*) 'Matrix loaded ... Land points extracted...'
 write(stdout,*) 
 write(stdout,700) ' Size of the matrix used in DINEOF: ', m,' x ', n
 write(stdout,*) '***'
 !--------------------------------------------

 call flush(stdout,istat)
 

 !--------Read time if B_DIFF option is activated in ppdef.h
#ifdef B_DIFF
 
  if(alpha.ne.0) then
write(*,*)'reading time file',timefile
    call uload(trim(timefile),time, valextime) 
write(*,*)'time file read'
    dt = minval(time(2:N)-time(1:N-1))

   if(alpha > dt**2/2) then
     write(*,*)''
     write(*,*)'******** Attention! ************************************************'
     write(*,*)'********************************************************************'
     write(*,800)'Alpha = ',alpha
     write(*,800)'Alpha should not exceed ',dt**2/2
     write(*,*)''
     write(*,*)'Alpha value is too large: the filter is unstable, DINEOF will stop'  
     write(*,*)'********************************************************************' 
     write(*,*)'********************************************************************'  
     stop
   endif

   write(*,*)'Diffusion of covariance matrix activated'
   write(*,803)'with parameters: alpha = ',alpha,' and number of iterations = ',numit

   
  else
   allocate(time(N))
   do i=1,N   
    time(i) = i
   enddo 
  endif   
  
#endif
 !-------------------------------------------------------------------------------------

 call flush(stdout,istat)
 maxm = M
 maxn=N

 ldu = maxm
 ldv=maxn

 
 ! Transpose matrix if M le N----------------------------------------------
 istrans=.false. 
 if (M.lt.N) then

    allocate(X2(N,M))
    X2=transpose(X)
    deallocate(X)
    allocate(X(N,M))
    X=X2
    deallocate(X2)
  
    I=M 
    M=N 
    N=I 
    maxm = M
    maxn=N
    ldu = maxm
    ldv=maxn
    allocate(meanA(maxm),& 
      w(ldu), V(ldv,maxncv), U(ldu, maxnev),&
      workl(maxncv*(maxncv+8)), workd(3*maxn),&
      S(maxncv,2), resid(maxn))

    ISTRANS=.true. 
!   stop "PLEASE PROVIDE TRANSPOSED MATRIX SINCE M GT N" 
 else



    allocate(meanA(maxm),& 
      w(ldu), V(ldv,maxncv), U(ldu, maxnev),&
      workl(maxncv*(maxncv+8)), workd(3*maxn),&
      S(maxncv,2), resid(maxn))
 end if 

 !----------------------------------------------------------------
 !--------------------------------------------------------

 !     ----------------------------------------------------------
 !     First statistics on reshaped matrix
 !     Calculate mean, missing data, and cross validation points
 !     ---------------------------------------------------------- 

 XMEAN=0 
 IMEAN=0 

 do J=1,N 
   do I=1,M
     if(X(I,J).ne.VALEX) then 
       IMEAN=IMEAN+1 
       XMEAN=XMEAN+X(I,J) 
     end if
   end do
 end do

 XMEAN=XMEAN/IMEAN 
 IMISS=M*N-IMEAN 
 write(stdout,*)
 write(stdout,600) ' Missing data: ','',IMISS,' out of ', M*N,' (',IMISS*100./(M*N),'%)'
 write(stdout,*) 

 if (.not.presentInitValue(initfilename,'clouds')) then
   ICROSS=min(0.01*M*N+40,0.03*M*N)
 else
   ! user chooses the cross-validation points by the entry "clouds"
   call getInitValue(initfilename,'clouds',cloudsfname)
   call uload(cloudsfname,cvp_indexes,valexc)
   write(stdout,*) 'User-specified cross-validation points: ',trim(cloudsfname)
   icross = size(cvp_indexes,1)
 end if

 write(stdout,802) ' Number of cross validation points ',ICROSS 
 write(stdout,*) '********************************************************************'
 write(stdout,*)
 
 call flush(stdout,istat)
 
 IMISST=IMISS+ICROSS

 if(IMISST.gt.(M*N)) stop 'Not enough data' 

 !     ------------------------------------------------------
 !     Declare and initialise IEX(IMISST),JEX(IMISST)
 !     data actually mising from the input matrix
 !     ------------------------------------------------------ 
 
 allocate(XEX(IMISST,0:1),IEX(IMISST),JEX(IMISST), stat=istat)

 if (istat /= 0) then
   write(stderr,*) 'allocation error: XEX(IMISST,0:1),IEX(IMISST),JEX(IMISST)'
 endif

 call flush(stdout,istat)

 ! M*N OPS 
 IM=1
 do J=1,N 
   do I=1,M 
     if(X(I,J).eq.VALEX) then 
       IEX(IM)=I 
       JEX(IM)=J            
       IM=IM+1 

     else 
       !     -----------------
       !     substract mean 
       !     -----------------
       X(I,J)=X(I,J)-XMEAN 
     end if
   end do
 end do

 if(IM-1.ne.IMISS) stop '???' 

 !initial variance of the matrix X
 varini=sum(X**2,X.ne.valex)/count(X.ne.valex)

 stdvini=sqrt(varini)
 sumVarIni = M * N * varini

 !-------------------------------------------------------------------------------
 ! Add pseudo missing data points for cross validation, but save the values into 
 ! an array 
 ! declare XEX(IMISST,0:N) 
 ! XEX(II,NN) contains the interpolation at point IEX(II),JEX(II) when NN modes are used 
 ! For NN=O initial data are stored 



 if (.not.presentInitValue(initfilename,'clouds')) then
   ! choose cross-validation points internally

   !----random data set (random cross-validation points) ----------------------------------

   ! ipid=getpid() 
   call random_seed(SIZE = R)
   allocate(seed(R))
   seed = ipid
   call random_seed(put = seed(1:R))
   deallocate(seed)

   K=0 

   do while (k.lt.icross)

     call random_number(Y)
     I=floor(Y*M)+1
     call random_number(Y)  
     J=floor(Y*N)+1 



     if(X(I,J).ne.VALEX) then 
       ! for cross validation

       K=K+1 
       IEX(IMISS+K)=I 
       JEX(IMISS+K)=J 
       XEX(IMISS+K,0)=X(I,J)
       X(i,j) = valex 
     end if

   end do

 else
   ! user chooses the cross-validation points by the entry "clouds"
   do k=1,icross
     i = cvp_indexes(k,1)
     j = cvp_indexes(k,2)

     if(X(I,J).ne.VALEX) then 

       IEX(IMISS+K) = i
       JEX(IMISS+K) = j
       XEX(IMISS+K,0)=X(I,J)
       !write(stderr,*)'xex',X(i,j)
       ! for debuging
       !if (k.le.10) write(stdout,*) 'first cross validation point ',i,j,X(i,j)+XMEAN

       X(i,j) = valex
     else
       write(stderr,*) 'Error in ',trim(cloudsfname),i,j,k,X(I,J)
       stop
     end if
   end do

   deallocate(cvp_indexes)
 end if

 !initialise missing values (=VALEX) to zero
 do J=1,N 
   do I=1,M
     if(X(I,J).eq.VALEX) then 
       X(I,J)=0
     end if
   end do
 end do

!write(stderr,*) 'icross',k

 !-------------------------------------------------------------
 !    Estimates the total variance

 call svariExp&
      &     (X, M, N, maxm, maxn, sumSV, w)
 
!------------------------------------------------------------

 !-------------------------------------------------------------
 ! Now first guess of EOF 1 U,V,D can be calculated from XCR 
 !     -----------
 !-------------------------------------------------------------
 ! call LANCZOS 

 MAKENIT=1
 mytol=0

 call cpu_time(t1)
 IT=0
 call ssvd_lancz&
      &       (X, M, N, maxm, maxn, nev, ncv, maxnev, maxncv, ldv,& 
      &        ldu, w, V, U, workl, workd, S, resid, select, tol, IT,&
#ifdef B_DIFF
      &        time,alpha,numit,&
#endif
      &        MAKENIT,mytol,varini)

 call cpu_time(t2)

 write(stdout,800)' Time (in seconds) for 1 EOF mode calculation in DINEOF',t2-t1
 write(stdout,*)

 ! VALSVD: should give value at point I,J for P modes 
 ! Fill in missing data points with first guess

! do K=1,IMISST

   NN=1
   call valsvd(U,S,V,X,IEX,JEX,VAL,NN,IMISST)
!   call valsvd(U,S,V,IEX(K),JEX(K),VAL,M,N,maxnev,maxncv,ldu,ldv,NN) 
!   X(IEX(K),JEX(K))=VAL

! end do



 ! now loop on 'nev' modes asked for the reconstruction 
 !--------------------------------------------------------------------
 !--------------------------------------------------------------------
 ! 

 valc(0)=0
 valcopt=10000
 allocate(valosc(nev,nitemax),Xlastit(maxm,maxn),nitedone(nev),nitedonef(nev), stat=istat)

 if (istat /= 0) then
   write(stderr,*) 'allocation error: valosc(nev,nitemax),Xlastit(maxm,maxn),nitedone(nev),nitedonef(nev)'
 endif

 valosc(:,:)=-1
 nitedone(:)=0.0
 Xlastit(:,:)=0.0

 write(stdout,300) '# EOF modes asked: ',nev,'        Convergence level required: ',toliter
 write(stdout,*)
 write(stdout,*)'Starting with the EOF mode calculation...'
 write(stdout,*)
 write(stdout,*)'EOF mode    Expected Error    Iterations made   Convergence achieved'
 write(stdout,*)'________    ______________    _______________   ____________________'
 write(stdout,*)
 
 do P=neini,nev

   ! Iterations.

!   write(stdout,*) 'computing eigenvector',P,' (of',nev,')'
   call flush(stdout,istat)

   do NITE=1,nitemax

     !test no exploit of precedent iterations:it=0
     IT=0
     VAL=0.0
     call ssvd_lancz&
          &       (X, M, N, maxm, maxn,   P, ncv, maxnev, maxncv, ldv,& 
          &        ldu, w, V, U, workl, workd, S, resid, select, tol, IT,&
#ifdef B_DIFF
          &        time,alpha,numit,&
#endif
          &        MAKENIT,mytol,varini)

!     do K=1,IMISST 
       NN=P 
       call valsvd(U,S,V,X,IEX,JEX,VAL,NN,IMISST)
!       call valsvd(U,S,V,IEX(K),JEX(K),VAL,M,N,maxnev,maxncv,ldu,ldv,NN) 
!       X(IEX(K),JEX(K))=VAL 
!     end do
     
     VAL = 0

     do K=1,IMISST 
       DIFF=X(IEX(K),JEX(K))-Xlastit(IEX(K),JEX(K)) 
       VAL=VAL+DIFF*DIFF
       Xlastit(IEX(K),JEX(K))=X(IEX(K),JEX(K)) 
     end do   

     VAL=sqrt(VAL/FLOAT(IMISST))/stdvini      
 
     valosc(P,NITE)=VAL
     !take out this write out of this nite loop
     !call usave(trim(DirOutput)//'/valosc.dat',valosc,valex) 
     
     nitedone(P)=NITE

     ! add call flush here also?
     
     if (VAL < toliter) then
	exit
     end if

   end do

!   write(stdout,*) 'stopped after iteration ',nitedone(p),' , with valosc: ',valosc(p,nitedone(p)),' <? tolerance :',toliter
   
   nitedonef=FLOAT(nitedone)
!   call usave(trim(DirOutput)//'/nitedone.dat',nitedonef,valex)   

   ! Calculate Cross validator 
   VAL=0 

   do K=1,ICROSS 
     DIFF=X(IEX(IMISS+K),JEX(IMISS+K))-XEX(IMISS+K,0) 
     VAL=VAL+DIFF*DIFF 
   end do

   VALC(P)=sqrt(VAL/FLOAT(ICROSS)) 

   ! Store best guess with P EOFs only if valc(P) better than valc(P-1) ; add 
   if(VALC(P).lt.VALCopt) then 
     VAL=VALC(P)
     VALCopt=VALC(P) 
     Popt=P 
     do K=1,IMISST
     	XEX(K,1)=X(IEX(K),JEX(K)) 
     end do
   end if 
     
   !Output "expected error" for each number of modes selected, for diagnosis during the run 
!   write(stdout,*) '"expected error" with ',p,' modes selected: ',valc(p) 

   write(stdout,400)P,valc(p),nitedone(p),valosc(p,nitedone(p))

   call flush(stdout,istat)

!new cross-validation error

#ifdef ALL_POINTS_CROSS_VALIDATION
   val = 0

   do i=1,M
     do j=1,N
!change this?
       call valsvd(U,S,V,X,IEX,JEX,VAL,NN,IMISST)
!       call valsvd(U,S,V,i,j,Xrec,M,N,maxnev,maxncv,ldu,ldv,P) 
       DIFF=X(i,j)-Xrec
       VAL=VAL+DIFF*DIFF 
     end do
   end do

   VALC_all(P) = sqrt(VAL/(M*N))
   VALC_present(P) = sqrt((M*N*VALC_all(P)**2 - ICROSS * VALC(P)**2)/(M*N-ICROSS))

#endif

!end of new cross-validation error

   mytol=(1./100.)*(valc(p)-valc(p-1))/valc(p)

   !stop after three increasing CV values
   
   if ((P.gt.(neini+2))) then
     if ((VALC(P).gt.VALC(P-1)) .and.(VALC(P-1).gt.VALC(P-2)) .and. (VALC(P-2).gt.VALC(P-3))) then
      exit
     end if
   end if

 end do
 call cpu_time(t3)

 call usave(trim(DirOutput)//'/valc.dat',valc,valex)   
#ifdef ALL_POINTS_CROSS_VALIDATION
 call usave(trim(DirOutput)//'/valc_all.dat',valc_all,valex)   
 call usave(trim(DirOutput)//'/valc_present.dat',valc_present,valex)   
#endif

 !-----------------------------------------------------------------
 ! Finished all values for P 
 !-----------------------------------------------------------------

 ! Find value of P which corresponds to the minimum value of VALC 
 P=Popt   ! already detected in loop on P modes
 

#ifdef DIAG_CROSSVAL				
!! not yet adapted to autom. stop criteria
 ! TO DO!!!
 !transpose is not implemented in this ifdef section!

 ! Recalculate the EOFs for optimal number of modes (P) without the cross-validation points 
 ! for diagnostic purposes

 allocate(Xp(M,N), stat=istat)

 if (istat /= 0) then
   write(stderr,*) 'allocation error: Xp(M,N)'
 endif

 do NITE=1,nitemax

   call ssvd_lancz&
        &       (X, M, N, maxm, maxn,   P, ncv, maxnev, maxncv, ldv,& 
        &        ldu, w, V, U, workl, workd, S, resid, select, tol, IT,&
#ifdef B_DIFF
        &        time,alpha,numit,&
#endif
        &        MAKENIT,mytol,varini)

!   do K=1,IMISS 
     NN=P 
     call valsvd(U,S,V,X,IEX,JEX,VAL,NN,IMISST)
!     call valsvd(U,S,V,IEX(K),JEX(K),VAL,M,N,maxnev,maxncv,ldu,ldv,NN) 
!     X(IEX(K),JEX(K))=VAL 
!   end do

 end do

   ! ReCalculate Cross validator 
   VAL=0 

   do K=1,ICROSS 
     DIFF=X(IEX(IMISS+K),JEX(IMISS+K))-XEX(IMISS+K,0) 
     VAL=VAL+DIFF*DIFF 
   end do

   VAL=sqrt(VAL/FLOAT(ICROSS)) 

   write(stdout,*) 'check cross-validator ',val


 if(rec.eq.1) then
   do i=1,M
     do j=1,N
       Xp(i,j)=Xmean
       do k=1,P
         Xp(i,j)=Xp(i,j)+U(i,k)*S(k,1)*V(j,k)
       end do
     end do
   end do
 end if

 call writeSVD(DirDiag,s,nev,u,n,m,v,sumVarIni,SumVarEnd)

 call writeMatrix(Xp,xmean,valex,diagfnames,maskfile,norma,M,N,imax,jmax,first,DirDiag,valc,fileMean,fileStd)

#endif

 !--------------------------------------------------------------

 ! Now optimal number P of EOFs is known: 
 ! Retrieve best guess from that number 

 do K=1,IMISS 
   X(IEX(K),JEX(K))=XEX(K,1) 
 end do
!save best reconstruction estimate for CV befor re-including them in the matrix X
call usave(trim(DirOutput)//'/CVpoints_best_estimate.dat',XEX(IMISS+1:IMISST,1),valex)
call usave(trim(DirOutput)//'/CVpoints_initial.dat',XEX(IMISS+1:IMISST,0),valex)

 !--------------------------------------------------------------
 ! Use now the points put aside from initial data set 

 do K=1,ICROSS 
   X(IEX(IMISS+K),JEX(IMISS+K))=XEX(IMISS+K,0) 
 end do

 allocate(valosclast(nitemax))
 valosclast(:)=-1
 nitedonelast=0
 Xlastit(:,:)=0.0

 IT=1
 write(stdout,*)
 write(stdout,*) ' Minimum reached in cross-validation'
 write(stdout,*) ' Number of optimal EOF modes: ',P
 write(stdout,*) 
 write(stdout,*) ' Make last reconstruction, including data put aside for cross-validation'
 NITE=1
 
 do while (NITE < nitemax+1) 

   VAL=0.0
   call ssvd_lancz&
        &       (X, M, N, maxm, maxn,   P, ncv, maxnev, maxncv, ldv,& 
        &        ldu, w, V, U, workl, workd, S, resid, select, tol, IT,&
#ifdef B_DIFF
        &        time,alpha,numit,&
#endif
        &        MAKENIT,mytol,varini)

!   do K=1,IMISS 
     NN=P 
     call valsvd(U,S,V,X,IEX,JEX,VAL,NN,IMISST)
!     call valsvd(U,S,V,IEX(K),JEX(K),VAL,M,N,maxnev,maxncv,ldu,ldv,NN) 
!     X(IEX(K),JEX(K))=VAL 
!   end do

   VAL = 0

   do K=1,IMISS 
       DIFF=X(IEX(K),JEX(K))-Xlastit(IEX(K),JEX(K)) 
       VAL=VAL+DIFF*DIFF
       Xlastit(IEX(K),JEX(K))=X(IEX(K),JEX(K)) 
   end do   

   VAL=sqrt(VAL/FLOAT(IMISS))/stdvini      
 
   valosclast(nite)=VAL
   call usave(trim(DirOutput)//'/valosclast.dat',valosclast,valex) 
     
   nitedonelast=NITE

   ! add call flush here also?

   NITE=NITE+1
     
   if (VAL < toliter) then
	NITE=nitemax+1
   end if

 end do
 
 write(stdout,*)
 write(stdout,400)P,valc(p),nitedone(p),valosc(p,nitedone(p))
! write(stdout,*) 'stopped after iteration ',nitedonelast,' , with valosc: ',valosclast(nitedonelast),' <? tolerance :',toliter
 
 SumVarEnd =  M * N * sum(X**2,X.ne.valex)/count(X.ne.valex)

 ! Reconstruct the matrix with the EOFs retained as optimal
 if (rec.eq.1) then
   do i=1,M
     do j=1,N
       X(i,j)=0
       do k=1,P
         X(i,j)=X(i,j)+U(i,k)*S(k,1)*V(j,k)
       end do
     end do
   end do
 end if



 write(stdout,*)
 write(stdout,*)'DINEOF finished!'
 write(stdout,*)
 write(stdout,801)'number of eigenvalues retained for the reconstruction        ',P
 write(stdout,800)'expected error calculated by cross-validation ',valc(P)
 write(stdout,800)'total time (in seconds) in lanczos process ', t3-t1
 write(stdout,*)
 write(stdout,*)'Now writing data...'
 write(stdout,*)
 write(stdout,*)
 write(stdout,*) 'total variance of the initial matrix ',sumVarIni/(M*N)
 write(stdout,*) 'total variance of the reconstructed matrix ',sumVarEnd/(M*N)

 if(eofno.eq.1) then
   

      do q=1,nbvars
        name1 = eofufnames(q)
        i = index(name1,'#',.true.)
        name1q = name1(1:i-1)
        call clear_files(name1q)      
      enddo

      i = index(eofvfname,'#',.true.)
      name2 = eofvfname(1:i-1)
      call clear_files(name2)
      
      i = index(eofsfname,'#',.true.)
      name3 = eofsfname(1:i-1)
      call clear_files(name3)


      if (istrans) then

          ! u and v have to be swaped if matrix has been transposed
          call writeSVD(DirOutput,eofvfname,eofsfname,s(1:P,1),u(:,1:P),valc(1:),sumVarEnd)
          call writeMatrix(v(:,1:P),0.,valex,eofufnames,maskfile,0,imax,jmax,first,DirOutput)
      else
          call writeSVD(DirOutput,eofvfname,eofsfname,s(1:P,1),v(:,1:P),valc(1:),sumVarEnd)
          call writeMatrix(u(:,1:P),0.,valex,eofufnames,maskfile,0,imax,jmax,first,DirOutput)
     endif
endif


if (istrans) then

 allocate(X2(N,M))
 X2=transpose(X)
 deallocate(X)
 allocate(X(N,M))
 X=X2
 deallocate(X2)
    I=M 
    M=N 
    N=I 
    maxm = M
    maxn=N
endif

 where (X.ne.valex) X=X+xmean 



 call writeMatrix(X,xmean,valex,resultfnames,maskfile,norma,imax,jmax,first,DirOutput,fileMean,fileStd)

 if(associated(maskfile)) then
   deallocate(maskfile)
 endif

 deallocate(X,imax,jmax,first,fileMean,fileStd,sstfile,resultfnames)

#ifdef B_DIFF
 deallocate(time)
#endif

 open(1,file=trim(DirOutput)//'/neofretained.val')
 write(1,500) P
 close(1)
 
 open(1,file=trim(DirOutput)//'/meandata.val')
 write(1,501) xmean
 close(1)

500 format(i3)
501 format(e10.4)
300 format(a20,i3,a35,e8.1)
400 format(i7,f20.4,i19,e23.4)
600 format(a14,a15,i11,a8,i11,a2,f6.2,a2)
700 format(a35,i25,a3,i6)
800 format(a55,f15.4)
801 format(a54,i3)
802 format(a34,i35)
803 format(a26,f15.2,a31,i3)

 write(stdout,*)
 write(stdout,*)'...done!'
 write(stdout,*)

 if (stdout.ne.6) close(stdout)

end program dineof






