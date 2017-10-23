#include "ppdef.h"      

      subroutine ssvd_lancz
     &   (A, m, n, maxm, maxn, nev, ncv, maxnev, maxncv, ldv, ldu,
     &    w, v, u, workl, workd, s, resid, select, tol,IT,
#ifdef B_DIFF
     &    time,alpha,numit,
#endif
     &    MAKENIT,mytol,varini)
       use ufileformat
       use dineof_utils      
       implicit none
c---------------------------------------------------------------------c
c     "iterativeEOF", Version 2.0, December, 2002.                    c
c      A Krylov-based code for computing EOF                          c
c      Vincent TOUMAZOU                                               c
c      Vincent.Toumazou@cnes.fr                                       c
c      CNES/MERCATOR-Ocean                                            c
c      18, Avenue Edouard Belin                                       c
c      F-31401 TOULOUSE Cedex 4                                       c
c      URL: http://www.mercator.com.fr                                c
c      URL: http://www.cnes.fr                                        c
c---------------------------------------------------------------------c
c     This subroutine computes NEV singular values for statistical    c
c     study of phenomenom. Left and right associated singular vectors c
c     are also computed.                                              c
c                                                                     c
c********************** CAUTION **************************************c
c                                                                     c
c     The computation of the singular elements is performed using     c
c     the freeware ARPACK. Please respect the copyright of ARPACK     c
c     and read the documentation(s) provided with the EOF code before c
c     any use.                                                        c
c                                                                     c
c********************** END OF CAUTION *******************************c
c                                                                     c
c     The forthcoming lines of code are an adaptation of the "Example c
c     Test Program" called "dsvd" from the ARPACK Freeware.           c
c     In order to tune this program (and to transform it into a       c
c     subroutine) some modifications have been done both in the       c
c     comments and in the statements and/or declarations.             c
c                                                                     c
c---------------------------------------------------------------------c
c
c     This subroutine calls ARPACK to find a few of the
c     largest singular values(sigma) and corresponding right singular
c     vectors (v) for the the matrix A by solving the symmetric problem:
c
c                        (A'*A)*v = sigma*v
c 
c     where A is an m by n real matrix.
c
c     This formulation is appropriate when  m  .ge.  n.
c     Reverse the roles of A and A' in the case that  m .le. n.
c
c     The only thing that must be supplied in order to use this
c     routine on your problem is to change the array dimensions 
c     appropriately, to specify WHICH singular values you want to 
c     compute and to supply a the matrix-vector products 
c
c                         w <-  Ax
c                         y <-  A'w
c
c     Further documentation is available in the header of DSAUPD
c     which may be found in the SRC directory.
c
c Routines called:
c ----------------
c     ssaupd  ARPACK reverse communication interface routine.
c     sseupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     snrm2   Level 1 BLAS that computes the norm of a vector.
c     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     sscal   Level 1 BLAS thst computes x <- x*alpha.
c     scopy   Level 1 BLAS thst computes y <- x.
c
c Authors of the ARPACK kernel:
c -----------------------------
c     Richard Lehoucq
c     Danny Sorensen
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c     For more informations about ARPACK (Use, Updates, Copyright), please
c     read first the documentation provided with the EOF code.
c
c-----------------------------------------------------------------------
c
c     %------------------------------------------------------%
c     | Storage Declarations:                                |
c     |                                                      |
c     | It is assumed that A is M by N with M .ge. N.        |
c     |                                                      |
c     | The maximum dimensions for all arrays are            |
c     | set here to accommodate a problem size of            |
c     | M .le. MAXM  and  N .le. MAXN                        |
c     |                                                      |
c     | The NEV right singular vectors will be computed in   |
c     | the N by NCV array V.                                |
c     |                                                      |
c     | The NEV left singular vectors will be computed in    |
c     | the M by NEV array U.                                |
c     |                                                      |
c     | NEV is the number of singular values requested.      |
c     |     See specifications for ARPACK usage below.       |
c     |                                                      |
c     | NCV is the largest number of basis vectors that will |
c     |     be used in the Implicitly Restarted Arnoldi      |
c     |     Process.  Work per major iteration is            |
c     |     proportional to N*NCV*NCV.                       |
c     |                                                      |
c     | You must set:                                        |
c     |                                                      |
c     | MAXM:   Maximum number of rows of the A allowed.     |
c     | MAXN:   Maximum number of columns of the A allowed.  |
c     | MAXNEV: Maximum NEV allowed                          |
c     | MAXNCV: Maximum NCV allowed                          |
c     %------------------------------------------------------%
c
      integer      maxm, maxn, maxnev, maxncv, ldv, ldu,IT,MAKENIT
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      Real
     &             A(maxm,maxn), w(ldu), v(ldv,maxncv),u(ldu,maxnev), 
     &             workl(maxncv*(maxncv+8)), workd(3*maxn), 
     &             s(maxncv,2), resid(maxn),frigo,x(maxn),
     &             varA,varini,weighB,w1,w2
#ifdef B_DIFF
      Real        :: time(n),alpha,AF(maxm,maxn),A2(maxn)
      integer     :: numit
#endif
      real,pointer     :: ticons(:)
      Real,allocatable :: B(:,:),Bp(:,:),Bp2(:),new_w(:,:)
      logical      select(maxncv)
      integer      iparam(11), ipntr(11)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character    bmat*1, which*2
      character(len=200)   :: dirOutput
      integer      ido, m, n, nev, ncv, lworkl, info, ierr,
     &             j, ishfts, maxitr, mode1, nconv, i,T_ticons,
     &             k
      logical      rvec
      Real      
     &             tol, sigma, temp, mytol, newtol,valex
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &             one, zero
      parameter    (one = 1.0E+0, zero = 0.0E+0)
c  
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Real           
     &             snrm2
      external     snrm2, saxpy, scopy, sscal, sgemv
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------------------------%
c     | The following include statement and assignments |
c     | initiate trace output from the internal         |
c     | actions of ARPACK.  See debug.doc in the        |
c     | DOCUMENTS directory for usage.  Initially, the  |
c     | most useful information will be a breakdown of  |
c     | time spent in the various stages of computation |
c     | given by setting msaupd = 1.                    |
c     %-------------------------------------------------%
c
      include 'includes/debug.h'
      ndigit = -3
      logfil = 6
      msgets = 0
      msaitr = 0 
      msapps = 0
      msaupd = 0
      msaup2 = 0
      mseigt = 0
      mseupd = 0
c
      bmat  = 'I'
      which = 'LA'
c
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SVD: N is greater than MAXN '
         go to 9000
      else if ( m .gt. maxm ) then
         print *, ' ERROR with _SVD: M is greater than MAXM '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SVD: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SVD: NCV is greater than MAXNCV '
         go to 9000
      end if

c Calculate matrix B = A'A
       
       allocate(B(n,n))
      
       B = matmul(transpose(A),A)
c       call usave('B.dat',B,valex);
       
       
#ifdef B_TICONS
       call uload('time_constraint.dat',ticons,valex)
       T_ticons = 30
c       varA = sum(A**2)/(size(A,1)*size(A,2))
       alpha = m * varini
c       alpha = 2e6
       do i = 1,size(B,1)
          do j = 1,size(B,2)
       B(i,j) = B(i,j) +
     &       alpha * exp(-(ticons(i)-ticons(j))**2/T_ticons**2)
          end do
       end do
#endif
#ifdef B_FILTER
       allocate(Bp(n,n),new_w(n,n))
       Bp = B
       T_ticons = 7
       call uload('time_constraint.dat',ticons,valex)
       weighB = 0.001
       new_w = 0.001
c       do i = 1,10
c          B(2:n-1,2:n-1) = 
c     &                 (1-2*weighB) *B(2:n-1,2:n-1) + 
c     &                      weighB  *B(1:n-2,1:n-2) + 
c     &                      weighB  *B(3:n,3:n)
c       end do

c explicit loops
c       do i = 1,10
c          do j = 2,size(B,1)-1
c             do k = 2,size(B,1)-1
c              Bp(j,k) = (1-2*weighB)*B(j,k) + 
c     &                   weighB*B(j-1,k-1)  + 
c     &                   weighB*B(j+1,k+1)
c             end do
c          end do
c          B = Bp
c       end do

       do j = 2,size(B,1)-1
          do k = 2,size(B,1)-1
             new_w(j,k) = 0.02 * 
     &                    exp(-(ticons(j)-ticons(k))**2/T_ticons**2)
          end do
       end do
c       call usave('w.dat',new_w,valex);
 
       do i = 1,10
          do j = 2,size(B,1)-1
             do k = 2,size(B,1)-1
              w1 = new_w(j,j-1)*new_w(k,k-1)
              w2 = new_w(j,j+1)*new_w(k,k+1)
              Bp(j,k) = (1-w1-w2)*B(j,k)      + 
     &                        w1 *B(j-1,k-1)  + 
     &                        w2 *B(j+1,k+1)
             end do
          end do       
          B = Bp
       end do

       
#endif

#ifdef B_DIFF
        allocate(Bp2(n),Bp(n,n))
        
c        alpha = 0.01
!        alpha = 0.001
c        numit = 10
c        call uload('time_constraint.dat',ticons,valex)
        do i=1,n
        Bp2 = B(i,:)        
        call dindiff(time,Bp2,alpha,numit)
        Bp(i,:)=Bp2
        end do
        B = transpose(Bp)
        do i=1,n
        Bp2 = B(i,:)        
        call dindiff(time,Bp2,alpha,numit)
        Bp(i,:)=Bp2
        end do
        B = transpose(Bp)
        deallocate(Bp2,Bp)
c        do i=1,m
c        A2 = A(i,:)
c        call dindiff(time,A2,alpha,numit)
c        AF(i,:)=A2
c        end do
c        A=AF

#endif
c       call usave('B.dat',B,valex);
c       call usave('AF.dat',AF,valex);
       
c       

c
c     %-----------------------------------------------------%
c     | Specification of stopping rules and initial         |
c     | conditions before calling SSAUPD                    |
c     |                                                     |
c     |           abs(sigmaC - sigmaT) < TOL*abs(sigmaC)    |
c     |               computed   true                       |
c     |                                                     |
c     |      If TOL .le. 0,  then TOL <- macheps            |
c     |              (machine precision) is used.           |
c     |                                                     |
c     | IDO  is the REVERSE COMMUNICATION parameter         |
c     |      used to specify actions to be taken on return  |
c     |      from SSAUPD. (See usage below.)                |
c     |                                                     |
c     |      It MUST initially be set to 0 before the first |
c     |      call to SSAUPD.                                | 
c     |                                                     |
c     | INFO on entry specifies starting vector information |
c     |      and on return indicates error codes            |
c     |                                                     |
c     |      Initially, setting INFO=0 indicates that a     | 
c     |      random starting vector is requested to         |
c     |      start the ARNOLDI iteration.  Setting INFO to  |
c     |      a nonzero value on the initial call is used    |
c     |      if you want to specify your own starting       |
c     |      vector (This vector must be placed in RESID.)  | 
c     |                                                     |
c     | The work array WORKL is used in SSAUPD as           | 
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.                                  |
c     %-----------------------------------------------------%
c
      lworkl = ncv*(ncv+8)
      info = 0
      ido = 0
c
c     %---------------------------------------------------%
c     | Specification of Algorithm Mode:                  |
c     |                                                   |
c     | This program uses the exact shift strategy        |
c     | (indicated by setting IPARAM(1) = 1.)             |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of SSAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | SSAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = n

c$$$      if(MAKENIT.eq.0) then
c$$$         if(nev**2.le.n) then
c$$$            maxitr = nev**2
c$$$         else 
c$$$            maxitr=n
c$$$         end if
c$$$      elseif(MAKENIT.eq.1) then
c$$$         maxitr = n
c$$$      end if

      mode1 = 1
c
      iparam(1) = ishfts
c                
      iparam(3) = maxitr
c                  
      iparam(7) = mode1
c
c     %------------------------------------------------%
c     | M A I N   L O O P (Reverse communication loop) |
c     %------------------------------------------------%
c

c  nev era 1
       if(IT.EQ.1) then
        do i=1,ldv
         resid(i)=v(i,1)
c          resid(i)=0
c         write(*,*)'aqui',v(i,1)
        end do
        info=1
        
       end if
       


 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine SSAUPD and take | 
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
      newtol=max(tol,mytol)
         call ssaupd ( ido, bmat, n, which, nev, newtol, resid, 
     &                 ncv, v, ldv, iparam, ipntr, workd, workl,
     &                 lworkl, info )
c
      
         if(IT.EQ.1) then
           !IT=0 ?????
           IT=0
           info=0
         end if

         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %---------------------------------------%
c           | Perform matrix vector multiplications |
c           |              w <--- A*x       (av())  |
c           |              y <--- A'*w      (atv()) |
c           | The user should supply his/her own    |
c           | matrix vector multiplication routines |
c           | here that takes workd(ipntr(1)) as    |
c           | the input, and returns the result in  |
c           | workd(ipntr(2)).                      |
c           %---------------------------------------%
c

          x=workd(ipntr(1):ipntr(1)+n-1)

          if (MAKENIT .eq. 0) then
            call sgemv('N',m,n,one,A,maxm,workd(ipntr(1)),1,zero,
     &                w,1)
            call sgemv('T',m,n,one,A,maxm,w,1,zero,
     &                workd(ipntr(2)),1)
 
          else

            workd(ipntr(2):ipntr(2)+n-1) = 
     &         matmul(B,x)

c            workd(ipntr(2):ipntr(2)+n-1) =    
c     &       matmul(matmul(transpose(A),A),workd(ipntr(1):ipntr(1)+n-1))
          end if     
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call SSAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         end if 
c
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message. Check the |
c        | documentation in SSAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
c
      else 
c
c        %--------------------------------------------%
c        | No fatal errors occurred.                  |
c        | Post-Process using SSEUPD.                 |
c        |                                            |
c        | Computed singular values may be extracted. |  
c        |                                            |
c        | Singular vectors may also be computed now  |
c        | if desired.  (indicated by rvec = .true.)  | 
c        |                                            |
c        | The routine SSEUPD now called to do this   |
c        | post processing                            | 
c        %--------------------------------------------%
c           
         rvec = .true.
c
         call sseupd ( rvec, 'All', select, s, v, ldv, sigma, 
     &        bmat, n, which, nev, newtol, resid, ncv, v, ldv, 
     &        iparam, ipntr, workd, workl, lworkl, ierr )
        
c
c        %-----------------------------------------------%
c        | Singular values are returned in the first     |
c        | column of the two dimensional array S         |
c        | and the corresponding right singular vectors  | 
c        | are returned in the first NEV columns of the  |
c        | two dimensional array V as requested here.    |
c        %-----------------------------------------------%
c
         if ( ierr .ne. 0) then
c
c           %------------------------------------%
c           | Error condition:                   |
c           | Check the documentation of SSEUPD. |
c           %------------------------------------%
c
            print *, ' '
            print *, ' Error with _seupd, info = ', ierr
            print *, ' Check the documentation of _seupd. '
            print *, ' '
c
         else
c
            nconv =  iparam(5)
c            write(*,*)'eof#,sum(lambda)/m/n',nev,sum(s(1:nev,1))/m/n
c            write(*,*)'s factor',varini*m*n/sum(s(1:nev,1))
 
c            write(*,*)'var A',varA,varini
            do 20 j=1, nconv
c
#ifdef B_TICONS
               s(j,1) = s(j,1)*varini*m*n/sum(s(1:nev,1))
#endif
               s(j,1) = sqrt(s(j,1))
c
c              %-----------------------------%
c              | Compute the left singular   |
c              | vectors from the formula    |
c              |                             |
c              |     u = Av/sigma            |
c              |                             |
c              | u should have norm 1 so     |
c              | divide by norm(Av) instead. |
c              %-----------------------------%
c
               call sgemv('N',m,n,one,A,maxm,v(1,j),1,
     &                   zero,u(1,j),1)
               temp = one/s(j,1)
               call sscal(m, temp, u(1,j), 1)
               temp = one/snrm2(m, u(1,j), 1)
               call sscal(m, temp, u(1,j), 1)
               
c
c              %---------------------------%
c              |                           |
c              | Compute the residual norm |
c              |                           |
c              |   ||  A*v - sigma*u ||    |
c              |                           |
c              | for the NCONV accurately  |
c              | computed singular values  |
c              | and vectors.  (iparam(5)  |
c              | indicates how many are    |
c              | accurate to the requested |
c              | tolerance).               |
c              | Store the result in 2nd   |
c              | column of array S.        |
c              %---------------------------%
c
               call scopy(m, u(1,j), 1, w, 1)
               call sgemv('N',m,n,one,A,maxm,v(1,j),1,-s(j,1),w,1)
               s(j,2) = snrm2(m, w, 1)
c
 20         continue
c
c           %-------------------------------%
c           | Display computed residuals    |
c           %-------------------------------%
c
            if (msaupd.gt.0) then
            call smout(6, nconv, 2, s, maxncv, -6,
     &                'Singular values and direct residuals')
            end if
         end if
c
c        %------------------------------------------%
c        | Print additional convergence information |
c        %------------------------------------------%
c
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' ' 
            print *, ' No shifts could be applied during implicit
     &                 Arnoldi update, try increasing NCV.'
            print *, ' '
         end if      

         if (msaupd.gt.0) then 
         print *, ' '
         print *, ' _SVD '
         print *, ' ==== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated',
     &            ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', 
     &              nconv 
         print *, ' The number of Implicit Arnoldi update',
     &            ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol, newtol
         print *, ' '
         end if
           

      end if
c
c     %-------------------------%
c     | Done with program ssvd. |
c     %-------------------------%
c
         deallocate(B)       
 9000 continue
c v(ldv,maxncv), u(ldu, maxnev) s(maxncv,2)

          do i=1,ldu
           do j=1,nev/2
            frigo=u(i,j)
            u(i,j)=u(i,nev+1-j)
            u(i,nev+1-j)=frigo
           enddo
          enddo

          do i=1,ldv
           do j=1,nev/2
            frigo=v(i,j)
            v(i,j)=v(i,nev+1-j)
            v(i,nev+1-j)=frigo
           enddo
          enddo

          do i=1,nev/2     
           do j = 1,2
            frigo=s(i,j)      
            s(i,j)=s(nev+1-i,j)
            s(nev+1-i,j)=frigo
           enddo
          enddo     
          
      end
