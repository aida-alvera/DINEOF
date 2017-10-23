      program siterativeEof
      use ufileformat
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

c$$$       interface
c$$$       subroutine readMatrix(fname,matrix,imax,jmax)
c$$$       character(len=*),intent(in) :: fname
c$$$       real,intent(out) :: matrix(:,:)
c$$$       integer,intent(out) :: imax,jmax
c$$$       end interface

      integer          maxm, maxn, maxnev, maxncv, ldv, ldu
      parameter       (maxm = 40000, maxn=250, maxnev=20, maxncv=25,
     &                 ldu = maxm, ldv=maxn, nlines=22)
      Real
     &                 meanA(maxm), 
     &                 w(ldu), v(ldv,maxncv), u(ldu, maxnev),
     &                 workl(maxncv*(maxncv+8)), workd(3*maxn),
     &                 s(maxncv,2), resid(maxn), 
     &                 tol, t, sumSV

      real, pointer :: A(:,:)
      logical          select(maxncv)
      integer          M, N, nev, ncv
      character*72     comment
      character(len=50) list(22)
     
c 
c

      write(*,*) '********************************************'
      write(*,*) '**********  starting program  **************'
      write(*,*) '********************************************'

c
c     Numerical parameters
c
      open(1,file='Input/numerical.data')
      read(1,500) comment
      read(1,500) comment
      read(1,500) comment
      read(1,500) comment
      read(1,500) comment
      read(1,*) nev
      read(1,500) comment
      read(1,500) comment
      read(1,*) ncv
      read(1,500) comment
      read(1,500) comment
      read(1,*) tol
      close(1)
      write(*,*) '***'
      write(*,*) 'Numerical data read'
      write(*,*) 'Number of required modes         = ', nev
      write(*,*) 'Size of projection               = ', ncv 
      write(*,*) 'Treshold for Lanczos convergence = ', tol
      write(*,*) '***'
c
c   Reads the matrix of data
c

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CAUTION : Call your subroutine for reading     C
C            the matrix of data.                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     call readMatrix
C    &     (Add your parameters here, A, m and n
C    &      must be the outputs of your subroutine
C    &      and maxm, maxn the leading dimension of A
C    &      have to be in the inputs of the subroutine)
C Example : call readMatrix(A, maxm, maxn, m, n)

      
      open(77,file = '/u/alvera/Datafiles/Analist/analistTEM.dat')
      
      do i=1,nlines
      read(77,'(a50)') list(i)
      end do
     
      call readMatrix(list,A, M, N)

      write(*,*) '***'
      write(*,*) 'Matrix loaded'
      write(*,*) 'Size of the matrix = ', m,' by ', n
      write(*,*) '***'


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CAUTION : Uncomment the 2 next calls if you    C
C            want to center your data with the    C
C            mean on each point with respect      C
C            to the total period.                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c    Computes the Row-means of the matrix
c
      call smeanByRow
     &     (A, m, n, meanA, maxm, maxn)

c
c    Transforms the matrix such that its raws have a zero mean
c
      call smeanToZero
     &     (A, m, n, meanA, maxm, maxn)

c
c    Estimates the total variance
c
      call svariExp
     &     (A, m, n, maxm, maxn, sumSV, w)

c
c   Computes the NEV singular values and the modes of interest
c

      write(*,*) '********************************************'
      write(*,*) '*****  Starting the Lanczos Process  *******'
      write(*,*) '********************************************'
      
      call ssvd_lancz
     &       (A, m, n, maxm, maxn, nev, ncv, maxnev, maxncv, ldv, 
     &        ldu, w, v, u, workl, workd, s, resid, select, tol)
     
      write(*,*) '********************************************'
      write(*,*) '*******  End of the Lanczos Process  *******'
      write(*,*) '********************************************'

c
c Uncomment the next call if you want to normalize your outputs
c
c     call snormal_output
c    &        (s, u, v, m, n, nev, maxncv, ldu, ldv,
c    &         maxnev, maxm)

c
c   Extracts singular values and associated modes, prepares the 
c   output file for graphical exploitation.
c

c
c   Singular values
c   _______________
      open(1,file='Output/outputEof.vlsng')
      do i=nev,1,-1
         write(1,501) s(i,1)
      enddo
      close(1)
c
c   Variance detailed
c   _______________
      open(1,file='Output/outputEof.varEx')
      do i=nev,1,-1
         write(1,504) 'Mode ', nev+1-i, ' = ',
     &               100.0*s(i,1)*s(i,1)/sumSV, ' %'
      enddo
      close(1)

c
c   Spatial modes 
c   _______________
      open(1,file='Output/outputEof.lftvec')
      do j=1,m
         write(1,503) (u(j,i),i=nev,1,-1)
      enddo
      close(1)
c
c   Temporal modes 
c   _______________
      open(1,file='Output/outputEof.rghvec')
      do j=1,n
         write(1,503) (v(j,i),i=nev,1,-1)
      enddo
      close(1)
      close(77)
 500  format(a72)
 501  format(e10.4)
 502  format(3e10.4)
 503  format(300(e10.4,1X))
 504  format(a5,i3,a3,f8.4,a2)
      end





