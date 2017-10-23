      subroutine smeanByRow
     &     (A, m, n, meanA, maxm, maxn,VALEX)

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
c      This subroutine computes the raw mean of the matrix A          c
c                                                                     c
c      INPUT  :                                                       c
c      ------                                                         c
c         - A        : SINGLE PRECISION                               c
c                      m by n matrix                                  c
c         - m, n     : INTEGER                                        c
c                      Dimensions of A                                c
c         - maxm, maxn are the leading dimension of A                 c
c                                                                     c
c      OUTPUT :                                                       c
c      ------                                                         c
c         - meanA    : SINGLE PRECISION                               c
c                      vector of length m whose components are the    c
c                      means of the raws of A.                        c
c---------------------------------------------------------------------c

      integer          m, n, maxm, maxn
      real             A(maxm, maxn), meanA(maxm),VALEX

      integer          i,j,n1
      n1=0
      do i=1,m
         meanA(i)=0.0
         do j=1,n
            if(A(i,j).ne.VALEX) then
            meanA(i)=meanA(i)+A(i,j)
            n1=n1+1
            endif
         enddo
         meanA(i)=meanA(i)/n
      enddo
      end
