#include "ppdef.h"

      subroutine svariExp
     &     (A, m, n, maxm, maxn, sumSV, w)

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
c      This subroutine computes the total variance by computing       c
c      trace(A^T.A).                                                  c
c                                                                     c
c      INPUT  :                                                       c
c      ------                                                         c
c         - A        : SINGLE PRECISION                               c
c                      m by n matrix                                  c
c         - m, n     : INTEGER                                        c
c                      Dimensions of A                                c
c         - maxm, maxn are the leading dimension of A                 c
c         - w        : SINGLE PRECISION                               c
c                      Vector of work.
c                                                                     c
c      OUTPUT :                                                       c
c      ------                                                         c
c         - sumSV    : SINGLE PRECISION                               c
c                      Sum of the square of the Singular Values of A  c
c---------------------------------------------------------------------c

      integer          m, n, maxm, maxn
      real             A(maxm, maxn), sumSV, w(*)

      integer          i,j
      real             tmp1, tmp2, snrm2
      external         snrm2

      tmp1=0.0
      tmp2=0.0
      
      do i=1,n
         do j=1,m
            w(j)=A(j,i)
         enddo
         tmp2=snrm2(m,w,1)
         tmp1=tmp2*tmp2+tmp1
      enddo
      sumSV=tmp1

      end
