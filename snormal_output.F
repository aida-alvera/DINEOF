      subroutine snormal_output
     &        (s, u, v, m, n, nev, maxncv, ldu, ldv,
     &         maxnev, maxm)
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
c      This subroutine normalizes the temporal and spatial modes.     c
c      Each temporal mode is scaled by its component of largest       c
c      modulus, therefore max(|v|)=1. The spatial mode is divided     c
c      by the factor used to scale v and multiplied by the associated c
c      singular value.                                                c
c                                                                     c
c      INPUT  :                                                       c
c      ------                                                         c
c         - s     : DOUBLE PRECISION                                  c
c                   Array containing the computed singular values.    c
c         - u     : DOUBLE PRECISION                                  c
c                   Left singular vectors.                            c
c         - v     : DOUBLE PRECISION                                  c
c                   Right singular vectors.                           c    
c         - m, n  : INTEGER                                           c
c                   Size of the problem, m points, n spatial steps.   c
c         - nev   : INTEGER                                           c
c                   Number of required modes.                         c
c         - maxncv, ldu, ldv,maxnev, maxm are the leading dimension   c
c           of the different arrays and vectors.                      c
c                                                                     c
c      OUTPUT :                                                       c
c      ------                                                         c
c         - u     : left singular vectors (spatial modes) normalized. c
c         - v     : right ones (temporal modes) normalized.           c
c---------------------------------------------------------------------c

      implicit         none
      Integer
     &                 j, k, m, n, nev, maxncv, maxnev, maxm, ldu,
     &                 ldv, idamax  
      real
     &                 s(maxncv,2), u(ldu, maxnev), v(ldv,maxncv)
      real 	       tmp1, tmp2
      integer         isamax
c
      do j=1,nev
c
c  Cherche la composante de + grand module du mode temporel
c
         k=isamax(n,v(1,j),1)
         tmp1=v(k,j)
         tmp2=1.0/tmp1
c
c  Normalise le mode temporel
c
         call sscal(n,tmp2,v(1,j),1)
         tmp2=tmp1*s(j,1)
c
c  Normalise le mode spatial
c
         call sscal(m,tmp2,u(1,j),1)
      enddo
      
 101  format(i6)
 102  format(e16.8)
 103  format(a6)
 104  format(10e16.8)
 105  format(12e16.8)
      end
