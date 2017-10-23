subroutine valsvd(U,S,V,IEX,JEX,VAL,M,N,maxnev,maxncv,ldu,ldv,NN) 

implicit none
 integer,intent(in) :: maxnev,maxncv,ldu,ldv
 real, intent(in) :: U(ldu, maxnev),V(ldv,maxncv),S(maxncv,2)
! real, intent(in) :: U(ldu, nev),V(ldv,ncv),S(ncv,2)
 integer,intent(in) :: M,N,NN,IEX,JEX
 real,intent(out) :: VAL
 integer :: K


VAL=0 
! k=1
DO K=1,NN 

VAL=VAL+U(IEX,K)*S(K,1)*V(JEX,K) 

ENDDO


 
!write(6,*) 'u',U(IEX,k) 
end subroutine
