module maths
contains
function fac(n) result(fak)
implicit none
integer::i,n
double precision::lfak,fak
if (n==0) then
 fak=1.d0
else
 lfak = 0.d0
 do i=1,n
 lfak = lfak + log(dble(i))
 enddo
 fak = exp(lfak)
endif
end function fac
function ffac(n) result (ffak)
implicit none
integer::n,k
double precision::ffak
if (mod(n,2) == 0) then
  k=n/2
  ffak = 2.d0**(dble(k))*fac(k)
 else
  k=(n+1)/2
  ffak = fac(2*k)/(2**(dble(k))*fac(k))
endif
end function ffac

integer function delta(i,j)
  integer :: i, j
  delta=0; if(i==j) delta=1
end function delta

!recursive function ffac(n) result (ffak)
!implicit none
!integer::n,ffak
!if (n .gt. -2 .and. n .le. 0) then
!    ffak = 1
!else
!    ffak = n*ffac(n-2)
!endif
!end function ffac

!function gausslag(n,func) result(inte)
!use lag_pol
!use constants
!implicit none
!double precision,external::func
!double precision::inte,wi,xxi,xi
!integer::i,n
!inte=0.d0
!do i=1,n
!wi=lag_w(i)
!xxi=(lag_zeros(i))
!inte = inte + wi*func(xxi)
!enddo
!end function gausslag
 SUBROUTINE SORT(N,A)
  IMPLICIT NONE
  INTEGER N,I,J
        DOUBLE PRECISION A(N),X
        DO 30 I=2,N
        X=A(I)
        J=I
   10   J=J-1
        IF(J.EQ.0 .OR. A(J).LE.X) GO TO 20
        A(J+1)=A(J)
        GO TO 10
   20   A(J+1)=X
   30   CONTINUE
        END subroutine

subroutine sorteigv(n,indexarray,arraytosort)
implicit none
integer::n,i,j,k,l
double precision::val1,val2
double precision::indexarray(n),arraytosort(n,n),tempar(n),te(n,n)
integer::sorted(n)

tempar=indexarray
sorted =0
call sort(n,tempar)
do i=1,n
 val1=indexarray(i)
 do j=1,n
         val2=tempar(j)
 if (val1 .eq. val2) sorted(j)=i
 enddo
enddo
do i=1,n
 do j=1,n
  te(i,j) = arraytosort(sorted(i),j)
 enddo
enddo
arraytosort= 0.d0
arraytosort=te
call sort(n,indexarray)


end subroutine





end module maths
