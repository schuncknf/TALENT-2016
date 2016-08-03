!> A maths module wich contains basic tools and mathematical functions
module maths
contains
!> A logarithmic implementation of the Factorial
!> \param n integer n!
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
!> A logarithmic recursive implementation of the Double factorial
!> \param n integer n!!
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
!> An array sorting function
!> \param n the size of the array
!> \param a the array to sort
subroutine sort(n,a)
  implicit none
  integer n,i,j
        double precision a(n),x
        do 30 i=2,n
        x=a(i)
        j=i
   10   j=j-1
        if(j.eq.0 .or. a(j).le.x) go to 20
        a(j+1)=a(j)
        go to 10
   20   a(j+1)=x
   30   continue
        end subroutine

!> An simultaneous array sorting function
!> \param n the size of the array
!> \param indexarray the array used to index the array to sort 
!> \param arraytosor the array to be sort
subroutine sorteigv(n,indexarray,arraytosort)
implicit none
integer::n,i,j,k,l
double precision::val1,val2
double precision::indexarray(n),arraytosort(n),tempar(n),te(n)
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
  te(i) = arraytosort(sorted(i))
enddo
arraytosort= 0.d0
arraytosort=te
call sort(n,indexarray)
end subroutine
!>This Function compute the trace of a squared matrix
!> \param n the size of the matrix
!> \param M the input matrix
function trace(M,n) result(tr)
implicit none
integer::n,i
double precision::M(n,n),tr
tr = 0.d0
do i = 1,n
tr = tr + M(i,i)
enddo
end function
end module maths
