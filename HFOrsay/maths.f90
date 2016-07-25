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
