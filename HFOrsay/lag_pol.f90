!Recursive defintion of laguerre polynomial
module lag_pol
double precision,dimension(:),allocatable::lag_zeros,lag_w
double precision,dimension(:,:),allocatable::temp
contains
subroutine lag_roots(n,m,pr)
use maths
use lag
implicit none
!double precision,external::laguerre
double precision::e,step,val,x_max,xi,x,x2,d,clag
integer,intent(in)::n
double precision,intent(in)::m
logical::t,pr
integer::i,j,limit
  write(*,*) '**** Laguerre Roots ****'
  allocate(lag_zeros(n))
  allocate(lag_w(n))
  allocate(temp(1,n+2))
  temp = 0.d0
  lag_zeros=0.d0
  lag_w=0.d0
  call lf_function_zeros(n,m,lag_zeros)
  do i=1,n
  xi=lag_zeros(i)
  if (pr) write(*,*) "Laguerre root",i,":",xi
  clag = 0.d0 
  call lf_function (1, n+1, m, xi, temp)
  clag = temp(1,n+2)
  if (m==0.d0) then
  !lag_w(i)=xi/((n+1)**2*laguerre(n+1,m,xi)**2)
  lag_w(i)=xi/((n+1)**2*clag**2)
  else
  lag_w(i)=xi*Gamma(dble(n+m+1))/(fac(n)*(n+1)**2*clag**2)
  endif
  if (pr) write(*,*) "Laguerre weight",i,":",lag_w(i)

  enddo
  write(*,*) '** End Laguerre Roots **'
end subroutine lag_roots
end module lag_pol




