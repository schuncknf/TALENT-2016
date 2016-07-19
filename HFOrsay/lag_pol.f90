!Recursive defintion of laguerre polynomial
module lag_pol
double precision,dimension(:),allocatable::lag_zeros,lag_w
contains
subroutine lag_roots(n,m,pr)
use lag
implicit none
!double precision,external::laguerre
double precision::e,x,step,val,x_max,xi
integer,intent(in)::n
double precision,intent(in)::m
logical::t,pr
integer::i,j
  write(*,*) '**** Laguerre Roots ****'
  open(10,file='roots_lag.dat')
  x = 0.0 ; step = 1.0e-6 ; e = 1.0e-9;x_max=1.0e4
  t = (laguerre(n,m,x) > 0)
  j = 0
  allocate(lag_zeros(n))
  allocate(lag_w(n))
  lag_zeros=0.d0
  lag_w=0.d0
  do while (x < x_max)
    val = laguerre(n,m,x)
    if (abs(val) < e) then
      j = j + 1
      if (pr) then
      WRITE(*,*) "Laguerre root",j," =", x
      write(10,*) x
      lag_zeros(j) = x
      endif
      t = .not. t
    else if ((val > 0) .neqv. t) then
      j = j + 1
      if (pr) then
      WRITE(*,*) "Laguerre root",j," =", x
      write(10,*) x
      lag_zeros(j) = x
      endif
      t = .not. t
    end if
    x = x + step
    if (j .eq. n) exit
  end do
  close(10)
  do i=1,n
  xi=lag_zeros(i)
  lag_w(i)=xi/((n+1)**2*laguerre(n+1,m,xi)**2)
  enddo
  write(*,*) '** End Laguerre Roots **'
end subroutine lag_roots
end module lag_pol




