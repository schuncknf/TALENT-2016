subroutine lag_roots(n)
implicit none
double precision,external::laguerre
double precision::e,x,step,val
double precision,dimension(:),allocatable::lag_zero
logical::t
integer::i,n
 n=2
 !lag_zero = 0.d0
 x = -1.0 ; step = 1.0e-6 ; e = 1.0e-9
  t = (laguerre(n,x) > 0)
  do while (x < 3.0)
    val = laguerre(n,x)
    if (abs(val) < e) then
      WRITE(*,"(A,F12.9)") "Root found at x =", x
      t = .not. t
    else if ((val > 0) .neqv. t) then
      write(*,"(A,F12.9)") "Root found near x = ", x
      t = .not. t
    end if
    x = x + step
  end do
!  lag_zero(i) = x
end
