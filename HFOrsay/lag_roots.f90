subroutine lag_roots(n,pr)
implicit none
double precision,external::laguerre
double precision::e,x,step,val,x_max
integer,intent(in)::n
logical::t,pr
integer::i,j
  write(*,*) '**** Laguerre Roots ****'
  open(10,file='roots_lag.dat')
  x = 0.0 ; step = 1.0e-6 ; e = 1.0e-9;x_max=1.0e4
  t = (laguerre(n,x) > 0)
  j = 0
  do while (x < x_max)
    val = laguerre(n,x)
    if (abs(val) < e) then
      j = j + 1
      if (pr) then
      WRITE(*,*) "Laguerre root",j," =", x
      write(10,*) x
      endif
      t = .not. t
    else if ((val > 0) .neqv. t) then
      j = j + 1
      if (pr) then
      WRITE(*,*) "Laguerre root",j," =", x
      write(10,*) x
      endif
      t = .not. t
    end if
    x = x + step
    if (j .eq. n) exit
  end do
  close(10)
  write(*,*) '** End Laguerre Roots **'
end
