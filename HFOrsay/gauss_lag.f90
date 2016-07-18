

!subroutine gauss_lag(func,xi)
!implicit none
!double precision::func,xi,wi
!integer::i,j,n



!end subroutine



!Recursive defintion of laguerre polynomial
recursive function laguerre(n,x) result (lag)
implicit none
integer::n
double precision::x,lag
if (n == 0) then
    lag = 1
else if (n == 1) then
    lag = 1.d0-x
else
    lag = ((2*n-1-x)*laguerre(n-1,x)-(n-1)*laguerre(n-2,x))/n
endif
end



