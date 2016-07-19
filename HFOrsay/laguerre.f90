module lag
contains
recursive function laguerre(n,m,x) result (lag)
implicit none
integer::n
double precision::x,lag,m
if (n == 0) then
    lag = 1
else if (n == 1) then
    lag = 1.d0+m-x
else
    lag = ((2*n-1+m-x)*laguerre(n-1,m,x)-(n-1+m)*laguerre(n-2,m,x))/n
endif
end function laguerre
end module lag

