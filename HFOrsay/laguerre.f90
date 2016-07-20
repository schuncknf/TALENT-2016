module lag
contains
function laguerre(n,m,x) result (lag)
implicit none
integer,intent(in)::n
double precision,intent(in)::x,m
double precision::l,lag
double precision::temp_p(n,n)
!allocate(temp_p(1,1:n+1))
temp_p=0.d0
lag = 0.d0
l=m
call lf_function (1, n+1, l, x, temp_p)
lag = temp_p(1,n)
!deallocate(temp_p)
!if (n == 0) then
!    lag = 1
!else if (n == 1) then
!    lag = 1.d0+m-x
!else
!    lag = ((2*n-1+m-x)*laguerre(n-1,m,x)-(n-1+m)*laguerre(n-2,m,x))/n
!endif
end function laguerre
end module lag

