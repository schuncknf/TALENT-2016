program test
implicit none
double precision::x
double precision,external::laguerre
double precision,dimension(:),allocatable::lag_zero
integer::n,i,nsize,j
logical::chk_lag
inquire(file="roots_lag.dat", exist=chk_lag)
if (chk_lag) open(11,file="roots_lag.dat")
! do j=1,10
!  read(11,*) lag_zero(j)
! enddo

write(*,*) "n?"
read(*,*) n
allocate(lag_zero(n))
lag_zero = 0.d0
call lag_roots(n,.true.)
do i=1,n
 read(11,*) lag_zero(i)
enddo
write(*,*) lag_zero
!enddo
end

