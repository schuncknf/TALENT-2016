program test
use lag
use lag_pol
use constants
use maths
use ho
implicit none
double precision::x,ri

!double precision,dimension(:),allocatable::lag_zero
integer::n,m,i,nsize,j,l
logical::chk_lag
!inquire(file="roots_lag.dat", exist=chk_lag)
!if (chk_lag) open(11,file="roots_lag.dat")
! do j=1,10
!  read(11,*) lag_zero(j)
! enddoho_rad_wf(n,l,r)

!write(*,*) "n,m,r?"
!read(*,*) n,l,x
n=0
l=0
do i=-1000,1000
ri = 0.01*i
!write(11,*) ri,laguerre(4,0.d0,ri)
!write(12,*) ri,laguerre(4,5.d0,ri)
!write(13,*) ri,laguerre(n+2,0.d0,ri)
write(11,*)  ri,ho_rad_wf(n,l,ri)
write(12,*)  ri,ho_rad_wf(n+1,l,ri)
write(13,*)  ri,ho_rad_wf(n+2,l,ri)
enddo

!write(*,*) dsqrt(2**(n+l+2)*fac(n)/(dsqrt(pi)*ffac(2*n+2*l+1)))
write(*,*) fac(n),ffac(n)
!allocate(lag_zero(n))
!lag_zero = 0.d0
call lag_roots(n,dble(l),.true.)
!do i=1,n
 !read(11,*) lag_zero(i)
!enddo
write(*,*) lag_zeros
write(*,*) 'Weights'
write(*,*) lag_w
!enddo
end

