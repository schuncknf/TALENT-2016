program test
use lag
use lag_pol
use constants
use maths
use ho
use pot
implicit none
double precision::x,ri,rj,l
real(8) ::fx
integer(4)::n

double precision,dimension(:,:),allocatable::z
integer::i,nsize,j,m
logical::chk_lag
!inquire(file="roots_lag.dat", exist=chk_lag)
!if (chk_lag) open(11,file="roots_lag.dat")
! do j=1,10
!  read(11,*) lag_zero(j)
! enddoho_rad_wf(n,l,r)

!allocate(z(1,2+1))
!call lf_function (1, 2, 0.5d0, 2.1d0, z)
write(*,*) 'n,m ?'
read(*,*) n,l
!allocate(z(n))
!call lf_function_zeros(n,m,z)
l=1.d0/2.d0
do i=-100,100
ri = 0.1*i
 do j=-100,100
 rj = 0.1*j
write(11,*) ri,rj,potential(ri,rj,one,half)
!write(12,*) ri,laguerre(4,5.d0,ri)
!write(13,*) ri,laguerre(n+2,0.d0,ri)
!write(11,*)  ri,ho_rad_wf(n,l,ri)
!write(12,*)  ri,ho_rad_wf(n+1,l,ri)
!write(13,*)  ri,ho_rad_wf(n+2,l,ri)
enddo
write(11,*) 
enddo

!write(*,*) dsqrt(2**(n+l+2)*fac(n)/(dsqrt(pi)*ffac(2*n+2*l+1)))
!write(*,*) fac(n),ffac(n)
!allocate(lag_zero(n))
!lag_zero = 0.d0
call lag_roots(n,0.d0,.true.)
write(*,*) 'ok'
call gausslag(n,m)
!do i=1,n
 !read(11,*) lag_zero(i)
!enddo
!enddo
end

