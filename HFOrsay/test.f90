function f(n,x) result (te)
implicit none
integer::n
double precision::te,x
te=n*x
end function
program test
use lag
use lag_pol
use constants
use maths
use ho
use pot
implicit none
double precision::x,ri,rj,l
double precision,dimension(:,:),allocatable::z
integer::n,i,nsize,j,m
!external f
!write(*,*) 'n,m ?'
!read(*,*) n,l
n = 10
!l=1.d0/2.d0
l=0.5d0
do i=-100,100
ri = 0.1*i
 do j=-100,100
 rj = 0.1*j
write(11,*) ri,rj,potential(ri,rj,one,half)
write(12,*) ri,rj,minnesota(ri,rj,v0r,kr)
enddo
write(11,*) 
write(12,*)
enddo

!write(*,*) dsqrt(2**(n+l+2)*fac(n)/(dsqrt(pi)*ffac(2*n+2*l+1)))
!write(*,*) fac(n),ffac(n)
!allocate(lag_zero(n))
!lag_zero = 0.d0


call lag_roots(n,l,.true.)
write(*,*) 'ok'
!f(x) = laguerre(1,1/2.d0,x)*laguerre(1,1/2.d0,x)
call gausslag(n,1,2,x)
write(*,*) x
!do i=1,n
 !read(11,*) lag_zero(i)
!enddo
!enddo
end
