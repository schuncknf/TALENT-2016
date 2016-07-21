program hforsay
use constants
use lag_pol
use maths
use pot
use ho

implicit none
double precision::ml,x
double precision::test
integer::i,nl,n1,n2
double precision::test11(4)


call reader()
call lag_roots(ngauss,0.5d0,.true.)
!call gausslag(2*nbase,1,1,x)
!call lag_roots(nbase,0.0d0,.true.)
!read(*,*)
call hfsolver()
!do i=1,10
!write(*,*) 'nl,ml,x ?'
!read(*,*) nl,ml,x
!write(*,*) 'Laguerre',laguerre(nl,ml,x)
!enddo





end program
