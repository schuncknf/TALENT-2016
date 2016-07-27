subroutine tbme(n1,n2,n3,n4,resu,pr,iflag)
use ho
use lag_pol
use constants
use pot
use maths
use basis
implicit none
integer::i,n1,n2,n3,n4,j,iflag
integer::l1,l2,l3,l4
integer::m1,m2,m3,m4
integer::j1,j2,j3,j4
double precision::inte1,inte2,resu
double precision::coeffi,coeffj
double precision::xxi,xi,xxj,xj
double precision::wi,wj
double precision::orth,lag1,lag2
double precision::a1,a2,a3,a4
double precision::norm_lag,norm_lag1
double precision::testw,ri,rj
double precision::nosc1,nosc2,nosc3,nosc4
integer :: stat
integer :: q1,q2,q3,q4
character(len=100) :: buf
logical::pr,fex
if (iflag == 0) then
l1=0;l2=0;l3=0;l4=0
a1=dble(l1)+0.5d0
a2=dble(l2)+0.5d0
a3=dble(l3)+0.5d0
a4=dble(l4)+0.5d0
inte2=0.d0
orth=0.d0
nosc1=ho_norm(n1,0.d0)
nosc2=ho_norm(n2,0.d0)
nosc3=ho_norm(n3,0.d0)
nosc4=ho_norm(n4,0.d0)
do i=1,ngauss
    wi=lag_w(i)
    xxi=(lag_zeros(i))
    lag1 = laguerre(n1,a1,xxi)*laguerre(n2,a2,xxi)
 do j=1,ngauss
    wj=lag_w(j)
    xxj=(lag_zeros(j))
    lag2 = laguerre(n3,a3,xxj)*laguerre(n4,a4,xxj)
    orth = orth + wi*wj*lag1*lag2
    inte2 = inte2 + wi*wj*(minessota(xxi,xxj))*lag1*lag2
        enddo
enddo
orth = orth*nosc1*nosc2*nosc3*nosc4/4.d0
resu = inte2
resu = resu*nosc1*nosc2*nosc3*nosc4/4.d0
if (pr .and. orth .gt. 0.001d0) then
write(*,'(a,4i3,f20.14)') "n1,n2,n3,n4",n1,n2,n3,n4,orth
endif
elseif (iflag == 1) then
resu = tbme_ext(n1,n2,n3,n4)
endif
end
