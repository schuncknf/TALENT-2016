subroutine tbme(n1,n2,n3,n4,resu,pr)
use ho
use lag_pol
use constants
use pot
use maths
implicit none
integer::i,n1,n2,n3,n4,n,j,l1,l2,l3,l4
double precision::inte1,inte2,resu
double precision::coeffi,coeffj
double precision::xxi,xi,xxj,xj
double precision::wi,wj
double precision::orth,lag1,lag2
double precision::a1,a2,a3,a4
double precision::norm_lag,norm_lag1
double precision::testw,ri,rj
double precision::nosc1,nosc2,nosc3,nosc4
logical::pr
l1=0;l2=0;l3=0;l4=0
a1=dble(l1)+0.5d0
a2=dble(l2)+0.5d0
a3=dble(l3)+0.5d0
a4=dble(l4)+0.5d0
inte1=0.d0
inte2=0.d0
orth=0.d0
norm_lag=Gamma(dble(n1)+1.d0+a1)/fac(n1)
norm_lag1=Gamma(dble(n4)+1.d0+a1)/fac(n4)
nosc1=ho_norm(n1,0.d0)
nosc2=ho_norm(n2,0.d0)
nosc3=ho_norm(n3,0.d0)
nosc4=ho_norm(n4,0.d0)
testw=0.d0
do i=1,ngauss
    wi=lag_w(i)
    xxi=(lag_zeros(i))
    ri = dsqrt(xxi)*bosc
    lag1 = laguerre(n1,a1,xxi)*laguerre(n2,a3,xxi)
 do j=1,ngauss
    wj=lag_w(j)
    testw = testw + wj
    xxj=(lag_zeros(j))
    rj = dsqrt(xxj)*bosc
    lag2 = laguerre(n3,a2,xxj)*laguerre(n4,a4,xxj)
    orth = orth + wi*wj*lag1*lag2
    !inte1 = inte1 + wi*wj*lag1*lag2*minnesota(ri,rj)
    inte2 = inte2 + wi*wj*(potential(xxi,xxj,v0r,kr)-potential(xxi,xxj,v0s,ks))*lag1*lag2
        enddo
enddo
orth = orth*nosc1*nosc2*nosc3*nosc3/4.d0
resu = inte2
resu = resu*nosc1*nosc2*nosc3*nosc3/4.d0
!endif
if (pr .and. orth .gt. 0.001d0) then
write(*,'(a,4i3,f20.14)') "n1,n2,n3,n4",n1,n2,n3,n4,orth
endif
end
