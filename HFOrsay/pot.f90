module pot
  contains
function minessota(r1,r2) result(v)
use constants
implicit none
double precision::v1,v2,r1,r2,v,ri,rj
ri=dsqrt(r1)*bosc
rj=dsqrt(r2)*bosc
v1 =v0r/kr*(exp(-kr*(ri**2+rj**2-2*ri*rj))-exp(-kr*(ri**2+rj**2+2*ri*rj)))/(4.d0*ri*rj) 
v2=-v0s/ks*(exp(-ks*(ri**2+rj**2-2*ri*rj))-exp(-ks*(ri**2+rj**2+2*ri*rj)))/(4.d0*ri*rj)
v=(v1+v2)*half
end function

subroutine kinetic(n,nr,nl,kin)
 use constants
 use basis
 implicit none
 double precision::kin(n,n)
 integer::nr(n),nl(n)
 integer::n,i
 integer::n1,l1
 kin=0.d0
 do i=1,n
!   kin(i,i) = (2.d0*nr(i)+nl(i)+1.5d0)*ama*2.d0/(bosc**2)
   !kin(i,i) = (2.d0*(i-1)+1.5d0)*ama*2.d0/(bosc**2)
    n1 = n_red(i)
    l1 =l_red(i)
     kin(i,i) = (2.d0*n1+l1+1.5d0)*ama*2.d0/(bosc**2)
 enddo

! ---- Petar
! kin = 0.d0
! do n1 = 1, nt (number of states) 
! read quantum numbers (n1,l1,j1,...) of state i
!	 kin(i,i) = (2*n1+l1+1.5d0)*ama*2.d0/(bosc**2) ! there is no (i-1) as before because state i=1 can (should) contain n=0 state
! enddo
! ---- Petar


!
 end subroutine
end module
