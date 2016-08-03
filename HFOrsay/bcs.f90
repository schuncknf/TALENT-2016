!> A module wich contains variables for BCS calculation
module bcs
double precision::x2,gap
contains
!> This pairing routine allow the computation of the new occ. number after BCS convergency
!> \param n The number of the states considered for routine
!> \param spener an array containing the n ordered single particles energies
!> \param v2 The array containing the <b>normalized</b> occupation numbers 
!> \param lpr a logical flag to print the convergency pattern of the chemical potential
subroutine pairing(n,spener,v2,lpr)
use constants
use basis
implicit none
integer::i,j,k,l,n
integer,parameter::limit=100,itbcs=100
double precision::x1,x2,x
double precision::gap_next,e
double precision::d,step,value
double precision::spener(n),v2(n)
logical::lpr,s
e  =1.d-11
x1 = 10.d0;x2=40.d0
gap = 1.d0 !MeV
i=0
do k=1,red_size
enddo
do k = 1,itbcs
if (k .eq. 1) then
  x2 = 17.d0
 else
 do 
  if (i > limit) then
    exit
  end if
  d = (x2 - x1) / (bcs_lambda(red_size,spener,gap,x2) - bcs_lambda(red_size,spener,gap,x1)) * bcs_lambda(red_size,spener,gap,x2)
  if (abs(d) < e) then
    exit    
  end if
  x1 = x2
  x2 = x2 - d
  i = i + 1
  enddo
  endif
  if(lpr) write(*,*) "It: ",k,"Lambda",x2


gap_next = 0.d0
do j=1,red_size
gap_next = gap_next + g_pair*half*gap/(dsqrt((spener(j)-x2)**2+gap**2))
enddo
if (abs(gap - gap_next) .le. e) then
exit
else
gap = gap_next
endif
enddo
do j=1,red_size
 v2(j) = half*(one-(spener(j)-x2)/(dsqrt((spener(j)-x2)**2+gap**2)))
enddo
end subroutine


!> A function computing the BCS kernel (The number equation)
!> \param n The number of BCS states
!> \param spener an array containing the n ordered single particles energies
!> \param gap a variable giving the BCS gap 
!> \param the too be computed checmical potential
function bcs_lambda(n,spener,gap,lambda) result(nbcs)
use basis
use maths
use constants
implicit none 
integer::i,n
double precision::gap
double precision::spener(n)
double precision::nbcs,lambda
nbcs = 0.d0
do i = 1,red_size
    nbcs = nbcs + (one - (spener(i) - lambda)/(dsqrt((spener(i)-lambda)**2+gap**2))) 
enddo 
nbcs = npart - nbcs
end function

end
