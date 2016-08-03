!> A module which contains variables for BCS calculation.
module bcs
double precision::x2,gap
contains
!> This pairing routine allows the computation of the new occupational number after BCS convergency.
!> \param n is the number of the states considered for routine.
!> \param spener is he array containing n ordered singleparticle energies.
!> \param v2 is the array containing the <b>normalized</b> occupation numbers. 
!> \param rhobcs The rho matrix computed after BCS  
!> \param lpr is the logical flag for printing the convergency pattern of the chemical potential.
subroutine pairing(n,spener,v2,rhobcs,lpr)
use constants
use basis
use maths
implicit none
integer::i,j,k,l,n,jf,q,sizebloc
integer,parameter::limit=100,itbcs=100
double precision::x1,x2,x
double precision::gap_next,e
double precision::d,step,value
double precision::spener(n),v2(n)
double precision::rhobcs(lmin:lmax,jmin:jmax,maxval(n_red) + 1,maxval(n_red) + 1)
logical::lpr,s
sizebloc=maxval(n_red) + 1
e  =1.d-11
x1 = 10.d0;x2=40.d0
gap = 1.d0 !MeV
rhobcs = 0.d0
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
open(78,file="occ.dat")
call sorteigv(red_size,spener,v2)
do j=1,red_size
 write(78,'(i3,2f18.12)') j,spener(j),v2(j)
enddo
close(78)

! --------- Construction of the matrix density wih BCS solutions
 do l = lmin,lmax !Construction of the density matrix for all the possible blocks
  if (l == 0) then
     jf = 0
  else
     jf = 1
  endif
  do k=0,jf
  if (l == 0) then
    j = 1
   elseif (k==0) then
    j = 2*l - 1
   else
    j = 2*l+1
  endif
    do q=1,sizebloc
     rhobcs(l,j,q,q) = two*v2(tag_hf(q-1,l,j))
   enddo !q
  enddo !j
  enddo !l

end subroutine


!> A function computing the BCS kernel (the number equation).
!> \param n is the number of BCS states.
!> \param spener is an array containing n ordered singleparticle energies.
!> \param gap is a variable determining the BCS gap. 
!> \param lambda is chemical potential to be computed.
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
