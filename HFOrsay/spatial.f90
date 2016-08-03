!> A subroutine ploting spatial density in r-space
!> \param rho The density matrix
subroutine spatial_rho(rho)
use ho
use basis
use constants
use maths
use lag_pol
implicit none
integer::i,jf,l,j,n,ir,spatial_mesh
integer::n1,n2,k,sizebloc
double precision::r,part_spat,step
double precision::rho(lmin:lmax,jmin:jmax,maxval(n_red) + 1,maxval(n_red) + 1)
double precision,allocatable::rhor(:)
sizebloc=maxval(n_red) + 1

!------- Spatial Mesh
spatial_mesh = 60
step = 0.1d0
allocate(rhor(0:spatial_mesh))
rhor = 0.d0

!------- Loop over the Mesh
do ir=0,spatial_mesh
 r= step*ir
!------- Loop over blocks
 do l = lmin,lmax 
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
  do n1 = 1,sizebloc
   do n2 = 1,sizebloc
       rhor(ir) = rhor(ir) + rho(l,j,n1,n2)*ho_rad_wf(n1-1,l,r)*ho_rad_wf(n2-1,l,r) 
   enddo !Quantum number n1
   enddo !Quantum number n2
  enddo !j
 enddo !l
 enddo !Mesh

!-------- Particules number conservation and output
open(66,file='density.dat')
part_spat = 0.d0
do i=0,spatial_mesh
 r = i*step
 part_spat = part_spat + r**2*rhor(i)*step
 write(66,'(2f18.14)') r,rhor(i)
enddo
 write(*,*) "Particule Spatial-Density",part_spat
close(66)
end subroutine
