!Subroutine to arrange the WoodsSaxon eigenvalues in ascending order
! 

subroutine arrage_energy(matrix,wavef)
   use globals
   implicit none
   integer, parameter :: dm=kind(1.d0)
   real(kind=dm) :: temp(4)
   real(kind=dm) :: tempwave(0:Nmesh)
   integer(8) :: iii, jjj
   real(kind=dm) , dimension(orbital,4)::  matrix
   real(kind=dm) , dimension(0:Nmesh,orbital) :: wavef
    
   do iii=1, orbital-1
      do jjj=1, orbital-iii
         if (matrix(jjj, 1).gt. matrix(jjj+1,1)) then
           temp(:)= matrix(jjj, :)
           matrix(jjj, :) = matrix(jjj+1, :)
           matrix(jjj+1, :) = temp(:)
                  
           tempwave(:)= wavef(:, jjj)
           wavef(:,jjj) = wavef(:,jjj+1)
           wavef(:,jjj+1) = tempwave(:)                  
         endif
      enddo
    enddo
end subroutine 


