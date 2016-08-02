!================================ arrange energy in ascending order =======================================

      subroutine arrage_energy(matrix, numorb, wavef)
      use globals
      implicit none
      real(kind=dm) :: temp(4)
      real(kind=dm) :: tempwave(Nmesh)
      integer(8) :: iii, jjj ,numorb
      real(kind=dm) , dimension(numorb,4)::  matrix
      real(kind=dm) , dimension(numorb,Nmesh) :: wavef

      do iii=1, numorb-1
            do jjj=1, numorb-iii
                  if (matrix(jjj, 1).gt. matrix(jjj+1,1)) then
                  temp(:)= matrix(jjj, :)
                  matrix(jjj, :) = matrix(jjj+1, :)
                  matrix(jjj+1, :) = temp(:)
                  tempwave(:)= wavef(jjj, :)
                  wavef(jjj, :) = wavef(jjj+1, :)
                  wavef(jjj+1, :) = tempwave(:)                  
                  endif
            enddo
      enddo
      end subroutine 

