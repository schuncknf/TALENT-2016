!================================ arrange energy in ascending order =======================================

      subroutine arrage_energy(matrix, numorb)
      use globals
      implicit none
      real(kind=dm) :: temp(4)
      integer(8) :: iii, jjj ,numorb
      real(kind=dm) , dimension(numorb,4)::  matrix

      do iii=1, numorb-1
            do jjj=1, numorb-iii
                  if (matrix(jjj, 1).gt. matrix(jjj+1,1)) then
                  temp(:)= matrix(jjj, :)
                  matrix(jjj, :) = matrix(jjj+1, :)
                  matrix(jjj+1, :) = temp(:) 
                  endif
            enddo
      enddo
      end subroutine 

