      module globals
      implicit none
      integer :: dm=kind(1.d0)
      real (kind=dm), parameter :: h2m= 20.75, eps=1e-6
      real (kind=dm) :: R_max=10., R_min=0., R_box=6.,h=0.1,Em
      integer (kind=dm) :: Nmesh
      ! psi(grid_index)
      real (kind=dm),allocatable :: Psi(:) 
      logical :: comparison= .true., bisloop= .true.
      real (kind=dm) :: E_minus=0., E_plus=100.
      Nmesh=nint((R_max-R_min)/h)
      end module globals
