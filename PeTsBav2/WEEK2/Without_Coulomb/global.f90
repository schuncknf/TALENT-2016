!     Global variables for the finite square well potential in main.f90
      module globals
      implicit none
      integer, parameter :: dm=kind(1.d0)
      real (kind=dm) ::h2m= 20.75, epsi=1e-6
      real (kind=dm) :: R_max=6., R_min=0., R_box=6.,h=0.1,Em, nmfactor, a=0.67, r0=1.67
      integer  :: Nmesh, NN=126, ZZ=82
      real (kind=dm),allocatable :: rho(:) ,Psi(:)! psi(grid_index)
      logical :: comparison= .true., bisloop= .true.
      real (kind=dm) :: E_minus=-150., E_plus=0., E_left, E_right
      integer(kind=dm) :: numnodes ! the number of excited state, that has to be the same as number of nodes in wavefunction
      integer(kind=dm) :: cnodes ! nodes of wavefuntion 
      integer(kind=dm) :: i,ii, uu
      integer(kind=dm) :: n, l, orbital
      integer(kind=dm) :: n_max, l_max,nl2j
      real(kind=dm)    :: sm, j
      real(kind=dm),allocatable :: k_sq(:), Enl2j(:,:)   
      character(99),allocatable :: filename_wave_func(:)
      real(kind=dm) :: x 
      end module globals



