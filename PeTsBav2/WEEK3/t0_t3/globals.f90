! module globals

module globals
   implicit none
   ! basic const
   integer :: Nmesh
   integer :: orbital
   integer :: n_max,l_max
   integer :: proton, neutron, nucleon
real(8),parameter :: pi=4.d0*atan(1.d0) ! hbar: (MeV*s), unit:(MeV),(fm),(s)
real(8),parameter :: t0=-1132.4d0, t3=23610.4d0, x0=0d0, x3=0d0, alpha=1d0 ! skyrme parameters
integer,parameter :: iter=10000
real(8) :: h2m, R_max, e2,  charge, epsi, dx, CMh2m, V_Co_d, V_Co_e
real(8) :: r0=1.27d0, a=0.67d0
logical :: bisloop=.true. , output_wave_func, CMcorrection
end module globals
