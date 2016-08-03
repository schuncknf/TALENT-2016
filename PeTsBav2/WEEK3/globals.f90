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
real(8) :: h2m, R_max, e2,  charge, epsi, dx, CMh2m
real(8) :: r0=1.27, a=0.67
logical :: bisloop=.true. , output_wave_func, CMcorrection, flag
 integer :: magic_n, magic_p, orbital_down, orbital_up,  NNN, PPP
end module globals
