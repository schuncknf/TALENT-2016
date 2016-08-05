
! module globals

module globals
   implicit none
   ! basic const
   integer :: Nmesh
   integer :: orbital
   integer :: n_max,l_max
   integer :: proton, neutron, nucleon
real(8),parameter :: pi=4.d0*atan(1.d0) ! hbar: (MeV*s), unit:(MeV),(fm),(s)
real(8) :: t0=-1132.4d0, t1=0d0, t2=0d0, t3=23610.4d0,&
     x0=0d0, x1=0d0, x2=0d0, x3=0d0, W0=125d0, alpha=1d0 ! skyrme parameters
!real(8) :: t0=-2483.45d0, t1=484.23d0, t2=-556.69d0, t3=13757d0,&
!     x0=0.776d0, x1=-0.317d0, x2=-1d0, x3=1.2630d0, W0=125, alpha=1/6d0 ! skyrme parameters
integer,parameter :: iter=10000
real(8) :: h2m, R_max, e2,  charge, epsi, dx, CMh2m, corr
real(8) :: r0=1.27d0, a=0.67d0, g=0.2d0
logical :: bisloop=.true. , output_wave_func, CMcorrection, flag, BCS_cal
integer :: magic_n, magic_p, orbital_down, orbital_up,  NNN, PPP, BCS_level
character(10) :: skforce
end module globals
