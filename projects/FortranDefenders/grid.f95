module grid
implicit none
     integer, parameter :: wp=kind(1.0d0)
     real(wp), parameter :: pi = 3.14159265358979_wp
     real(wp), parameter :: e2 = 1.439646_wp
     real(wp), parameter :: hbar = 6.582119E-22_wp
     real(wp) :: h,conv,hbar22m,v0,nrad,vpb,r0,small
     real(wp), allocatable,dimension(:) :: meshpoints
     real(wp), allocatable, dimension(:,:,:,:) :: wavefunctions,wfl,wfr
     integer :: nbox, nodes, radius, lmax
     integer :: nn,np,nt
contains
     subroutine init_params

          namelist /input/ nbox,h,nodes,v0,radius,r0,conv,hbar22m,nn,np,lmax
          read(5,input)
          nt = np+nn
          nrad = r0 * (nt)**(1._wp/3._wp)
          vpb = -51.+33.*(nn-np)/nt

     end subroutine init_params

     subroutine init_grids

          integer :: i
          small = 1E-20_wp
          allocate(meshpoints(0:nbox))
          meshpoints = (/ (real(i)*h,i=0,nbox) /)

     end subroutine init_grids

     subroutine init_wavefunctions

          allocate(wavefunctions(0:nbox,0:lmax,2,2),wfr(0:nbox,0:lmax,2,2),wfl(0:nbox,0:lmax,2,2))

     end subroutine


end module grid
