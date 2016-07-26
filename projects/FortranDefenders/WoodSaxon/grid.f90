module grid
implicit none
     integer, parameter :: wp=kind(1.0d0)
     real(wp), parameter :: pi = 3.14159265358979_wp
     real(wp), parameter :: e2 = 1.439646_wp
     real(wp), parameter :: hbar = 6.582119E-22_wp
     real(wp) :: h,conv,hbar22m,v0,nrad,vpb(2),r0,small,spin(2)
     real(wp), allocatable,dimension(:) :: meshpoints
     real(wp), allocatable, dimension(:,:) :: rho
     real(wp), allocatable, dimension(:,:,:,:,:) :: wavefunctions,wfl,wfr
     integer :: nbox, nodes, radius, lmax, welltype,nmax,njoin
     integer :: nn,np,nt
contains
     subroutine init_params

          namelist /box/ nbox,h
          namelist /squarewell/ welltype,nodes,v0,radius
          namelist /params/ r0,conv,hbar22m
          namelist /nucleus/ nn,np,lmax
          read(5,box)
          read(5,squarewell)
          read(5,params)
          read(5,nucleus)
          if ((welltype /= 1) .AND. (welltype /= 2)) write (*,*) "Put in a proper welltype!"
          nt = np+nn
          nrad = r0 * (nt)**(1._wp/3._wp)
          vpb(1) = -51.+33.*(nn-np)/nt
          vpb(2) = -51.-33.*(nn-np)/nt
          spin(1) = -0.5
          spin(2) = 0.5
          njoin = 200
          nmax = nn - np

          if (nmax.ge.0) then
          nmax = nn
          else
          nmax = np
          end if

     end subroutine init_params

     subroutine init_grids

          integer :: i
          small = 1E-20_wp
          allocate(meshpoints(0:nbox))
          meshpoints = (/ (real(i)*h,i=0,nbox) /)
          meshpoints(0) = small

     end subroutine init_grids

     subroutine init_wavefunctions

          allocate(wavefunctions(0:nbox,lmax,0:lmax,2,2),wfr(0:nbox,lmax,0:lmax,2,2),wfl(0:nbox,lmax,0:lmax,2,2),rho(0:nbox,4))

     end subroutine


end module grid
