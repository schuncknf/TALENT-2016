module grid
implicit none
     integer, parameter :: wp=kind(1.0d0)
     real(wp), parameter :: pi = 3.14159265358979_wp
     real(wp), parameter :: e2 = 1.439646_wp
     real(wp), parameter :: hbar = 6.582119E-22_wp
     real(wp), parameter :: a = 0.67_wp
     real(wp), parameter :: vso = 23_wp
     real(wp) :: t0,x0,t1,x1,t2,x2,t3,x3,sig,w0
     real(wp) :: cmcorr
     real(wp) :: h,conv,hbar22m,v0,nrad,vpb(2),r0,small,spin(2)
     real(wp) :: a0r0,a1r1,a0s0,a1s1,a0tau0,a1tau1,a0t0,a1t1,a0r0p,a1r1p,& 
                & a0s0p,a1s1p,cddr0,cddr1,cdds0,cdds1,cso0,cso1
     real(wp), allocatable,dimension(:) :: meshpoints,ucoul
     real(wp), allocatable, dimension(:,:) :: rho,tau,jsc,drho,dtau,djsc,ddrho
     real(wp), allocatable, dimension(:,:) ::uc,umr,udd,uso
     real(wp), allocatable, dimension(:,:,:,:,:) :: wavefunctions,wfl,wfr
     integer :: nbox, nodes, radius, lmax, welltype,nmax,njoin
     integer :: nn,np,nt,icoul,icm
     logical :: j2terms
contains
     subroutine init_params

          namelist /box/ nbox,h
          namelist /squarewell/ welltype,nodes,v0,radius
          namelist /params/ r0,conv,hbar22m
          namelist /nucleus/ nn,np,lmax
          namelist /interaction/ t0,x0,t1,x1,t2,x2,t3,x3,sig,w0, &
                   j2terms,icoul,icm
          read(5,box)
          read(5,squarewell)
          read(5,params)
          read(5,nucleus)
          read(5,interaction)
          if ((welltype /= 1) .AND. (welltype /= 2)) write (*,*) "Put in a proper welltype!"
          nt = np+nn
          nrad = r0 * (nt)**(1._wp/3._wp)
          vpb(1) = -51.+33.*(nn-np)/nt
          vpb(2) = -51.-33.*(nn-np)/nt
          cmcorr = 1 - (1/nt)
          spin(1) = -0.5
          spin(2) = 0.5
          njoin = 200
          nmax = nn - np

          if (nmax.ge.0) then
          nmax = nn
          else
          nmax = np
          end if
          !!
          !setting the coupling constants
          !!
          a0r0 = 3._wp/8._wp * t0 
          a1r1 = - 1._wp/4._wp * t0 * ( 1._wp/2._wp + x0 )
          a0s0 = - 1._wp/4._wp * t0 * ( 1._wp/2._wp - x0 )
          a1s1 = - 1._wp/8._wp * t0 
          !
          a0tau0 = 3._wp/16._wp * t1 + 1._wp/4._wp * t2 * ( 5._wp/4._wp + x2 ) 
          a1tau1 = - 1._wp/8._wp * t1 * ( 1._wp/2._wp + x1 ) + &
                 & 1._wp/8._wp * t2 * ( 1._wp/2._wp + x2 )
          !          
          a0t0 = - 1._wp/8._wp *t1* ( 1._wp/2._wp - x1 ) + &
                 & 1._wp/8._wp * t2 * ( 1._wp/2._wp + x2 )
          a1t1 =  1._wp/16._wp * (t2-t1) 
          !
          a0r0p = - 9._wp/64._wp * t1 + 1._wp/16._wp * t2 *( 5._wp/4._wp + x2 )
          a1r1p = 3._wp/32._wp * t1 * ( 1._wp/2._wp + x1 ) + &
                 & 1._wp/32._wp * t2 * ( 1._wp/2._wp + x2 )
            
          a0s0p = 3._wp/32._wp * t1 * ( 1._wp/2._wp - x1 ) + &
                 & 1._wp/32._wp * t2 * ( 1._wp/2._wp + x2 )
          a1s1p = 3._wp/64._wp * t1 + 1._wp/64._wp * t2
    	  !
          cddr0 =  t3 / 16.0_wp
          cddr1 = - 1._wp/24._wp * t3 * ( 1._wp/2._wp + x3 )
          cdds0 = - 1._wp/24._wp * t3 * ( 1._wp/2._wp - x3 )
          cdds1 = - 1._wp/48._wp * t3
          !
          cso0 = - 3._wp/4._wp * w0
          cso1 = - 1._wp/4._wp * w0

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

     end subroutine init_wavefunctions

    subroutine init_fields

          allocate(uc(0:nbox,2),umr(0:nbox,2),udd(0:nbox,2),uso(0:nbox,2),ucoul(0:nbox))

    end subroutine init_fields

end module grid
