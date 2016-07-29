      module Skyrme_force

      function Skyrme (xy, rho, rhon, rhop, psin, psip)
      implicit none
      use globals
      real(8) :: xy, kinetic, Spot
      real (8), dimension(:), allocatable :: rho, rhon, rhop
      allocate(rho(Nmesh), rhon(Nmesh), rhop(Nmesh))

      call Kindens(xy, taun, psin)
      call kindens(xy, taup, psip)
      tau(:) = taun(:) + taup(:) 

      do k=1, Nmesh
      kinetic(k) = h2m*(tau_n(k) + tau_p(k))

      Spot(k) = rho(k)*(t0 + (2.d0+alpha)*t3/12.d0*rho**alpha) -  & 
      (rhon(k) + rhop(k))*(t0/2.d0 +t3/12.d0*rho**alpha)-           &
      alpha*rho**(alpha-1.d0)*t3/24.d0*(rhop(k)**2 + rhon(k)**2) 

      enddo

      VSkyrme(:) = kinetic(:)+ Spot(:)

      end function Skyrme



      subroutine Kindens (xy, ktau, wf)
      use globals
      implicit none
      real(8)  xy
      integer(8) k, j
      real(8), dimension(:), allocatable :: ktau, wf
      allocate(ktau(0:Nmesh))

      ! mirar quines daquestes variables aniran a lapartat globals i
      ! quines es queden dins de la subrutina

      ! Calculation kinetic density for protons, neutrons and total
      
      ! TODO: sha de fer el calcul per als orbitals ocupats NOMES, 
      ! llavors fas la suma


      !material necessari: necessito un vector amb els nombres quantics
      !NECESSITO UNA SUBRUTINA QUE EM FACI LA DERIVADA DE LA FUNCIO PSI
      call derivative(psi, dpsi)

      do k=0, Nmesh !loop sobre els punts de la xarxa r
            do j=1, orbitals !loop over orbitals
      ktau(k)= = ktau(k) + (2.d0*j+1.d0)/4.d0/pi/xy**2*&
      ((dpsin(k) - psin(k)/xy)**2 + l(l+1.d0)/xy**2*psin(k)**2)

            enddo !loop over orbitals
      enddo ! loop over r
      return
      end function Tau





      subroutine derivative(array, darray)
      use globals
      implicit none
      real(8), dimension(:) , allocatable :: array
      integer(8) :: dime, i

      darray(0) = 0.d0 
      do i=1, Nmesh-2
      darray(i)= (2.d0*array(i+1) - 3.d0*array(i)/2.d0 - array(i+2)/2)/(2.d0*h)      
      enddo
      darray(Nmesh-1) = 0.d0
      darray(Nmesh) =0.d0

      return      
      end subroutine








      end module Skyrme_force
