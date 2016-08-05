!! Function to calculate the Coulomb potential
!!
!!*******************************************************************************
!! N  : number of the point in the mesh such as dr*N=r, until where to integrate
!!      the two first integrals (in)
!! rho_p : proton density (in)
!!********************************************************************************
      real(8) function Vcoulomb(N, rho_p)
      use globals
      implicit none
      integer i, N
      real(8) x, int_a, int_c, exchange, rho_p(0:Nmesh)
      int_a=0.d0
      int_c=0.d0

      do i=1, N
      x= i*dx
      int_a= int_a+ rho_p(i)*x**2*dx
      !int_b = int_b + rho_p(x)*x*dx
      enddo
      do i=N+1, Nmesh
      x=i*dx
      int_c = int_c + rho_p(i)*x*dx
      enddo
      exchange = e2*(3.d0/pi)**(1./3.)*rho_p(N)**(1./3.)

      V_Co_d(N) = 4.d0*pi*e2*(1.d0/(N*dx)*int_a +int_c)

      Vcoulomb =  4.d0*pi*e2*(1.d0/(N*dx)*int_a +int_c) -exchange
      end function
