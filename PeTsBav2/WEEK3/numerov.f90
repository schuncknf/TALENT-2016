module      numerov_module

!primer parametre que entra es matriu, sestora a la posicio E(orbital, 4)
! segon parametre que entra es matriu(orbital, r)

      subroutine get_eigenvv (orbi)
      use globals
      implicit none
      real(8) :: E_left, E_right, Em, cnodes , nmfactor
      real(8) :: x
      real(8) , allocatable :: densi(Nmesh), k_sq(Nmesh)
      integer(8) :: i

                  E_left=E_minus
                  E_right= E_plus
                  psi(0) = 0.
                  psi(1) = 0.1
                  !begin bissection method
                  do while (bisloop) ! begin loop bissection method
                      if (abs(E_right-E_left)<epsi) exit
                      Em=(E_right+E_left)/2.                  
                      cnodes=0     
                      !begin numerov method
                      Do i = 1,Nmesh-1
                              x=(i)*h
                              k_sq(i) = (Em-SKYRMEEEE)/h2m -1.*(l*(l+1))/x**2
                      Enddo

                      Do i = 2,Nmesh
                              psi(i) = (2.*(1.-5.*h**2*k_sq(i-1)/12.)*psi(i-1) &
                              -(1.+h**2*k_sq(i-2)/12.)*psi(i-2))/(1.+h**2*k_sq(i)/12.)
                              if (psi(i-1)*psi(i)<0)  cnodes=cnodes+1 ! we find the number of nodes
                      Enddo
                      !end numerov method
 
                      if(cnodes.Gt.n) Then 
                              E_right = Em
                      else 
                              E_left = Em
                      Endif
                  end do! end loop bissection method

                  Enl2j(orbi, 1)=E_right !eigenvalue
                  ! normalization of the eigenfunction
                  ! We recalculate psi because we want to find the function that does
                  ! not go to infinite SHOOTING METHOD WITH RIGHT HAND SIDE EQUAL TO 0

                  cnodes = 0
                  Do i = 1,Nmesh-1
                        x=(i)*h
                        Vpot = Skyrme (xy)
                        k_sq(i) = (Enl2j(orbi, 1)-Vpot)/h2m - 1.*((l*(l+1))/x**2) 
                  Enddo
            
                 psi(0) = 0.
                 psi(1) = 0.1
      
                 Do i = 2,Nmesh
                        psi(i) = (2.*(1.-5.*h**2*k_sq(i-1)/12.)*psi(i-1) &
                        -(1.+h**2*k_sq(i-2)/12.)*psi(i-2))/(1.+h**2*k_sq(i)/12.)
                        if (psi(i-1)*psi(i)<0)  cnodes=cnodes+1 ! we find the number of nodes  
                        if (psi(i)*psi(i-1)<=0d0.and.cnodes>n) then
                              psi(i) = 0.
                        end if
                 end do

                 ! normalization
                 nmfactor = 0d0
                 densi(:) = abs(psi( :))**2
                 call simpson(Nmesh, h, densi, nmfactor)
                 psi(:) = psi(:)/sqrt(nmfactor)    !wave func
      end subroutine get_eigenvv



      subroutine calculate_dens()


      end subroutine








end module numerov_module
