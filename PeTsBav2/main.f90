      program infinite 
      use globals
      implicit none
      integer(kind=dm) :: numnodes ! the number of excited state, that has to be the same as number of nodes in wavefunction
      integer(kind=dm) :: cnodes ! nodes of wavefuntion 
      real(kind=dm),allocatable :: k_sq(:)
      allocate(k_sq(0:Nmesh),psi(0:Nmesh)
      real(kind=dm) :: x  
      numnodes= 0
      Em=E_minus
      psi(0) = 0.
      psi(1) = 0.1
      do while (comparison)!loop to compare if eigenvalues are between 0-100MeV.
            E_left=Em
            E_right= E_plus
            !begin bissection method
    
            do while (bisloop) ! begin loop bissection method
                  
                  Em=(E_right+E_left)/2.
                  cnodes=0     
                 !begin numerov method
                  Do i = 0,Nmesh
                   x=R_min+float(i)*h
                   k_sq(i) = (Em-vpot(x))/h2m
                  Enddo
                  Do i = 2,Nmesh
                   psi(i) = (2*(1.-5.*h**2*k_sq(i-1)/12.)*psi(i-1)-(1.+h**2*k_sq(i-2)/12.)*psi(i-2))/(1.+h**2*k_sq(i)/12.)
                  if (psi(i-1)*psi(i)<0)  cnodes=cnodes+1 ! we find the number of nodes
                  Enddo
      !end numerov method
                  if(cnodes.Gt.numnodes) Then !if nodes> number of excited state  E_right=Em
                     E_right = Em
                  else !(cnodesnumnodes) Then !otherwise E_left=Em
                     E_left = Em
                  Endif
                  if (E_right-E_left<eps) bisloop= .false.
      ! if Em-E_plus<epsilon comparison=false ->exit
            end do! end loop bissection method  
        !end bissection method    
            Em=E_right !eigenvalue
            if (Em-E_plus.lt. eps) comparison=.false.

      



      numnodes=numnodes+1
      enddo ! end loop to compare if eigenvalues are between 0-100MeV. 
      end program infinite



      function vpot(xy)
      use globals
      implicit none
      real(kind=dm)::xy
      if(xy.Lt.R_box) Then
         vpot = 0.
      else 
                
      endif
      end function


