      program infinite 
      use globals
      implicit none
      real(kind=dm) :: vpot
      character(50) :: filename
     
      NAMELIST / parameters / epsi, h2m
      NAMELIST / sizebox / R_max, R_min,  R_box, h   
      NAMELIST / energy / E_minus, E_plus
      open (30, file='infinite.input', status='old')
      read(30, NML=parameters)
      read(30, NML=sizebox)
      read(30, NML=energy)
      close(30)

      Nmesh=nint((R_max-R_min)/h)
      allocate(k_sq(0:Nmesh),psi(0:Nmesh), rho(0:Nmesh))
      numnodes= 0
      Em=E_minus
      do while (comparison)!loop to compare if eigenvalues are between 0-100MeV.
      write(filename,'(I3.3,a4)') numnodes, '.dat'
      filename = trim(filename)
            E_left=Em
            E_right= E_plus
            psi(0) = 0.
            psi(1) = 0.1
            !begin bissection method
            do while (bisloop) ! begin loop bissection method

                  if (E_right-E_left<1e-6) exit
                  Em=(E_right+E_left)/2.
                  cnodes=0     
                 !begin numerov method
                  Do i = 0,Nmesh
                   x=R_min+i*h
                   k_sq(i) = (Em-vpot(x))/h2m
                  
                  Enddo
                  Do i = 2,Nmesh
                   psi(i) = (2.*(1.-5.*h**2*k_sq(i-1)/12.)*psi(i-1)-(1.+h**2*k_sq(i-2)/12.)*psi(i-2))/(1.+h**2*k_sq(i)/12.)
                  if (psi(i-1)*psi(i)<0)  cnodes=cnodes+1 ! we find the number of nodes
                  Enddo
      !end numerov method
!                  print*, cnodes, Em, E_left, E_right
                  if(cnodes.Gt.numnodes) Then !if nodes> number of excited state  E_right=Em
                     E_right = Em
                  else !(cnodesnumnodes) Then !otherwise E_left=Em
                     E_left = Em
                  Endif
 !                 print*, cnodes, numnodes
                 if ((E_plus-E_left)<epsi) bisloop= .false.
      ! if Em-E_plus<epsilon comparison=false ->exit
            end do! end loop bissection method  
        !end bissection method    
            Em=E_right !eigenvalue


!normalization of the eigenfunction

! We recalculate psi because we want to find the function that does
! not go to infinite SHOOTING METHOD WITH RIGHT HAND SIDE EQUAL TO 0


      cnodes = 0
      Do i = 0,Nmesh
       x=R_min+i*h
       k_sq(i) = (Em-vpot(x))/h2m  
      Enddo
            psi(0) = 0.
            psi(1) = 0.1
      Do i = 2,Nmesh
      psi(i) = (2.*(1.-5.*h**2*k_sq(i-1)/12.)*psi(i-1)-(1.+h**2*k_sq(i-2)/12.)*psi(i-2))/(1.+h**2*k_sq(i)/12.)
      if (psi(i-1)*psi(i)<0)  cnodes=cnodes+1 ! we find the number of nodes
        
      if (psi(i)*psi(i-1)<=0d0.and.cnodes>numnodes) then
      psi(i) = 0.
      end if
      
      end do
    
! normalization
     rho(:) = abs(psi(:))**2
      call simpson(Nmesh, h, psi, nmfactor)
      nmfactor = sum(rho(:))*h   
      psi(:) = psi(:)/sqrt(nmfactor)    !wave func
       

      if (E_plus-Em.lt. epsi) exit
      open(10,file=filename)
      do i=0,Nmesh
         x = R_min + i*h
         write(10,*) x, psi(i)
      end do
      close(10)
           
      numnodes=numnodes+1
      enddo ! end loop to compare if eigenvalues are between 0-100MeV. 
      end program infinite


      real (kind=dm) function vpot(xy)
      use globals
      implicit none
      real(kind=dm)::xy
      vpot=0. 
          
      end function

!***************************************************
!*      this subroutine is simpson's formula       *
!*      N: the number of division                  *
!*      dx: axis step size                         *
!*      f: array of function                       *
!*      int: integrated value                      *
!***************************************************

      subroutine simpson(N, dx, f, res)
      implicit none
      real(8) :: res, even, odd, x, dx
      integer :: N
      real(8) :: f(0:N)
      integer :: i
      res = 0d0
      odd = 0d0
      even = 0d0
      do i=1,N, 2
      odd = odd + f(i)
      end do
      odd = 4.*odd
      do i=2,N, 2
      even = even + f(i)
      end do
      even = 2*even
      res = f(0) + f(N) + odd + even
      res = dx/3*res
      end subroutine simpson



