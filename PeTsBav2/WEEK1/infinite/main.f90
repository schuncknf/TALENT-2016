      ! Program to calculate the eigenvalues and 
      ! eigenvectors in a infinite well potential
      ! Parts: main.f90, globals.f90, Simpson.f90, infinite.input
      !      
      !Authors: Fang Ni, Kai Wang, Claudia Gonzalez Boquera
      !TALENT COURSE 4 2016 
      program infinite 
      use globals
      implicit none

      real(kind=dm) :: vpot
      character(50) :: filename

      !------------------Read parameters from input---------------------   
      NAMELIST / parameters / epsi, h2m
      NAMELIST / sizebox / R_max, R_min,  R_box, h   
      NAMELIST / energy / E_minus, E_plus, V0
      open (30, file='infinite.input', status='old')
      read(30, NML=parameters)
      read(30, NML=sizebox)
      read(30, NML=energy)
      close(30)
      !-----------------------------------------------------------------

      write(*,'(3a20)') 'nodes', 'eigenvalue', 'filename eigvec'

      Nmesh=nint((R_box-R_min)/h)

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

                  if (E_right-E_left<epsi) exit
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

                  if(cnodes.Gt.numnodes) Then 
                     E_right = Em
                  else 
                     E_left = Em
                  Endif

                 if ((E_plus-E_left)<epsi) bisloop= .false.

            end do! end loop bissection method  
 
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
            call simpson(Nmesh, h, rho, nmfactor)
            psi(:) = psi(:)/sqrt(nmfactor)    !wave func
     
            if (E_plus-Em.lt. epsi) exit
            open(10,file=filename)
            do i=0,Nmesh
                  x = R_min + i*h
                  write(10,*) x, psi(i)
            end do
            close(10)

            write(*, *) numnodes, Em, filename
           
            numnodes=numnodes+1

      enddo ! end loop to compare if eigenvalues are between 0-100MeV. 
      end program infinite


      real (kind=dm) function vpot(xy)
      use globals
      implicit none
      real(kind=dm)::xy
      vpot=0.d0
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
      do i=1,N-1, 2
      odd = odd + f(i)
      end do
      odd = 4.*odd
      do i=2,N-2, 2
      even = even + f(i)
      end do
      even = 2.*even
      res = f(0) + f(N) + odd + even
      res = dx/3.*res
      end subroutine simpson



