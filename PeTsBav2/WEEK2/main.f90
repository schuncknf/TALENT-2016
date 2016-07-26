      ! Program to calculate the eigenvalues and 
      ! eigenvectors in a woods saxon potential 3D
      ! Parts: main.f90, globals.f90, Simpson.f90, WS.input
      !      
      !Authors: Fang Ni, Kai Wang, Claudia Gonzalez Boquera
      !TALENT COURSE 4 2016 
      program WoodsSaxon3D
      use globals  
      implicit none

      real(kind=dm) :: vpot    
      character(50) :: filename

      !------------------Read parameters from input---------------------
      NAMELIST / Nucleus / NN, ZZ
      NAMELIST / parameters / epsi, h2m
      NAMELIST / sizebox /  R_box, h 
      NAMELIST / energy / E_minus, E_plus
      NAMELIST / max_quantum_number / n_max, l_max
      NAMELIST / WoodSaxon / r0, a

      open (30, file='WS.input', status='old')
      read(30, NML= Nucleus)
      read(30, NML=parameters)
      read(30, NML=sizebox)
      read(30, NML=energy)
      read(30, NML=max_quantum_number)
      read(30, NML=WoodSaxon)
      
      close(30)
      !-----------------------------------------------------------------
      write(*,'(3a20)') 'nodes', 'eigenvalue', 'filename eigvec'

      Nmesh=nint((R_box)/h)
      print*, nmesh
      allocate(k_sq(0:Nmesh),psi(0:Nmesh), rho(0:Nmesh))

      numnodes= 0

      Em=E_minus

      sm=0.
      j=0.
      uu=1

orbital = 0
do n=0, n_max !loop n
      do l=0, l_max !loop l 
            do ii=1,-1, -2!loop sm--->j
                 if((l.eq.0).and.(ii.eq.-1)) exit 
                 orbital = orbital + 1            
            end do
      end do
end do 
write(*,*) orbital

allocate(Enl2j(2,orbital))

do n=0, n_max !loop n
      do l=0, l_max !loop l    
            do ii=1,-1, -2!loop sm--->j
                  sm=1./2.*ii
                  if((l.eq.0).and.(ii.eq.-1)) exit 
                  j=1.*l+sm
                        nl2j= 10000*(1+n)+100*l+int(2*j)
                        do while (comparison)!loop to compare if eigenvalues are greater than 0MeV. TODO
                              write(filename,'(I3.3,a4)') numnodes, '.dat'
                              filename = trim(filename)
                              E_left=E_minus
                              E_right= E_plus
                              psi(0) = 0.
                              psi(1) = 0.1
                              
                              !begin bissection method
                              do while (bisloop) ! begin loop bissection method

                                    if (abs(E_right-E_left)<1e-6) exit
                                    
                                    Em=(E_right+E_left)/2.
                                   ! print*, Em
                                    cnodes=0     
                                    !begin numerov method
                                    Do i = 0,Nmesh
                                          x=i*h
                                          k_sq(i) = (Em-Vpot(x, sm))/h2m -1.*(l*(l+1))/x**2 
                                         
                                    Enddo
                                    Do i = 2,Nmesh
                     psi(i) = (2.*(1.-5.*h**2*k_sq(i-1)/12.)*psi(i-1)-(1.+h**2*k_sq(i-2)/12.)*psi(i-2))/(1.+h**2*k_sq(i)/12.)
                  
                                          if (psi(i-1)*psi(i)<0)  cnodes=cnodes+1 ! we find the number of nodes
                                    Enddo
                                    !end numerov method

                  if(cnodes.Gt.n) Then 
                     E_right = Em
                  else 
                     E_left = Em
                  Endif

               !  if (abs(E_right-Em)<epsi) bisloop= .false.

            end do! end loop bissection method  
            Em=E_right !eigenvalue
            
            Enl2j(1, uu) = Em
            Enl2j(2, uu) = nl2j
            print*, Em, nl2j
            !normalization of the eigenfunction

            ! We recalculate psi because we want to find the function that does
            ! not go to infinite SHOOTING METHOD WITH RIGHT HAND SIDE EQUAL TO 0

            cnodes = 0
            Do i = 0,Nmesh
                  x=i*h
                  k_sq(i) = (Em-Vpot(x, sm))/h2m - 1.*((l*(l+1))/x**2) 
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
            nmfactor = 0d0
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
            write(*,'(3F10.8)' ) numnodes, Em, filename     
      enddo ! 
            uu=uu+1
            enddo!loop sm--->j
      enddo!loop l 

enddo    !loop n        


      
      end program WoodsSaxon3D

      ! Woods-Saxon potential
      real (8) function Vpot(xy, mm)
      use globals
      implicit none
      real(kind=dm)::V0, R, fx, xy, mm, ls, AA, VWS, VSO
      AA=NN+ZZ
      V0 = (-51. + 33.*(NN-ZZ)/AA)
      R= r0*AA**(1./3.)
      fx= 1./(1.+exp((xy-R)/a))
      VWS= V0*fx

        if (mm.eq. 0.5d0) then
           ls = 0.5*l      ! j=l+1/2
        else 
           ls = -0.5*(l+1) ! j=l-1/2
        end if
      VSO=0.44*r0**2*V0*exp((xy-R)/a)/(xy*a*(1.+exp((xy-R)/a))**2)*ls
      Vpot= VWS+VSO
      end function

      subroutine arrage_energy(matrix, numorb)
      use globals
      implicit none
      real(kind=dm) :: temp1, temp2
      real(kind=dm), allocatable :: matrix(:, :)
      integer(8) :: iii, jjj ,numorb
      allocate (matrix(1:2, numorb))
      do iii=1, numorb-1
            do jjj=1, numorb-iii
                  if (matrix(1, jjj).gt. matrix(1, jjj+1)) then
                  temp1= matrix(1,jjj)
                  temp2= matrix(2, jjj)
                  matrix(:,jjj) = matrix(:, jjj+1)
                  matrix(1, jjj+1)= temp1
                  matrix(2, jjj+1) = temp2
                  endif
            enddo
      enddo
      end subroutine 


