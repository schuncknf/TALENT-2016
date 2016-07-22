      ! Program to calculate the eigenvalues and 
      ! eigenvectors in a finite square well potential
      ! Parts: main.f90, globals.f90, Simpson.f90, finite.input
      !      
      !Authors: Fang Ni, Kai Wang, Claudia Gonzalez Boquera
      !TALENT COURSE 4 2016 
      program WoodsSaxon3D
      use globals  
      implicit none

      real(kind=dm) :: VWS,VSO,vpot    
      character(50) :: filename

      !------------------Read parameters from input---------------------
      NAMELIST / parameters / epsi, h2m
      NAMELIST / sizebox / R_max, R_min, R_box, h , a  
      NAMELIST / energy / E_minus, E_plus, V0

      open (30, file='WS.input', status='old')
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


      n=0
      l=0
      sm=0.
      j=0.


do while (comparison1)!loop n
      do l=0, n  !loop l    
            do ii=-1,1, 2!loop sm--->j
                  sm=1./2.*ii
                  j=float(l)+sm
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
                                          vpot= VWS(x)+VSO(x, j, l, sm)      
                                          k_sq(i) = (Em-vpot)/h2m - float(l*(l+1))/x**2 
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
      
            write(*,'(3F10.8)' ) numnodes, Em, filename
           
            numnodes=numnodes+1

      enddo ! end loop to compare if eigenvalues are between 0-100MeV. 






            enddo!loop sm--->j
      enddo!loop l 

enddo    !loop n        













      
      end program WoodsSaxon3D

      ! Woods-Saxon potential
      real (8) function VWS(xy)
      use globals
      implicit none
      real(kind=dm)::xy, V0, R, fx
      V0 = (-51. + 33.*(N-Z)/A)
      R= r0*A**(1./3.)
      fx= 1./(1.+exp((xy-R)/a))
      VWS= V0*fx
      end function

      real (8) function VSO(xy, j, l, sm)
      use globals
      implicit none
      real(kind=dm)::xy, VLS, r0, xy, dfxdr
      VLS       

      VSO= VLS*r0**2./xy*dfxdr
      end function
















