      ! Program to calculate the eigenvalues and 
      ! eigenvectors in a woods saxon potential 3D
      ! Parts: main.f90, globals.f90, Simpson.f90, WS.input
      !      
      ! Authors: Fang Ni, Kai Wang, Claudia Gonzalez Boquera
      ! TALENT COURSE 4 2016 
      program WoodsSaxon3D
      use globals  
      implicit none

      real(kind=dm) :: vpot
      integer(kind=dm):: Nn_l, Nn_tmp
      integer(kind=dm):: iorb,orbital_tmp

      !------------------Read parameters from input---------------------
      NAMELIST / Nucleus / NN, ZZ
      NAMELIST / parameters / epsi, h2m, e
      NAMELIST / sizebox /  R_box, h 
      NAMELIST / energy / E_minus, E_plus
      NAMELIST / max_quantum_number / n_max, l_max
      NAMELIST / WoodSaxon / r0, a

      open(30, file='WS.input', status='old')
      read(30, NML= Nucleus)
      read(30, NML=parameters)
      read(30, NML=sizebox)
      read(30, NML=energy)
      read(30, NML=max_quantum_number)
      read(30, NML=WoodSaxon)   
      close(30)
      !-----------------------------------------------------------------
      Nmesh=nint((R_box)/h)
      allocate(k_sq(0:Nmesh),psi(0:Nmesh), rho(0:Nmesh))
      orbital = 0
      do n=0, n_max !loop n
            do l=0, l_max !loop l 
                  do ii=1,-1, -2!loop sm--->j
                  if((l.eq.0).and.(ii.eq.-1)) exit 
                        orbital = orbital + 1            
                  end do
            end do
      end do 
      write(*,*) 'number of calculated orbitals',orbital
      allocate(Enl2j(1:orbital,1:4),filename_wave_func(1:orbital))
      allocate(all_wavefunction(1:orbital,0:Nmesh))
      allocate(dens_n(0:Nmesh),dens_t(0:Nmesh) )
      dens_t(:)=0.
do ipart=1, 2 !1=neutrons, 2=protons
      Em=E_minus
      uu=1

      do n=0, n_max-1 !loop n
            do l=0, l_max  !loop l   
                  do ii=1,-1, -2 !loop sm--->j
                  sm=1./2.*ii
                  if((l.eq.0).and.(ii.eq.-1)) exit 
                  j=l+sm                        
                  E_left=E_minus
                  E_right= E_plus
                  psi(0) = 0.
                  psi(1) = 0.1
                  !begin bissection method
                  do while (bisloop) ! begin loop bissection method
                      if (abs(E_right-E_left)<1e-6) exit
                      Em=(E_right+E_left)/2.                  
                      cnodes=0     
                      !begin numerov method
                      Do i = 1,Nmesh-1
                              x=(i)*h
                              k_sq(i) = (Em-Vpot(x, ii, ipart))/h2m -1.*(l*(l+1))/x**2
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
                  Em=E_right !eigenvalue
      
                  Enl2j(uu, 1)= Em
                  Enl2j(uu, 2)= n+1
                  Enl2j(uu, 3)= l
                  Enl2j(uu, 4)= sm

                  ! normalization of the eigenfunction
                  ! We recalculate psi because we want to find the function that does
                  ! not go to infinite SHOOTING METHOD WITH RIGHT HAND SIDE EQUAL TO 0

                  cnodes = 0
                  Do i = 1,Nmesh-1
                        x=(i)*h
                        k_sq(i) = (Em-Vpot(x, ii, ipart))/h2m - 1.*((l*(l+1))/x**2) 
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
                 rho(:) = abs(psi(:))**2
                 call simpson(Nmesh, h, rho, nmfactor)       
                 psi(:) = psi(:)/sqrt(nmfactor)    !wave func
                 
                 if (ipart.eq.1) then 
                 write(filename_wave_func(uu),'(a8,a1,I3.3,a1,I3.3,a4,2i2.2,i3.3,a4)') &
                  'Neutron_','p', ZZ, 'n', NN, 'nl2j', n+1, l, nint(2*j), '.dat'
                 else
                  write(filename_wave_func(uu),'(a7,a1,I3.3,a1,I3.3,a4,2i2.2,i3.3,a4)') &
                  'Proton_','p', ZZ, 'n', NN, 'nl2j', n+1, l, nint(2*j), '.dat'
                 endif                
                 open(10,file=filename_wave_func(uu))
                 write(10,'(a5,3i5)') '#nl2j=', n+1, l, nint(2*j)
           
                 do i=1,Nmesh-1
                        x = i*h
                        write(10,'(3f15.8)') x, psi(i), vpot(x, ii, ipart)
                 end do
                 close(10)
                 all_wavefunction(uu,:) = psi(:)
                 uu=uu+1
                 enddo!loop sm--->j
            enddo!loop l 
      enddo    !loop n  

      call arrage_energy(Enl2j, orbital, all_wavefunction)
      if (ipart.eq.1)then
      open(111,file='energy_level_n.output')
      write(111,*) '      energy      n    l         spin'
      do i=1, orbital
            write(111,100) Enl2j(i,1),nint(Enl2j(i,2)),nint(Enl2j(i,3)),Enl2j(i,4) 
      enddo
      else
      open(111,file='energy_level_p.output')
      write(111,*) '      energy      n    l         spin'
      do i=1, orbital
            write(111,100) Enl2j(i,1),nint(Enl2j(i,2)),nint(Enl2j(i,3)),Enl2j(i,4) 
      enddo
      endif

!DENSITY PROFILE

if(ipart.eq.1) then
open(12,file='density_N.dat')
Nn_tmp = NN
else
open(12,file='density_P.dat')
Nn_tmp = ZZ
endif 
     
dens_n = 0d0

do iorb=1,orbital
   Nn_l = Nn_tmp
   j=Enl2j(iorb, 3) + Enl2j(iorb,4) 
   Nn_tmp = Nn_tmp - nint(2d0*j+1.)
   if (Nn_tmp<=0) then
      orbital_tmp = iorb-1
      exit
   end if
   dens_n(:) = dens_n(:) + (2d0*j+1) * abs(all_wavefunction(iorb,:))**2/(4d0*pi)
end do
dens_n(:) = dens_n(:) + Nn_l*abs(all_wavefunction(orbital_tmp+1,:))**2/(4d0*pi)

   dens_t(:) = dens_t(:) + dens_n(:)
do i=1,Nmesh-1
   x = i*h
   write(12,'(2f15.8)') x, dens_n(i)/x**2 
end do
close(12)

enddo !do particles
open(13,file='density_T.dat')
do i=1,Nmesh-1
   x = i*h
   write(13,'(2f15.8)') x, dens_t(i)/x**2     
end do
close(13)

100   format(1f15.8, 2I5.2, 1f15.3)

      end program WoodsSaxon3D



!=========================== potential function ===============================


      ! Woods-Saxon +SO potential
      real (8) function Vpot(xy, mm, ip)
      use globals
      implicit none
      real(kind=dm)::V0, R, fx, xy, ls, AA, VWS, VSO, VCO, Rp
      INTEGER(8) :: mm, ip
      !WoodsSaxon
      AA=NN+ZZ
      V0 = (-51. + 33.*(NN-ZZ)/AA)
      R= r0*AA**(1./3.)
      fx= 1./(1.+exp((xy-R)/a))
      VWS= V0*fx
      !SO
      if (xy.eq.0.)then
            Vpot= VWS    
      else  
            if (l.eq.0) then
            Vpot= VWS
            else
                  if (mm.eq. +1) then
                        ls = 0.5*l      ! j=l+1/2
                  else 
                        ls = -0.5*(l+1) ! j=l-1/2
                  end if
                  VSO=0.44*r0**2*V0*exp((xy-R)/a)/(xy*a*(1.+exp((xy-R)/a))**2)*ls
                  Vpot= VWS+VSO
            endif
      !Coulomb
                  if (ip.eq.2) then 
          !        Rp= r0*ZZ**(1./3.)
                  Rp=R
                        if (xy.le.Rp) then
                              VCO= e**2*ZZ/2./Rp*(3.-((xy/Rp)**2))
                        else
                              VCO= e**2*ZZ/xy
                        endif
                  Vpot=Vpot +VCO
                  endif

      endif
      end function



