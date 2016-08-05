! Subroutine to look for the first phi and density trials.
! Ols WoodSaxon code used
subroutine WoodsSaxon3D(energy_n,energy_p,quant_n,quant_p,wf_n, wf_p)
    use globals
    implicit none
    integer, parameter :: dm=kind(1.d0)
      real (kind=dm) :: h, Em, nmfactor
      integer  :: NN, ZZ
      real (kind=dm),allocatable :: rho(:) ,Psi(:)! psi(grid_index)

      real (kind=dm) :: E_minus, E_plus=0., E_left, E_right

      integer(kind=dm) :: cnodes ! nodes of wavefuntion 
      integer(kind=dm) :: i,ii, uu, ipart
      integer(kind=dm) :: n, l
      real(kind=dm)    :: sm, j, x
      real(kind=dm) ::   wf_n(0:Nmesh, orbital), wf_p(0:Nmesh,orbital), energy_n(1:orbital), energy_p(1:orbital)
      real(kind=dm) ::   quant_n(1:orbital,1:4), quant_p(1:orbital,1:4)
      real(kind=dm),allocatable :: k_sq(:), Enl2j(:,:),dens_n(:), dens_t(:)
      real(kind=dm),allocatable :: all_wavefunction(:,:)
      real(kind=dm) :: vpot

      h=dx
      ZZ=proton
      NN=neutron

      E_minus=(-51. + 33.*(NN-ZZ)/(NN+ZZ))
      allocate(k_sq(0:Nmesh),psi(0:Nmesh), rho(0:Nmesh))
      allocate(Enl2j(1:orbital,1:4))
      allocate(all_wavefunction(0:Nmesh, 1:orbital))
      allocate(dens_n(0:Nmesh),dens_t(0:Nmesh) )
      
do ipart=1, 2 !1=neutrons, 2=protons
      Em=E_minus
      Enl2j(:, :) = 0.d0
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
                  do while (.true.) ! begin loop bissection method
                      if (abs(E_right-E_left)<epsi) exit
                      Em=(E_right+E_left)/2.
                      cnodes=0
                      !begin numerov method
                      Do i = 1,Nmesh
                              x=(i)*h     
                              k_sq(i) = (Em-Vpot(x,ii,ipart,l))/h2m -1.*l*(l+1)/x**2
                      Enddo
                               k_sq(0) = (Em-Vpot(0.d0,ii,ipart,l))/h2m 
  
                      Do i = 2,Nmesh
                              psi(i) = (2.*(1.-5.*h**2*k_sq(i-1)/12.)*psi(i-1) &
                              -(1.+h**2*k_sq(i-2)/12.)*psi(i-2))/(1.+h**2*k_sq(i)/12.)
                              if (psi(i-1)*psi(i)<0)  cnodes=cnodes+1 ! we find the number of nodes
                      Enddo
                      !end numerov method
                      !if((n.eq.0).and.(l.eq.7))  print*, cnodes, Em
                      if(cnodes.Gt.n) Then
                              E_right = Em
                      else 
                              E_left = Em
                      Endif
                  end do! end loop bissection method

                  Em=E_right !eigenvalue

                  if ((Em.gt.0.).or.(abs(Em-E_minus).lt.epsi).or.(abs(Em-E_plus).lt.epsi)) exit
                  
                  Enl2j(uu, 1)= Em
                  Enl2j(uu, 2)= n+1
                  Enl2j(uu, 3)= l
                  Enl2j(uu, 4)= sm
       
                  ! normalization of the eigenfunction
                  ! We recalculate psi because we want to find the function that does
                  ! not go to infinite SHOOTING METHOD WITH RIGHT HAND SIDE EQUAL TO 0

                  cnodes = 0
                      Do i = 1,Nmesh
                              x=(i)*h     
                              k_sq(i) = (Em-Vpot(x,ii,ipart,l))/h2m -1.*l*(l+1)/x**2
                      Enddo
                               k_sq(0) = (Em-Vpot(0.d0,ii,ipart,l))/h2m 

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
                 all_wavefunction(:,uu) = psi(:)
                 uu=uu+1
                 enddo!loop sm--->j
            enddo!loop l 
      enddo    !loop n
      
      call arrage_energy(Enl2j,all_wavefunction)
      
      if (ipart.eq.1) then
      quant_n(:,:) = Enl2j(:,:)
      energy_n(:) = Enl2j(:,1)
      wf_n(:,:)=all_wavefunction(:,:)
      else
      quant_p(:,:) =  Enl2j(:,:)
      energy_p(:) = Enl2j(:,1)
      wf_p(:,:)=all_wavefunction(:,:)
      endif
enddo !do particles

end subroutine WoodsSaxon3D

!=========================== potential function ===============================

! Woods-Saxon + SO potential
real (8) function Vpot(xy, mm, ip, l)
      use globals
      implicit none
      integer, parameter :: dm=kind(1.d0)
      integer  :: NN, ZZ
      real(kind=dm)::V0, R, fx, xy, ls, AA, VWS, VSO, VCO
      integer(kind=dm):: l
      INTEGER(8) :: mm, ip
      !WoodsSaxon
      NN=neutron
      ZZ=proton
      AA= NN+ZZ
      V0 = (-51. + 33.*(NN-ZZ)/AA)
      
      R= r0*AA**(1./3.)
      fx= 1./(1.+exp((xy-R)/a))
      !WS
      VWS= V0*fx
      !SO

      if (mm.eq. +1) then
         ls = 0.5*l      ! j=l+1/2
      else 
         ls = -0.5*(l+1) ! j=l-1/2
      end if
      VSO=0.44d0*r0**2*V0*exp((xy-R)/a)/(xy*a*(1.+exp((xy-R)/a))**2)*ls

      !CO
      if (xy.le.R) then
         VCO= e2*ZZ/2./R*(3.-((xy/R)**2))
      else
         VCO= e2*ZZ/xy
      endif
      !total
      if (ip==1) then
         if (xy.eq.0.)then
            Vpot= VWS 
         else  
            Vpot= VWS + VSO
         end if
      else if(ip==2) then
         if (xy.eq.0.)then
            Vpot= VWS + VCO
         else  
            Vpot= VWS + VSO + VCO
         endif   
      endif
end function

      
