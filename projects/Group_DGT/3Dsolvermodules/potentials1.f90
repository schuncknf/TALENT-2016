module potentials

!potential function
use variables

contains 

!===============================================================
! Calculates the potential function
	! kpot= 0 no potential      
	! kpot= 1 square well of size a and deepth Vvalue
	! kpot= 2 Wood-Saxon potential
	! kpot= 3 Coulomb
	! kpot= 4 spin-orbit W potential
	! kpot= 5 Radial potential Eq. 2.10
	!
	! Variables:
	! vr= potential as a result
	! i= point in a mesh where potential is calculated
	! meshsize = distance between two meshpoints
	! Rmax = right bound for potential
	! Vvalue = depth of infinite potential 
	! numN = number of neutrons
	! numZ = number of protons
	! rzero = a constant related to a radius of nucleus
	! rProt = proton radius
	! charge = 0 (neutron), 1 (proton)
	! l = quantum number l
	! jj = quantum number j*2 (=>integer)
!===============================================================
         function potV(i,posiR,points,meshsize, Rmax, a,kpot, Vvalue, &
                  numN, numZ,rzero,rProt,charge, l, jj, &
                  pre_dens, pre_denZ, pre_denN) result(vr)
         
	 use variables
         implicit none
  !       INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
         real (kind = dp) :: meshsize, Rmax, Rmin,a,Vvalue, rzero, bigR, ferre
         real (kind=dp) :: vr, pos, leftlim, rightlim, pos0, rProt
         real (kind=dp) :: Uq_pot, Wq_pot1,Wq_pot2, V_Cou, V_cfug, Uq_sk
         real (kind=dp) :: rho, rho_q, rho_p, rho_n
         REAL(kind=dp), DIMENSION(0:points) :: pre_dens, pre_denN, pre_denZ
         integer :: i, kpot, numZ, numN, points, l, jj, charge, alphas
         real (kind=dp) :: t0s =-1132.400E0_dp
	 real (kind=dp) :: t3s = 23610.40E0_dp
	 real (kind=dp), dimension(0:points) :: posiR
	 alphas=1
! Skyrme parameters (now potential is for x0=x3=t1=x1=t2=x2=w=0)
! rho=total density 
! t0s = -1132.400
! alpha= 1
! t3s= 23610.40
! rho_q = neutron/proton density
! rho_n = neutron density
! rho_p = proton density

         bigR= rzero*(numN+numZ)**(1.0E0_dp/3.0E0_dp)
! initialize the density to the previous value of the loop
         rho = pre_dens(i)
	 if(charge .eq. 0) rho_q = pre_denN(i)
	 if(charge .eq. 1) rho_q = pre_denZ(i)

!-------------------------
! NO POTENTIAL

         if(i .eq. 0) then
           vr=0.0d0

         else if (i .ne. 0) then

!---------------------------
! SQUARE WELL
         if(kpot .eq. 0) then
             vr=0        
         else if(kpot .eq. 1) then 
           leftlim=6. !random numbers
           rightlim=10. !random numbers
           if(pos .ge. leftlim .and. pos .le. rightlim) then
             vr= Vvalue
           else
             vr= 0
           end if

!-------------------------------------
! WOODS-SAXON
           else if(kpot .eq. 2) then
             ferre=1/(1+exp((abs(posiR(i))-bigR)/a))
             vr=-44.0192307692*ferre
 !            vr=-(51-33*(numN-numZ)/(numN+numZ))*ferre

!---------------------------------------
! COULOMB
         else if(kpot .eq. 3) then
             if(abs(posiR(i)) .le. rProt) then
             vr=numZ*esquare/2/rProt*(3-(posiR(i)/rProt)**2)         
             else if(abs(posiR(i)) .gt. rProt) then
             vr=numZ*esquare/rProt 
             end if 

!-------------------------------------------
! SPIN-ORBIT
         else if(kpot .eq. 4) then
             ferre=-0.44*exp((abs(posiR(i))-bigR)/a) * rzero**2 &
                   /a/pos/(1+exp((abs(posiR(i))-bigR)/a))**2
             vr=-(51-33*(numN-numZ)/(numN+numZ))*ferre   

!--------------------------------------------  
! RADIAL POTENTIAL          TO DO CLEAN THE CONDITION ON THE POTENTIALS
         else if(kpot .eq. 5) then
             ! Uq potential is the Wood-Saxon
             if(i .eq. 0.1E0_dp)then
                pos = meshsize
             end if
             Uq_pot=-(51-33*(numN-numZ)/(numN+numZ)) &
                   /(1+exp((abs(posiR(i))-bigR)/a))
                     !- 44.0192307692 &
                    !/(1+exp((abs(posiR(i))-bigR)/a))
!                   
             !jj=2j so jj is integer
             ferre=-0.44*exp((abs(posiR(i))-bigR)/a) * rzero**2 &
                   /a/(posiR(i))/(1+exp((abs(posiR(i))-bigR)/a))**2
!
             Wq_pot1=ferre*(dble(jj)/2*(dble(jj)/2+1)-l*(l+1)-3/4)/(posiR(i)) 
! 
             if(i .eq. 0) Wq_pot1=0.0E0_dp     
             Wq_pot2=-(51-33*(numN-numZ)/(numN+numZ))
                     !-44.0192307692
                     !
!
             V_cfug=hb2m*((l*(l+1))/(posiR(i))**2)
             if(i .eq. 0) V_cfug=0.0E0_dp  
             vr=Uq_pot+Wq_pot1*Wq_pot2+V_cfug
! Coulomb part
             if(charge .eq. 1) then
               if(abs(posiR(i)) .le. rProt) then
                V_Cou=numZ*esquare/2/rProt*(3-(posiR(i)/rProt)**2)         
               else if(abs(posiR(i)) .gt. rProt) then
                V_Cou=numZ*esquare/rProt 
               end if
             vr= vr+V_Cou
             end if
!         end if


!--------------------------------------
! SIMPLE SKYRME (manuale 3.3.1) with t0 and t3
	else if(kpot .eq. 6) then
          Uq_sk = rho*t0s + (2+alphas)*t3s*(rho**alphas)/12 + rho_q*(-t0s/2 -t3s/12*(rho**alphas)) &
                   + alphas*(rho**(alphas-1))*(-t3s*(rho_p**2+rho_n**2)/24)
          V_cfug=hb2m*((l*(l+1))/(posiR(i))**2)
 !         if(i .eq. 0) V_cfug=0.0E0_dp 
          vr=-Uq_sk+V_cfug
! Coulomb part for Skyrme
!             if(charge .eq. 1) then
!               if(abs(posiR(i)) .le. rProt) then
!                V_Cou=numZ*esquare/2/rProt*(3-(posiR(i)/rProt)**2)         
!               else if(abs(posiR(i)) .gt. rProt) then
!                V_Cou=numZ*esquare/rProt 
!               end if
!	      end if

!	     vr=Uq_sk+V_Cou
        end if
       end if  
end function potV

end module potentials
