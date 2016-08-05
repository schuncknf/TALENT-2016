module wavefunctions

contains

!============================================================================
! CALCULATION OF WAVEFUNCTION
!  
! points = maximum number of meshpoints
! meshsize = a distance between two meshpoints
! n_rad = radial quantum number
! l = quantum number l
! Eup0 = inital upper boundary for energy
! Edown0 = initial lower boundary for energy
! epsil = demanded precision
! pott = 
! trialU = trial wavefunction
! Etrial = energy of trial wavefunction
!=========================================================
!wave function calculation
         subroutine wavef(points,meshsize, n_rad,l, Eup0, Edown0, epsil, pott, &
                 trialU, Etrial) 
         use variables
         implicit none
         REAL(kind=dp), INTENT(IN) :: Eup0, Edown0, epsil
         REAL(kind=dp) :: Eup, Edown
         REAL(kind=dp), INTENT(IN)  :: meshsize
         integer, INTENT(IN)  :: n_rad, points, l
         REAL(kind=dp), INTENT(OUT) :: Etrial
         integer:: converg, kk
         integer :: Nodecount , ifail, i, idir, imin, imax, kpot
         REAL(kind=dp) :: a1, a2, a3, normal, Eexp, rcentr
         REAL(kind=dp), DIMENSION(0:points), INTENT(OUT) :: trialU
         REAL(kind=dp) :: pot, a, Vvalue
         REAL(kind=dp) , DIMENSION(0:points), INTENT(IN) :: pott
!
         !loading the potential
           Eup=Eup0
           Edown=Edown0
!
! If l is greater than 5, we start calculating the wavefunction from imin (not 
!from zero because of central potential)
           if (l .gt. 5) then
             rcentr=sqrt(Eup/(l*(l+1))/hb2m)
             imin=int(rcentr/meshsize)
           else
             imin=1
           end if
!
           idir=1
           do i=0,imin-1
              trialU(i)=0.0E0_dp!meshsize**(l+1)
           end do

!           if(l .eq. 0) then
!           trialU(0)=meshsize**(l+1)
!           end if
           trialU(imin)=meshsize**(l+1)
           imax=points-1

!loop until convergence
         converg=0
         do while (converg .ne. 1)
         Etrial=(Eup+Edown)/2
!         print*, Etrial

!building up the wave function up to Rmax

         do i=imin,imax,idir
           pot=(Etrial-pott(i))/hb2m !insert V if V=0
           a1= 2.0E0_dp*(1.0E0_dp-5.0E0_dp/12.0E0_dp*pot*meshsize**2)
           a2=(1.0E0_dp+1.0E0_dp/12.0E0_dp*pot*meshsize**2)
           a3=(1.0E0_dp+1.0E0_dp/12.0E0_dp*pot*meshsize**2)
           trialU(i+idir)= (a1*trialU(i)-a2*trialU(i-idir))/a3
         end do

! counting nodes
         Nodecount=0    
         do i=1,points
           if(trialU(i-1)*trialU(i) .lt. 0.0E0_dp) &
             Nodecount= Nodecount+1
         end do 
         if(trialU(points-1)*trialU(points) .lt. 0.0E0_dp) &
             Nodecount= Nodecount-1 
!         print*, Nodecount

!checking nodes    
         if(Nodecount .gt. n_rad) then
           Eup=Etrial
         else if(Nodecount .le. n_rad) then
           Edown=Etrial
         end if
! checking convergence
         if (abs(Eup-Edown) .lt. epsil) then
           converg=1
         end if
41       continue
        end do
 
         end subroutine wavef


!======================================================================
! CALCULATING THE WAVEFUNCTION IN TWO PARTS: LEFT AND RIGHT PART
! points = maximum number of meshpoints
! meshsize = a distance between two meshpoints
! n_rad = radial quantum number
! l = quantum number l
! numN = number of neutrons
! numZ = number of protons
! Eup0 = inital upper boundary for energy
! Edown0 = initial lower boundary for energy
! epsil = demanded precision
! pott = 
! trialU = trial wavefunction
! Etrial = energy of trial wavefunction
!============================================================
! 

         subroutine parwf(points,meshsize, n_rad, l, numN, numZ, Etrial, &
                    epsil, pott, trialU) 
         use variables
         implicit none
         REAL(kind=dp), INTENT(IN) :: Etrial, epsil
         REAL(kind=dp), INTENT(IN)  :: meshsize
         integer, INTENT(IN)  :: n_rad, l, numN, numZ
!         REAL(kind=dp):: Eup, Edown
         integer, INTENT(IN):: points
         integer :: Nodecount , ifail, i, idir, imin, imax, converg, i_sur
         REAL(kind=dp) :: a1, a2, a3, normal, normalr, normall, rcentr
         REAL(kind=dp), DIMENSION(0:points), INTENT(OUT) :: trialU
         REAL(kind=dp), DIMENSION(0:points) ::  triall, trialr
         REAL(kind=dp) :: pot, surf
         REAL(kind=dp), DIMENSION(0:points), INTENT(IN) :: pott
         REAL(kind=dp) :: derivr, derivl, contir, contil, coeff
!
! half left wavefunction

            imin=1
	    idir=1

!           do i=0,imin-1
              triall(0)=0.0E0_dp!meshsize**(l+1)
!           end do
           triall(imin)=meshsize**(l+1)
           imax=points-1!because f(x+2h) is necessary
         do i=imin,imax,idir
           pot=(Etrial-pott(i))/hb2m !insert V if V=0
           a1= 2.0E0_dp*(1.0E0_dp-5.0E0_dp/12.0E0_dp*pot*meshsize**2)
           a2=(1.0E0_dp+1.0E0_dp/12.0E0_dp*pot*meshsize**2)
           a3=(1.0E0_dp+1.0E0_dp/12.0E0_dp*pot*meshsize**2)
           triall(i+idir)= (a1*triall(i)-a2*triall(i-idir))/a3
         end do
!
! half right wavefunction
           imin=1
           idir=-1

!           do i=0,imin-1
              trialr(0)=0.0E0_dp!meshsize**(l+1)
!           end do

           imax=points-1
           trialr(points)=0.0E0_dp
           trialr(points-1)= (100.**n_rad)*1.0E-8_dp

         do i=imax,imin+1,idir
           pot=(Etrial-pott(i))/hb2m !insert V if V=0
           a1= 2.0E0_dp*(1.0E0_dp-5.0E0_dp/12.0E0_dp*pot*meshsize**2)
           a2=(1.0E0_dp+1.0E0_dp/12.0E0_dp*pot*meshsize**2)
           a3=(1.0E0_dp+1.0E0_dp/12.0E0_dp*pot*meshsize**2)
           trialr(i+idir)= (a1*trialr(i)-a2*trialr(i-idir))/a3
         end do

! calculating derivatives and wave function values at "surface" point
           surf=1.25E0_dp*(numN+numZ)**(1.0E0_dp/3.0E0_dp)
!           print*, surf
           i_sur=int(surf/meshsize)
           derivl= (triall(i_sur+2)/12-2*triall(i_sur+1)+2*triall(i_sur-1)/3 &
                 -triall(i_sur-2)/12)/4/meshsize
           contil=triall(i_sur)
           derivr= (trialr(i_sur+2)/12-2*trialr(i_sur+1)+2*trialr(i_sur-1)/3 &
                 -trialr(i_sur-2)/12)/4/meshsize
           contir=trialr(i_sur)
! 
! matching condition
!
         converg=0
         do while(converg .ne. 1)

!
!
            !if(derivr .lt. 1.0E-10_dp) then
!              print*, "left function =", contil
!              print*, "right function =", contir
!              print*, "left derivative =", derivl
!              print*, "right derivative =", derivr
          !  end if

            coeff= contil/contir
!              print*, "coefficient =", coeff

            trialr(:)= coeff*trialr(:)
!            triall(:)=derivr*triall(:)
!            trialr(:)=derivl*trialr(:)
!

!
!check values
           derivl= (triall(i_sur+2)/12-2*triall(i_sur+1)+2*triall(i_sur-1)/3 &
                 -triall(i_sur-2)/12)/4/meshsize
           contil=triall(i_sur)
           derivr= (trialr(i_sur+2)/12-2*trialr(i_sur+1)+2*trialr(i_sur-1)/3 &
                 -trialr(i_sur-2)/12)/4/meshsize
           contir=trialr(i_sur)
!           
!           if((abs(contir-contil)/abs(contir+contil)) .lt. 1.0E-6_dp .and.  &
!             (abs(derivr-derivl)/abs(derivr+derivl)) .lt. 1.0E-6_dp) then 
           if((abs(contir-contil)/abs(contir+contil)) .lt. 1.0E-6_dp) then
 !            print '("Merging successfull")'
!             print *, "continuity: ", contil, contir
!             print *, "derivative: ", derivl, derivr
!             print*, "normal L=", normall
!             print*, "normal R=",normalr
             converg=1
           else      
             print '("Error in merging left and right wavefunctions")'
!             print *, "continuity: ", contil, contir
!             print *, "derivative: ", derivl, derivr
           end if
!
          end do
!
!         allocate(trialU(0:points), stat = ifail)
!         if (ifail .ne. 0) STOP           
!
           do i=0,i_sur
              trialU(i)= triall(i)
           end do
           do i=i_sur+1, points
              trialU(i)= trialr(i)
           end do
         end subroutine parwf


end module wavefunctions
