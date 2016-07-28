     program numerov
      use variables
      implicit none
!     program to solve 1D Schroedinger equation
!      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
      REAL(kind=dp) :: Eup, Edown, epsil
      REAL(kind=dp) :: Rmin, Rmax, meshsize, a, Vvalue, rzero
      integer :: Nodemax, kpot
      REAL(kind=dp) :: Etrial, Eexp
      integer :: points, ifail, converg, i, charge, l, jj, iii, jmin, jmax
      integer :: Nodecount, numN, numZ, strval, n_rad, n_prin
      REAL(kind=dp) :: a1, a2, a3, pot, normal, deltae
      REAL(kind=dp) :: potV, rProt
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: trialwf
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: pott
      CHARACTER(LEN=2):: Nodecc
      CHARACTER(LEN=5) :: string
!==========     set initial values     ==============================
      Eup=1000.E0_dp
      Edown=-100.E0_dp
      epsil=1.E-8_dp
!      Nodemax=5
      Rmin=0.E0_dp
      Rmax=30.E0_dp !fermi
      meshsize=0.1E0_dp
! kpot= 0 no potential      
! kpot= 1 square well of size a and deepth Vvalue
! kpot= 2 Wood-Saxon potential
! kpot= 3 Coulomb
! kpot= 4 spin-orbit W potential
! kpot= 5 Radial potential Eq. 2.10
      numN= 126
      numZ= 82 
      kpot= 5
      a=0.67E0_dp
      Vvalue=-100.E0_dp 
      rzero=1.27E0_dp 

      rProt= rzero*(numZ**(1.0E0_dp/3.0E0_dp))
!     single particle inputs
!      n_prin= n_rad+l+1
      n_prin=4
!==================================================================     
!
      points=(Rmax-Rmin)/meshsize
         allocate(trialwf(0:points), stat = ifail)
         if (ifail .ne. 0) STOP

!         trialwf(:)=0.0E0_dp
         allocate(pott(0:points), stat = ifail)
         if (ifail .ne. 0) STOP

      iii=0
      OPEN(UNIT=8, FILE='summary_table.dat', status='unknown')
      do charge=0,0!1
         do n_rad= 0,n_prin
            do l= 0, n_prin
               if((2*n_rad+l) .gt. n_prin) then
                go to 31
               end if 
 !               if(l .eq. 0) then
 !                jmin=1
 !                jmax=1
 !               else
 !                jmin=2*l-1
 !                jmax=2*l+1
 !               end if
               jmin=abs(2*l-1)
               do jj= jmin,2*l+1,2
               iii=iii+1
               print*,iii

!         OPEN(UNIT=9, FILE='potential.dat', status='unknown')
         do i=0,points
         pott(i)=potV(i,meshsize,Rmax,Rmin,a,kpot, Vvalue, &
                  numN, numZ,rzero,rProt,charge, l, jj)
         write(9,*) Rmin+i*meshsize, pott(i)
         end do
!         CLOSE(unit=9)
!

           call wavef(points,meshsize, n_rad, l, Eup, Edown, epsil, pott, &
                 trialwf, Etrial) 

           call parwf(points,meshsize, n_rad, l, numN, numZ, Etrial, &
                    epsil, pott, trialwf) 
!renormalization

         normal=0.0E0_dp
         do i=0,points
            normal = normal + trialwf(i)**2*meshsize
         end do
         trialwf(:)=trialwf(:)/sqrt(normal)
!
!         print '("Wavefunction normalized")'
!         print*, normal
! output number of nodes excluded a possible node in the last point of the box
! normalized wavefunction
         Nodecount=0    
         do i=1,points
           if(trialwf(i-1)*trialwf(i) .lt. 0.0E0_dp) &
             Nodecount= Nodecount+1
         end do 
         if(trialwf(points-1)*trialwf(points) .lt. 0.0E0_dp) &
             Nodecount= Nodecount-1 
         print*, "Number of nodes= ", Nodecount
         print*, "Energy= ", Etrial
!
         if ( iii < 10 ) then
           write( Nodecc, '("0",i1)' ) iii
         else
           write( Nodecc, '(i2)' ) iii
         end if

         OPEN(UNIT=7, FILE='wf'//Nodecc//'.dat', status='unknown')
         
         write(7,*) "Number of nodes=", Nodecount
         write(7,*) "Energy= ", Etrial
         do i=0,points
           write(7,*) Rmin+i*meshsize, trialwf(i)
         end do
         CLOSE(unit=7)

         print*, "------------------------------------- "
       write(8,'("n_rad",i3,2x,"l",i3,2x,"2j",i3,4x,"wf",i4,4x,"Energy ",1f18.10)') &
                n_rad, l, jj, iii, Etrial
!

                end do !loop j
              end do ! loop l
31            continue
            end do !  loop n_rad
         end do ! charge             
         CLOSE(unit=8)
         deallocate(trialwf)
         deallocate(pott)
      end program numerov
! end
!
!
!
 

!wave function calculation
         subroutine wavef(points,meshsize, n_rad,l, Eup0, Edown0, epsil, pott, &
                 trialwf, Etrial) 
         use variables
         implicit none
         REAL(kind=dp), INTENT(IN) :: Eup0, Edown0, epsil
         REAL(kind=dp) :: Eup, Edown
         REAL(kind=dp), INTENT(IN)  :: meshsize
         integer, INTENT(IN)  :: n_rad, points, l
         REAL(kind=dp), INTENT(OUT) :: Etrial
         integer:: converg
         integer :: Nodecount , ifail, i, idir, imin, imax, kpot
         REAL(kind=dp) :: a1, a2, a3, normal, Eexp
         REAL(kind=dp), DIMENSION(0:points), INTENT(OUT) :: trialwf
         REAL(kind=dp) :: pot, a, Vvalue
         REAL(kind=dp) , DIMENSION(0:points), INTENT(IN) :: pott
!

!         allocate(trialwf(0:points), stat = ifail)
!         if (ifail .ne. 0) then
!          print*, "failed allocation"
!          STOP
!         end if
         !loading the potential
           Eup=Eup0
           Edown=Edown0
!
           idir= 1
           trialwf(0)=0.0E0_dp
           trialwf(1)= meshsize*(l+1)
           imin=1
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
           trialwf(i+idir)= (a1*trialwf(i)-a2*trialwf(i-idir))/a3
         end do
! counting nodes
         Nodecount=0    
         do i=1,points
           if(trialwf(i-1)*trialwf(i) .lt. 0.0E0_dp) &
             Nodecount= Nodecount+1
         end do 
!         print*, Nodecount
!checking nodes    
         if(Nodecount .gt. n_rad) then
           Eup=Etrial
         else if(Nodecount .le. n_rad) then
           Edown=Etrial
         end if
! checking convergence
         if (abs(Eup-Edown) .lt. epsil) then
 !          print '("Convergence obtained for wavefunction")'
!           print *, "with energy= ", Etrial
! expected energy from the formula Eq. 1.14 in manuale_hf.pdf
!          if (kpot .eq. 0) then 
!            Eexp= ((Node+1)*pi)**2*hb2m/(Rmax-Rmin)**2
!            print *, "Expected energy= ", Eexp
!          end if
           converg=1
         end if
        end do
 !       do i=0,points
 !            print*, trialwf(i)
 !          end do
         end subroutine wavef
!
!
!
         subroutine parwf(points,meshsize, n_rad, l, numN, numZ, Etrial, &
                    epsil, pott, trialwf) 
         use variables
         implicit none
         REAL(kind=dp), INTENT(IN) :: Etrial, epsil
         REAL(kind=dp), INTENT(IN)  :: meshsize
         integer, INTENT(IN)  :: n_rad, l, numN, numZ
!         REAL(kind=dp):: Eup, Edown
         integer, INTENT(IN):: points
         integer :: Nodecount , ifail, i, idir, imin, imax, converg, i_sur
         REAL(kind=dp) :: a1, a2, a3, normal, normalr, normall
         REAL(kind=dp), DIMENSION(0:points), INTENT(OUT) :: trialwf
         REAL(kind=dp), DIMENSION(0:points) ::  triall, trialr
         REAL(kind=dp) :: pot, surf
         REAL(kind=dp), DIMENSION(0:points), INTENT(IN) :: pott
         REAL(kind=dp) :: derivr, derivl, contir, contil, coeff
!
! half left wavefunction
!         allocate(triall(0:points), stat = ifail)
!         if (ifail .ne. 0) STOP 
!           
           idir=1
           triall(0)=0.0E0_dp
           triall(1)= meshsize*(l+1)
           imin=1
           imax=points-1 !because f(x+2h) is necessary
         do i=imin,imax,idir
           pot=(Etrial-pott(i))/hb2m !insert V if V=0
           a1= 2.0E0_dp*(1.0E0_dp-5.0E0_dp/12.0E0_dp*pot*meshsize**2)
           a2=(1.0E0_dp+1.0E0_dp/12.0E0_dp*pot*meshsize**2)
           a3=(1.0E0_dp+1.0E0_dp/12.0E0_dp*pot*meshsize**2)
           triall(i+idir)= (a1*triall(i)-a2*triall(i-idir))/a3
         end do
!
! half right wavefunction
!         allocate(trialr(0:points), stat = ifail)
!         if (ifail .ne. 0) STOP
!
           idir=-1
           trialr(points)=0.0E0_dp
           trialr(points-1)= 1.0E-15_dp
           imin=1 !because f(x-2h) is necessary
           imax=(points-1)
         do i=imax,imin,idir
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
!         normall=0.0E0_dp
!         do i=1,(points+1)/2
!            normall = normall + triall(i)**2*meshsize
!         end do
!         print*, "normal L=", normall
!         triall(:)=triall(:)/sqrt(normall*2)
!         normalr=0.0E0_dp
!         do i=points-1,(points+1)/2,-1
!            normalr = normalr + trialr(i)**2*meshsize
!         end do
!         print*, "normal R=",normalr
!         trialr(:)=trialr(:)/sqrt(normalr*2)
!
!
            !if(derivr .lt. 1.0E-10_dp) then
              print*, "left function =", contil
              print*, "right function =", contir
              print*, "left derivative =", derivl
              print*, "right derivative =", derivr
          !  end if
            coeff= derivl/derivr
              print*, "coefficient =", coeff
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
           if(abs((contir-contil)/contir) .lt. 1.0E-6_dp .and.  &
             abs((derivr-derivl)/derivr) .lt. 1.0E-6_dp) then 
             print '("Merging successfull")'
             print *, "continuity: ", contil, contir
             print *, "derivative: ", derivl, derivr
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
!         allocate(trialwf(0:points), stat = ifail)
!         if (ifail .ne. 0) STOP           
!
           do i=0,i_sur
              trialwf(i)= triall(i)
           end do
           do i=i_sur+1, points
              trialwf(i)= trialr(i)
           end do
         end subroutine parwf
!
!
!
!potential function
         function potV(i,meshsize,Rmax,Rmin,a,kpot, Vvalue, &
                  numN, numZ,rzero,rProt,charge, l, jj) result(vr)
         use variables
         implicit none
 !        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
         real (kind = dp) :: meshsize, Rmax, Rmin,a,Vvalue, rzero, bigR, ferre
         real (kind=dp) :: vr, pos, leftlim, rightlim, pos0, rProt
         real (kind=dp) :: Uq_pot, Wq_pot1,Wq_pot2, V_Cou, V_cfug
         integer :: i, kpot, numZ, numN, points, l, jj, charge
         ! symmetric square well
         pos=Rmin+i*meshsize
         points=(Rmax-Rmin)/meshsize !redundant
         pos0=Rmin !redundant
         bigR= rzero*(numN+numZ)**(1.0E0_dp/3.0E0_dp)
!
         if(kpot .eq. 0) then
             vr=0         
         else if(kpot .eq. 1) then 
           leftlim=Rmin
           rightlim=a
           if(pos .ge. leftlim .and. pos .le. rightlim) then
             vr= Vvalue
           else
             vr= 0
           end if
           else if(kpot .eq. 2) then
             ferre=1/(1+exp((abs(pos-pos0)-bigR)/a))
             vr=-44.0192307692*ferre
 !            vr=-(51-33*(numN-numZ)/(numN+numZ))*ferre
         else if(kpot .eq. 3) then
             if(abs(pos-pos0) .le. rProt) then
             vr=numZ*esquare/2/rProt*(3-(pos/rProt)**2)         
             else if(abs(pos-pos0) .gt. rProt) then
             vr=numZ*esquare/rProt 
             end if 
         else if(kpot .eq. 4) then
             ferre=-0.44*exp((abs(pos-pos0)-bigR)/a) * rzero**2 &
                   /a/pos/(1+exp((abs(pos-pos0)-bigR)/a))**2
             vr=-(51-33*(numN-numZ)/(numN+numZ))*ferre   
         else if(kpot .eq. 5) then
             ! Uq potential is the Wood-Saxon
             if(i .eq. 0.1E0_dp)then
                pos = meshsize
             end if
             Uq_pot=- 44.0192307692 &
                    /(1+exp((abs(pos-pos0)-bigR)/a))
!                   -(51-33*(numN-numZ)/(numN+numZ)) &
!                   /(1+exp((abs(pos-pos0)-bigR)/a))
             !jj=2j so jj is integer
!             ferre=-0.44*exp((abs(pos-pos0)-bigR)/a) * rzero**2 &
!                   /a/(pos)/(1+exp((abs(pos-pos0)-bigR)/a))**2
!
             Wq_pot1=0.0E-0_dp!ferre*(dble(jj)/2*(dble(jj)/2+1)-l*(l+1)-3/4)/(pos) 
! 
             if(i .eq. 0) Wq_pot1=0.0E0_dp     
             Wq_pot2=-44.0192307692
                     !-(51-33*(numN-numZ)/(numN+numZ))
!
             V_cfug=hb2m*((l*(l+1))/(pos)**2)
             if(i .eq. 0) V_cfug=0.0E0_dp  
             vr=Uq_pot+Wq_pot1*Wq_pot2+V_cfug
! Coulomb part
             if(charge .eq. 1) then
               if(abs(pos-pos0) .le. rProt) then
                V_Cou=numZ*esquare/2/rProt*(3-(pos/rProt)**2)         
               else if(abs(pos-pos0) .gt. rProt) then
                V_Cou=numZ*esquare/rProt 
               end if
             vr= vr+V_Cou
             end if
         end if
         
end function potV
