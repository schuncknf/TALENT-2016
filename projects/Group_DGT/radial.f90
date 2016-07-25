      module variables
      implicit none
      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
      REAL(KIND=dp), PARAMETER :: pi = 3.1415926535897932384626433832795E0_dp
      REAL(KIND=dp), PARAMETER :: mc2=938.9059000E0_dp
      REAL(KIND=dp), PARAMETER :: hbarc = 197.32891000E0_dp
      REAL(KIND=dp), PARAMETER :: hb2m = 20.75E0_dp
      end module variables
!
!
!
      program numerov
      use variables
      implicit none
!     program to solve 1D Schroedinger equation
!      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
      REAL(kind=dp) :: Eup, Edown, epsil
      REAL(kind=dp) :: Rmin, Rmax, meshsize, a, Vvalue, rzero
      integer :: Nodemax, kpot, Node
      REAL(kind=dp) :: Etrial, Eexp
      integer :: points, ifail, converg, i
      integer :: Nodecount, numN, numZ
      REAL(kind=dp) :: a1, a2, a3, pot, normal, deltae
      REAL(kind=dp) :: potV
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: trialwf
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: pott
      CHARACTER(LEN=2):: Nodecc
!==========     set initial values     ==============================
      Eup=1000.E0_dp
      Edown=-100.E0_dp
      epsil=1.E-12_dp
      Nodemax=20
      Rmin=-20.E0_dp
      Rmax=20.E0_dp !fermi
      meshsize=0.001E0_dp
! kpot= 0 no potential      
! kpot= 1 square well of size a and deepth Vvalue
! kpot=2  Wood-Saxons potential
! kpot= 3 Coulomb
      numN= 8
      numZ= 8 
      kpot= 2
      a=0.67E0_dp
      Vvalue=-100.E0_dp 
      rzero=1.27E0_dp   
!==================================================================     
!
      points=(Rmax-Rmin)/meshsize
         allocate(trialwf(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
!         trialwf(:)=0.0E0_dp
         allocate(pott(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
         do i=0,points
         pott(i)=potV(i,meshsize,Rmax,Rmin,a,kpot, Vvalue, &
                  numN, numZ,rzero)
!         print *, pott(i)
         end do
!
         OPEN(UNIT=8, FILE='summary_table.dat', status='unknown')
         do Node=0, Nodemax
           call wavef(Rmax,Rmin,points,meshsize, Node, Eup, Edown, epsil, pott, &
                 trialwf, Etrial)

           call parwf(Rmax,Rmin,points,meshsize, Node, Etrial, &
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
         if ( Nodecount < 10 ) then
           write( Nodecc, '("0",i1)' ) Nodecount
         else
           write( Nodecc, '(i2)' ) Nodecount
         end if

         OPEN(UNIT=7, FILE='wavefunction'//Nodecc//'.dat', status='unknown')
         
         write(7,*) "Number of nodes=", Nodecount
         write(7,*) "Energy= ", Etrial
         do i=0,points
           write(7,*) Rmin+i*meshsize, trialwf(i)
         end do
         CLOSE(unit=7)
         print*, "------------------------------------- "
         write(8,'("Nod ",i3,10x,"Energy ", 1f18.10)') Nodecount, Etrial
         end do
         CLOSE(unit=8)
      end program numerov
! end
!
!
!
 

!wave function calculation
         subroutine wavef(Rmax,Rmin,points,meshsize, Node, Eup0, Edown0, epsil, pott, &
                 trialwf, Etrial) 
         use variables
         implicit none
         REAL(kind=dp), INTENT(IN) :: Eup0, Edown0, epsil
         REAL(kind=dp) :: Eup, Edown
         REAL(kind=dp), INTENT(IN)  :: Rmin, Rmax, meshsize
         integer, INTENT(IN)  :: Node
         REAL(kind=dp), INTENT(OUT) :: Etrial
         integer, INTENT(IN):: points
         integer:: converg
         integer :: Nodecount , ifail, i, idir, imin, imax, kpot
         REAL(kind=dp) :: a1, a2, a3, normal, Eexp
         REAL(kind=dp) , ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: trialwf
         REAL(kind=dp) :: pot, a, Vvalue
         REAL(kind=dp) , DIMENSION(0:points), INTENT(IN) :: pott
!
         allocate(trialwf(0:points), stat = ifail)
         if (ifail .ne. 0) STOP

         !loading the potential
           Eup=Eup0
           Edown=Edown0
         !
!         if(dir .eq. left) then
           idir= 1
           trialwf(0)=0.0E0_dp
           trialwf(1)= 0.1E0_dp
           imin=1
           imax=points-1
!         else if(dir .eq. hale) then
!           idir= 1
!           trialwf(0)=0.0E0_dp
!           trialwf(1)= 0.1E0_dp
!           imin=1
!           imax=(points+1)/2+2 !because f(x+2h) is necessary
!         else if (dir .eq. hari) then
!           idir=-1
!           trialwf(points)=0.0E0_dp
!           trialwf(points-1)= 0.1E0_dp
!           imin=(points+1)/2-2 !because f(x-2h) is necessary
!           imax=(points-1)
!         end if
!
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
         if(Nodecount .gt. Node) then
           Eup=Etrial
         else if(Nodecount .le. Node) then
           Edown=Etrial
         end if
! checking convergence
         if (abs(Eup-Edown) .lt. epsil) then
           print '("Convergence obtained for wavefunction")'
           print *, "with energy= ", Etrial
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
         subroutine parwf(Rmax,Rmin,points,meshsize, Node, Etrial, &
                    epsil, pott, trialwf) 
         use variables
         implicit none
         REAL(kind=dp), INTENT(IN) :: Etrial, epsil
         REAL(kind=dp), INTENT(IN)  :: Rmin, Rmax, meshsize
         integer, INTENT(IN)  :: Node
!         REAL(kind=dp):: Eup, Edown
         integer, INTENT(IN):: points
         integer :: Nodecount , ifail, i, idir, imin, imax, converg
         REAL(kind=dp) :: a1, a2, a3, normal, normalr, normall
         REAL(kind=dp), DIMENSION(0:points), INTENT(OUT) :: trialwf
         REAL(kind=dp), DIMENSION(0:points) ::  triall, trialr
         REAL(kind=dp) :: pot
         REAL(kind=dp), DIMENSION(0:points), INTENT(IN) :: pott
         REAL(kind=dp) :: derivr, derivl, contir, contil, coeff
!
! half left wavefunction
!         allocate(triall(0:points), stat = ifail)
!         if (ifail .ne. 0) STOP 
!           
           idir=1
           triall(0)=0.0E0_dp
           triall(1)= 0.1E0_dp
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
           trialr(points-1)= 0.1E0_dp
           imin=1 !because f(x-2h) is necessary
           imax=(points-1)
         do i=imax,imin,idir
           pot=(Etrial-pott(i))/hb2m !insert V if V=0
           a1= 2.0E0_dp*(1.0E0_dp-5.0E0_dp/12.0E0_dp*pot*meshsize**2)
           a2=(1.0E0_dp+1.0E0_dp/12.0E0_dp*pot*meshsize**2)
           a3=(1.0E0_dp+1.0E0_dp/12.0E0_dp*pot*meshsize**2)
           trialr(i+idir)= (a1*trialr(i)-a2*trialr(i-idir))/a3
         end do
! calculating derivatives and wave function values at central point
           i=(points+1)/2
           derivl= (triall(i+2)/12-2*triall(i+1)+2*triall(i-1)/3 &
                 -triall(i-2)/12)/4/meshsize
           contil=triall(i)
           derivr= (trialr(i+2)/12-2*trialr(i+1)+2*trialr(i-1)/3 &
                 -trialr(i-2)/12)/4/meshsize
           contir=trialr(i)
! 
! matching condition
!
         converg=0
         do while(converg .ne. 1)
!
         normall=0.0E0_dp
         do i=1,(points+1)/2
            normall = normall + triall(i)**2*meshsize
         end do
!         print*, "normal L=", normall
         triall(:)=triall(:)/sqrt(normall*2)
         normalr=0.0E0_dp
         do i=points-1,(points+1)/2,-1
            normalr = normalr + trialr(i)**2*meshsize
         end do
!         print*, "normal R=",normalr
         trialr(:)=trialr(:)/sqrt(normalr*2)
!
!
           if((-1)**Node .eq. -1) then
            if(derivr .lt. 1.0E-10_dp) then
              print*, "dividing by zero"
            end if
            coeff= derivl/derivr
            trialr(:)= coeff*trialr(:)
           else if((-1)**Node .eq. 1) then
            if(contir .lt. 1.0E-10_dp) then
              print*, "dividing by zero"
            end if
            coeff= contil/contir
            trialr(:)= coeff*trialr(:)
           end if
!

!
!check values
           i=(points+1)/2
           derivl= (triall(i+2)/12-2*triall(i+1)+2*triall(i-1)/3 &
                 -triall(i-2)/12)/4/meshsize
           contil=triall(i)
           derivr= (trialr(i+2)/12-2*trialr(i+1)+2*trialr(i-1)/3 &
                 -trialr(i-2)/12)/4/meshsize
           contir=trialr(i)
!           
           if(abs((contir-contil)) .lt. 1.0E-6_dp .and.  &
             abs((derivr-derivl)) .lt. 1.0E-6_dp) then 
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
           do i=0,(points+1)/2
              trialwf(i)= triall(i)
           end do
           do i=((points+1)/2)+1, points
              trialwf(i)= trialr(i)
           end do
         end subroutine parwf
!
!
!
!potential function
         function potV(i,meshsize,Rmax,Rmin,a,kpot, Vvalue, &
                  numN, numZ,rzero) result(vr)
         implicit none
         INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
         real (kind = dp) :: meshsize, Rmax, Rmin,a,Vvalue, rzero, bigR, ferre
         real (kind=dp) :: vr, pos, leftlim, rightlim, pos0
         integer :: i, kpot, numZ, numN, points
         ! symmetric square well
         pos=Rmin+i*meshsize
         points=(Rmax-Rmin)/meshsize !redundant
         pos0=Rmin+(points+1)/2*meshsize !redundant
         bigR= rzero *(numN+numZ)**(1.0E0_dp/3.0E0_dp)
!
         if(kpot .eq. 0) then
             vr=0         
         else if(kpot .eq. 1) then 
           leftlim= (Rmax+Rmin)/2-a/2
           rightlim=(Rmax+Rmin)/2+a/2
           if(pos .ge. leftlim .and. pos .le. rightlim) then
             vr= Vvalue
           else
             vr= 0
           end if
         else if(kpot .eq. 2) then
             ferre=1/(1+exp((abs(pos-pos0)-bigR)/a))
             vr=-(51-33*(numN-numZ)/(numN+numZ))*ferre
         end if
end function potV
