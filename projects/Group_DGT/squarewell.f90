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
      REAL(kind=dp) :: Rmin, Rmax, meshsize, a, Vvalue
      integer :: Node, kpot
      REAL(kind=dp) :: Etrial, Eexp
      integer :: points, ifail, converg, i
      integer :: Nodecount
      REAL(kind=dp) :: a1, a2, a3, pot, normal, deltae
      REAL(kind=dp) :: potV
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: trialwf
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: pott
!     CHARACTER(LEN=5), DIMENSION(7) :: char
!==========     set initial values     ==============================
      Eup=100.E0_dp
      Edown=-100.E0_dp
      epsil=1.E-12_dp
      Node=0
      Rmin=-20.E0_dp
      Rmax=20.E0_dp !fermi
      meshsize=0.001E0_dp
! kpot= 0 no potential      
! kpot= 1 square well of size a and deepth Vvalue
      kpot= 1
      a=6.E0_dp
      Vvalue=-100.E0_dp    
!==================================================================     
!
      points=(Rmax-Rmin)/meshsize
! allocate the trialwf
      allocate(trialwf(0:points), stat = ifail)
      if (ifail .ne. 0) STOP
!loading the potential
      allocate(pott(0:points), stat = ifail)
      if (ifail .ne. 0) STOP
      do i=0,points
         pott(i)=potV(i,meshsize,Rmax,Rmin,a,kpot, Vvalue)
!         print *, pott(i)
      end do
!loop until convergence
      converg=0
      do while (converg .ne. 1)
         Etrial=(Eup+Edown)/2
!         print*, Etrial
!building up the wave function up to Rmax
         trialwf(0)=0.0E0_dp
         trialwf(1)= 0.1E0_dp
         do i=1,points-1
           pot=(Etrial-pott(i))/hb2m !insert V if V=0
           a1= 2.0E0_dp*(1.0E0_dp-5.0E0_dp/12.0E0_dp*pot*meshsize**2)
           a2=(1.0E0_dp+1.0E0_dp/12.0E0_dp*pot*meshsize**2)
           a3=(1.0E0_dp+1.0E0_dp/12.0E0_dp*pot*meshsize**2)
           trialwf(i+1)= (a1*trialwf(i)-a2*trialwf(i-1))/a3
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
           print '("Convergence obtained")'
           print *, "Obtained energy= ", Etrial
! expected energy from the formula Eq. 1.14 in manuale_hf.pdf
          if (kpot .eq. 0) then 
            Eexp= ((Node+1)*pi)**2*hb2m/(Rmax-Rmin)**2
            print *, "Expected energy= ", Eexp
          end if
           converg=1
         end if
      end do
!renormalization
         normal=0.0E0_dp
         do i=0,points
            normal = normal + trialwf(i)**2*meshsize
         end do
         trialwf(:)=trialwf(:)/sqrt(normal)
!
         print '("Wavefunction normalized")'
! output number of nodes excluded a possible node in the last point of the box
! normalized wavefunction
         OPEN(UNIT=7, FILE='wavefunction.dat', status='unknown')
         if(trialwf(points-1)*trialwf(points) .lt. 0.0E0_dp) &
             Nodecount= Nodecount-1 
         print*, Nodecount
         write(7,*) "Number of nodes=", Nodecount
         do i=0,points
           write(7,*) Rmin+i*meshsize, trialwf(i)
         end do
         CLOSE(unit=7)
      end program numerov
! end
!
!
!
!potential function
         function potV(i,meshsize,Rmax,Rmin,a,kpot, Vvalue) result(vr)
         implicit none
         INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
         real (kind = dp) :: meshsize, Rmax, Rmin,a,Vvalue
         real (kind=dp) :: vr, pos, leftlim, rightlim
         integer :: i, kpot
         ! symmetric square well
         pos=Rmin+i*meshsize
         if(kpot .eq. 1) then 
           leftlim= (Rmax+Rmin)/2-a/2
           rightlim=(Rmax+Rmin)/2+a/2
           if(pos .ge. leftlim .and. pos .le. rightlim) then
             vr= Vvalue
           else
             vr= 0
           end if
         else if(kpot .eq. 0) then
             vr=0
         end if
         end function potV
