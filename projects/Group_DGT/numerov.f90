      program numerov
      implicit none
!     program to solve 1D Schroedinger equation
      REAL(kind=8) :: Eup, Edown, epsil, Etrial
      integer :: Node, try, points, ifail, converg, i
      integer :: Nodecount
      REAL(kind=8) :: a1, a2, a3, pot, normal, deltae
      REAL(kind=8) :: Rmin, Rmax, meshsize
      REAL(kind=8) , ALLOCATABLE, DIMENSION(:) :: trialwf
      REAL(KIND=8), PARAMETER :: pi = 3.1415926535897932384626433832795d0
      REAL(KIND=8), PARAMETER :: mc2=938.9059000d0
      REAL(KIND=8), PARAMETER :: hbarc = 197.32891000d0
      REAL(KIND=8), PARAMETER :: hb2m = 20.75d0
!     CHARACTER(LEN=5), DIMENSION(7) :: pname, qname
!     set initial values
      Eup=100.0d0
      Edown=0.0d0
      epsil=0.0000000000001d0
      Node=0
      Rmin=0
      Rmax=6 !fermi
      meshsize=0.1d0
!
      points=(Rmax-Rmin)/meshsize
! allocate the trialwf
      allocate(trialwf(0:points), stat = ifail)
      if (ifail .ne. 0) STOP
!loop until convergence
      converg=0
      do while (converg .ne. 1)
         Etrial=(Eup+Edown)/2
         print*, Etrial
!building up the wave function
         trialwf(0)=0.0d0
         trialwf(1)= 0.1d0
         do i=1,points-1
           pot=Etrial/hb2m !insert V if V=0
           a1= 2.0d0*(1.0d0-5.0d0/12.0d0*pot*meshsize**2)
           a2=(1.0d0+1.0d0/12.0d0*pot*meshsize**2)
           a3=(1.0d0+1.0d0/12.0d0*pot*meshsize**2)
           trialwf(i+1)= (a1*trialwf(i)-a2*trialwf(i-1))/a3
         end do
! counting nodes
         Nodecount=0    
         do i=1,points
           if(trialwf(i-1)*trialwf(i) .lt. 0.0d0) &
             Nodecount= Nodecount+1
         end do 
         print*, Nodecount
!checking nodes    
         if(Nodecount .gt. Node) then
           Eup=Etrial
         else if(Nodecount .le. Node) then
           Edown=Etrial
         end if
! checking convergence
         if (abs(Eup-Edown) .lt. epsil) then
           print '("Convergence obtained")'
           print *, Eup, Edown, Etrial
           converg=1
         end if
      end do
!renormalization
         normal=0.0d0
         do i=0,points
            normal = normal + trialwf(i)**2*meshsize
         end do
         trialwf(:)=trialwf(:)/sqrt(normal)
!
         print '("Wavefunction normalized")'
         OPEN(UNIT=7, FILE='wavefunction.dat', status='unknown')
         if(trialwf(points-1)*trialwf(points) .lt. 0.0d0) &
             Nodecount= Nodecount-1 
         print*, Nodecount
         write(7,*) "Number of nodes=", Nodecount
         do i=0,points
           write(7,*) Rmin+i*meshsize, trialwf(i)
         end do
         CLOSE(unit=7)
      end program numerov
