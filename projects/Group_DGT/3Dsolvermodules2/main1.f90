!===========================================================================
!This program solves Schrödinger equation by using numerov method and shooting.
!Schrödinger equation can be solved for infinite square well, finite square well, wood-saxon potential etc.
!
! run e.g.:
! gfortran -fno-whole-file main.f90 variables.f90 potentials.f90 calcwaves.f90 -o main
!
!------------------------------------------------------------
!MAIN PROGRAM:
!------------------------------------------------------------
	! Used variables:
	! Eup = Upper limit for energy in numerov algorithm
	! Edown= Lower limit for energy in numerov algorithm
	! epsil= Demanded precision for convergence
	! Rmax= Right hand side limit for the box in which we are solving SE
	! meshsize= Distance between meshpoints
	! a= diffusivity (0.67 fm)
	! Vvalue= the depth of the square well
	! rzero= a constant (1.27 fm)
	! Nodemax=the maximum number of nodes 
	! N_max = Max number of principal quantum number
	! N_maxN= Max number of principal q.n. for neutrons
	! N_maxZ= Max number of principal q.n. for protons
	! N_act= Active principal quantum number (for loops)
	! Num_par = number of particles (for loops)
	! Num_par_max= maximum number of particles
	! kpot= index for different potentials, see below
	! Etrial= Energy of the trial wavefunction
	! Eexp = expected energy eigenvalue (analytical solution, infinite square well)
	! points= maximum number of points for e.g. wavefunction
	! ifail= logical argument for checking if something is failed
	! converg= logical argument for checking convergence (0=not converged,1=converged)
	! i=index for loops
	! charge= separator for neutrons and protons (0=neutron, 1=proton)
	! l=orbital quantum number
	! jj= "2*j", so two times total angular momenta
	! iii= index for counting solutions
	! jmin= lower bound for the angular momenta
	! jmax= upper bound for the angular momenta
	! Nodecount = number of nodes
	! numN= number of neutrons
	! numZ= number of protons
	! strval= ??
	! n_rad = radial quantum number
	! st_post = ??
	! minN = index for searching minimum energy of neutron states
	! minZ= index for searching minimum energy of proton states
	! nstatN= number of neutron states
	! nstatZ= number of proton states
	! num_stat= total number of states
	! loop = ??
	! HF_conv = ??
	! a1= coefficient in numerov algorithm (manuale p.7) 
	! a2= coefficient in numerov algorithm (manuale p.7)
	! a3= coefficient in numeror algorithm (manuale p.7)
	! pot= potential
	! normal= variable for normalization (integral over wavefunction)
	! potV=value of potential at some point
	! rProt=radius (of protons)
	! trialwf= trial wafefunction
	! pott= the whole potential (array)
	! Nodecc= a variable for indexing all different solutions (different N,l..)
     program numerov
      use variables
      use potentials
      use wavefunctions

      implicit none

      REAL(kind=dp) :: Eup, Edown, epsil
      REAL(kind=dp) :: Rmax, meshsize, a, Vvalue, rzero
      integer :: Nodemax, kpot, N_max, N_maxN, N_maxZ, N_act, Num_par
      integer :: Num_par_max
      REAL(kind=dp) :: Etrial, Eexp
      integer :: points, ifail, converg, i,j, charge, l, jj, nstat, jmin, jmax
      integer :: Nodecount, numN, numZ, strval, n_rad, st_pos, minN
      integer :: minZ, nstatN, nstatZ, num_stat, loop, HF_conv
      REAL(kind=dp) :: a1, a2, a3, pot, normal, deltae, pos
      REAL(kind=dp) :: tot_E, t0, t3, alpha, E_before
      REAL(kind=dp) :: rProt, valueR0, den_int, denN_int, denZ_int
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: trialwf, posiR
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: pott, density, denN, denZ
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: pre_dens, pre_denN, pre_denZ
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:,:) :: states, states2
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: ordEneN, ordEneZ, occup
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: der_u, der_uN, der_uZ
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: kinden, kindenN, kindenZ
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: E_cou, E_cou1
      REAL(kind=dp) ::  Crho0,  Crho1, E_cou2
      CHARACTER(LEN=3):: char_nstat
      CHARACTER(LEN=5) :: string
      CHARACTER(LEN=3):: part
!==========     set initial values     ==============================
      Eup=100.E0_dp
      Edown=-100.E0_dp
      epsil=1.E-6_dp
      Rmax=30.E0_dp !fermi
      meshsize=0.08E0_dp
! kpot= 0 no potential      
! kpot= 1 square well of size a and deepth Vvalue
! kpot= 2 Wood-Saxon potential
! kpot= 3 Coulomb
! kpot= 4 spin-orbit W potential
! kpot= 5 Radial potential Eq. 2.10
!
      numN= 8
      numZ= 8
      kpot= 5
      a=0.67E0_dp
      Vvalue=-100.E0_dp 
      rzero=1.27E0_dp 

      rProt= rzero*(numZ**(1.0E0_dp/3.0E0_dp))
      N_maxN=2
      N_maxZ=2
!==================================================================     
!
! Calculates the number of points in a mesh and allocates the trial wavefunction
! the potential, position array, densities
! Gives error message if there are problems

      points=int(Rmax/meshsize)

         allocate(trialwf(0:points), stat = ifail)
         if (ifail .ne. 0) STOP

!         trialwf(:)=0.0E0_dp
         allocate(pott(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
!
         allocate(posiR(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
!
         allocate(density(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
!
         allocate(denN(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
!
         allocate(denZ(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
!
         allocate(der_u(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
!
         allocate(der_uN(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
!
         allocate(der_uZ(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
!
         allocate(kinden(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
!
         allocate(kindenN(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
!
         allocate(kindenZ(0:points), stat = ifail)
         if (ifail .ne. 0) STOP

         allocate(pre_dens(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
!
         allocate(pre_denN(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
!
         allocate(pre_denZ(0:points), stat = ifail)
         if (ifail .ne. 0) STOP

! Define array of R
         do i=0,points
         posiR(i)=i*meshsize
         end do
	 posiR(0) = 0.01    ! We don't want R=0 at all!

         num_stat=0
         nstatN=0
         nstatZ=0
         do charge=0,1
         if(charge .eq. 0) then
           N_max=N_maxN
         else if(charge .eq. 1) then
           N_max=N_maxZ
         end if


!Number of states
         do N_act=0,N_max   
          do n_rad=0,N_act
            do l= 0, N_act
               if((2*n_rad+l) .ne. N_act) then
                go to 21
               end if 
               jmin=abs(2*l-1)
               do jj= jmin,2*l+1,2
                  num_stat=num_stat+1
                if(charge .eq. 0) then
                 nstatN=nstatN+1
                else if(charge .eq. 1) then
                 nstatZ=nstatZ+1
                end if
               end do!j
21            continue
            end do!l
          end do!n_rad
          end do!N_act
          end do !charge



         allocate(states(1:points+6,1:num_stat), stat = ifail)
         if (ifail .ne. 0) STOP
         allocate(states2(1:points+6,1:num_stat), stat = ifail)
         if (ifail .ne. 0) STOP
         allocate(ordEneN(1:nstatN), stat = ifail)
         if (ifail .ne. 0) STOP
         allocate(ordEneZ(1:nstatZ), stat = ifail)
         if (ifail .ne. 0) STOP

         allocate(occup(1:num_stat), stat = ifail)
         if (ifail .ne. 0) STOP

!states contain Energy,charge, nrad, l, jj, trialwf(0:points) then index of state
!
         loop=0
         HF_conv=0
         E_before=0 !CHANGED
        !do loop =0,10
        do while (HF_conv .ne. 1)
      nstat=0
	
!

! This opens a summary table and begins the actual program, 
! solves the problem for different combination of l and n
      OPEN(UNIT=8, FILE='summary_table.dat', status='unknown')
      do charge=0,1
         if(charge .eq. 0) then
           N_max=N_maxN
           Num_par_max=numN
           Num_par=0
         else if(charge .eq. 1) then
           N_max=N_maxZ
           Num_par_max=numZ
           Num_par=0
         end if

!
!Starts to calculate the number of states which needed (for neutrons, protons, together)

         do N_act=0,N_max   
         do n_rad=0,N_act
            do l= 0, N_act
               if((2*n_rad+l) .ne. N_act) then
                go to 31
               end if 
               jmin=abs(2*l-1)
               do jj= jmin,2*l+1,2
               nstat=nstat+1

               
! Finds out the value of potential at a point                 
         OPEN(UNIT=9, FILE='potential.dat', status='unknown')

         if(loop .eq. 0) then
           kpot=5
         else if (loop .ne. 0) then
           kpot=6
         end if


         do i=0,points
         pott(i)=potV(i,posiR,points,meshsize, Rmax, a,kpot, Vvalue, &
                  numN, numZ,rzero,rProt,charge, l, jj, &
                  pre_dens, pre_denZ, pre_denN)
         write(9,*) posiR(i), pott(i)
         end do
         CLOSE(unit=9)
!

           call wavef(points,meshsize, n_rad, l, Eup, Edown, epsil, pott, &
                 trialwf, Etrial) ! WEDO: U_r

           call parwf(points,meshsize, n_rad, l, numN, numZ, Etrial, &
                    epsil, pott, trialwf) ! WEDO: U_r


! Calculates the integral over trial wavefunction and normalizes the wafefunction

         normal=0.0E0_dp
!         if(l .eq. 0 .and. kpot .eq. 5) then
!            trialwf(0)=trialwf(1)
!         end if
         do i=0,points-1
            normal = normal + trialwf(i)**2*meshsize !trialwf(r)= u_r(r)
         end do
         do i=0,points
         trialwf(i)=trialwf(i)/sqrt(normal)
         end do
!         if(l .eq. 0 .and. kpot .eq. 5) then
!            trialwf(0)=trialwf(1)
!         end if
!
!         print '("Wavefunction normalized")'
!         print*, normal


! Calculates the output number of nodes excluded a possible node in the last point of the box
         Nodecount=0    
         do i=1,points !WEDO:-1
           if(trialwf(i-1)*trialwf(i) .lt. 0.0E0_dp) &
             Nodecount= Nodecount+1
         end do 
         if(trialwf(points-1)*trialwf(points) .lt. 0.0E0_dp) &
             Nodecount= Nodecount-1 
!
! allocating super matrix
         states(1,nstat)=Etrial
         states(2,nstat)=charge
         states(3,nstat)=n_rad
         states(4,nstat)=l
         states(5,nstat)=jj
         do i= 0,points
            states(6+i,nstat)=trialwf(i) ! WEDO: U_r
         end do
                end do !loop j
31            continue
              end do ! loop l
            end do !  loop n_rad
         end do ! N_max
         end do ! charge             


! This defines an array which includes all the energies (for calculating the correct order)
         ordEneN(:)=states(1,1:nstatN)
         ordEneZ(:)=states(1,nstatN+1:num_stat)



! Finds out the correct order of states by finding the lowest value of energy in
! the array and then "deleting" that by inserting maximum energyvalue Eup in the 
! array of energies
! Correct order of energies is stored in the matrix "states2"
         do i=1,nstatN
            minN= MINLOC(ordEneN,dim=1)
            states2(:,i) = states(:,minN)
            ordEneN(minN)=Eup
         end do
         do i=nstatN+1, num_stat
            minZ= MINLOC(ordEneZ,dim=1)
            states2(:,i) = states(:,nstatN+minZ)
            ordEneZ(minZ)=Eup
         end do


! Modifies the variable "char_nstat" to be in a correct form, eg. "001", "013" 
! or "800"
          do i=1,num_stat
!            print*, states2(1,i)
         if ( i < 10 ) then
           write( char_nstat, '("00",i1)' ) i
         else if ( i .ge. 10 .and. i .lt. 100 ) then
           write( char_nstat, '("0",i2)' ) i
         else if ( i .ge. 100 .and. i .lt. 1000 ) then
           write( char_nstat, '(i3)' ) i
         else if ( i .gt. 10000 ) then
           write( char_nstat, '("###")' )
         end if

! Routine for writing meshpoint and square of "wavefunction"
!
         OPEN(UNIT=7, FILE='states'//char_nstat//'.dat', status='unknown')
!         
         write(7,*) "#Radial number=", states2(2,i)
         write(7,*) "#Energy= ", states2(1,i)
         do j=6,points
           write(7,*) posiR(j-6), states2(j,i), states2(j,i)/(posiR(j-6))
         end do
         CLOSE(unit=7)
         if (int(states2(2,i)) .eq. 0) then
            part="n"
         else if (int(states2(2,i)) .eq. 1) then
            part="p"
         end if
         print*, "------------------------------------------------------------- "
       print'(a, "  n_rad= ",i3, " l=",i3, " 2j=",i3,4x,"sp_en= ",1f18.10)', &
          part,  int(states2(3,i)), int(states2(4,i)), int(states2(5,i)), states2(1,i)
      write(8,'(a," n_rd_",i3,2x,"l",i3,2x,"2j",i2,4x,"u(r)_",a,4x,"sp_en= ",1f18.10)') &
                part, int(states2(3,i)), int(states2(4,i)),  &
                int(states2(5,i)), char_nstat, states2(1,i)

          end do
         CLOSE(unit=8)
!
      density(:)=0.0E0_dp
      denN(:)=0.0E0_dp
      denZ(:)=0.0E0_dp
      occup(:)=0
 
        kinden(:)=0_dp
	kindenN(:)=0_dp
	kindenZ(:)=0_dp
        do j=0,points-2 !WEDO

! HERE SOME COMMENTS
!loop on neutron density
          charge=0
           N_max=N_maxN
           Num_par_max=numN
           Num_par=0
         do i=1, nstatN
         if(Num_par+int(states2(5,i)+1) .lt. Num_par_max) then
           occup(i)=(states2(5,i)+1)
         else if (Num_par .le. Num_par_max) then
           occup(i)=Num_par_max-Num_par
         else if (Num_par .gt. Num_par_max) then
           occup(i)=0
         end if
         Num_par= Num_par+int(states2(5,i)+1)
!         IT CAN BE CALCULATED FROM RIGHT
!
            density(j)= density(j) + (occup(i))*states2(6+j,i)**2/4/pi &
                        /posiR(j)**2
!
          ! der_u(j)= (states2(6+j+1,i)-states2(6+j,i))/meshsize
	   der_u(j) = (2*states2(6+j+1,i)-3*states2(6+j,i)/2.-states(6+j+2,i)/2.)/(2.*meshsize)
!
            kinden(j)= kinden(j)+occup(i)/posiR(j)**2/4/pi*    &
                     ((der_u(j)-states2(6+j,i)/posiR(j))**2+ &
                     states2(4,i)*(states2(4,i)+1)/posiR(j)**2*states2(6+j,i)**2)

            denN(j)= denN(j) + (occup(i))*states2(6+j,i)**2/4/pi &
                        /posiR(j)**2
!
           !der_uN(j)= (states2(6+j+1,i)-states2(6+j,i))/meshsize
           der_uN(j) = (2*states2(6+j+1,i)-3*states2(6+j,i)/2.-states(6+j+2,i)/2.)/(2.*meshsize)
!
            kindenN(j)= kindenN(j)+occup(i)/posiR(j)**2/4/pi*    &
                     ((der_u(j)-states2(6+j,i)/posiR(j))**2+ &
                     states2(4,i)*(states2(4,i)+1)/posiR(j)**2*states2(6+j,i)**2)

          !  if((j .eq. 1 .or. j .eq. 2) .and. int(states2(4,i)) .eq. 0) then
          !  density(0)= density(0) +(dble(occup))*states2(10,i)/4/pi &
          !              /posiR(j)**2
          !  denN(0)= denN(0) +(dble(occup))*states2(10,i)/4/pi &
          !              /posiR(j)**2 WEDO
          !  end if
         end do
         end do
         
!loop on proton density
         do j=0,points-2
           charge=1
           N_max=N_maxZ
           Num_par_max=numZ
           Num_par=0
         do i=nstatN+1, num_stat
         if(Num_par+(states2(5,i)+1) .lt. Num_par_max) then
           occup(i)=(states2(5,i)+1)
         else if (Num_par .le. Num_par_max) then
           occup(i)=Num_par_max-Num_par
         else if (Num_par .gt. Num_par_max) then
           occup(i)=0
         end if
         Num_par= Num_par+(states2(5,i)+1)
!         print*, int(states2(5,i)+1),occup
!
            density(j)= density(j) + (occup(i))*states2(6+j,i)**2/4/pi &
                        /posiR(j)**2 
!
            !der_u(j)= (states2(6+j+1,i)-states2(6+j,i))/meshsize
	    der_u(j) = (2*states2(6+j+1,i)-3*states2(6+j,i)/2.-states(6+j+2,i)/2.)/(2.*meshsize)

!
            kinden(j)= kinden(j)+(occup(i))/posiR(j)**2/4/pi*    &
                     ((der_u(j)-states2(6+j,i)/posiR(j))**2+ &
                     states2(4,i)*(states2(4,i)+1)/posiR(j)**2*states2(6+j,i)**2)

            denZ(j)= denZ(j) + (occup(i))*states2(6+j,i)**2/4/pi &
                        /posiR(j)**2
!
            !der_uZ(j)= (states2(6+j+1,i)-states2(6+j,i))/meshsize
	     der_uZ(j) = (2*states2(6+j+1,i)-3*states2(6+j,i)/2.-states(6+j+2,i)/2.)/(2.*meshsize)
!
            kindenZ(j)= kindenZ(j)+occup(i)/posiR(j)**2/4/pi*    &
                     ((der_u(j)-states2(6+j,i)/posiR(j))**2+ &
                     states2(4,i)*(states2(4,i)+1)/posiR(j)**2*states2(6+j,i)**2)

         end do
         end do

         den_int=0.
         denN_int=0.
         denZ_int=0.
         do i=0, points-1
         den_int= den_int+ density(i)*meshsize*4*pi*posiR(i)**2
         denN_int= denN_int+ denN(i)*meshsize*4*pi*posiR(i)**2
         denZ_int= denZ_int+ denZ(i)*meshsize*4*pi*posiR(i)**2
         end do
         print*, "total density=", den_int
         print*, "N density=", denN_int
         print*, "Z density=", denZ_int
!
! print density matrix
!
         OPEN(UNIT=10, FILE='density.dat', status='unknown')
         do i=0,points
           write(10,*) posiR(i), density(i), denN(i), denZ(i)
         end do
         CLOSE(unit=10)

         OPEN(UNIT=10, FILE='TAU.dat', status='unknown')
         do i=0,points
           write(10,*) posiR(i), kinden(i), kindenN(i), kindenZ(i)
         end do
         CLOSE(unit=10)


!        Coulomb density energy  WE DO
!         do i=0, points
!              E_cou2=E_cou2+ meshsize*posiR(i)*denZ(i)
!         end do
!
!         allocate(E_cou(0:points), stat = ifail)
!         if (ifail .ne. 0) STOP
!         allocate(E_cou1(0:points), stat = ifail)
!         if (ifail .ne. 0) STOP
!
!
!         E_cou(:)=0.0d0
!         E_cou1(:)=0.0d0
!         do j=1,points
!            do i=0, j
!              E_cou1(j)=E_cou1(j)+ meshsize*(posiR(i)**2*denZ(i)/posiR(j) &
!                      -posiR(i)*2*denZ(i))
!            end do
!            E_cou(j)=4.0d0*pi*esquare*(E_cou1(j)+E_cou2)
!            print*, E_cou(j)
!         end do! 

         OPEN(UNIT=17, FILE='denkin.dat', status='unknown')
         do i=0,points
           write(17,*) posiR(i), kinden(i)
         end do
         CLOSE(unit=17)
!   
         alpha=1.
         tot_E=0.0d0
         do j=1,points
            Crho0=3./8.*(-1132.4) + 3./48.*23610.4* &
              density(j)**alpha
            Crho1=-1./8.*(-1132.4)-1./48.*23610.4*&
              (denN(j)-denZ(j))**alpha
            tot_E= tot_E + 4*pi*meshsize*(posiR(j))**2 &
                *( hb2m* kinden(j) + &
             (Crho0)*density(j)**2+ (Crho1)*(denN(j)-denZ(j))**2)!+ &
 !             E_cou(j)*denZ(j)/2.0d0)
!         print*, tot_E
         end do
         print*, "###############################"        
         print*, tot_E
         print*, "###############################"      
! store density for next run
         pre_dens(:)=density(:)
!
         pre_denN(:)=denN(:)
!
         pre_denZ(:)=denZ(:)
          
           if (abs(tot_E-E_before) .lt. 0.0001) then
              HF_conv=1
              print*, "We are happy"
           end if
            E_before=tot_E
          loop= loop+1
         end do

         deallocate(trialwf)
         deallocate(pott)
      end program numerov
! end
!
!
!
!


