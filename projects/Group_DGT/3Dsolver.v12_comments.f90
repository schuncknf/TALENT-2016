!===========================================================================
!This program solves Schrödinger equation by using numerov method and shooting.
!Schrödinger equation can be solved for infinite square well, finite square well, wood-saxon potential etc.
!------------------------------------------------------------
!MAIN PROGRAM:
!------------------------------------------------------------
	! Used variables:
	! Eup = Upper limit for energy in numerov algorithm
	! Edown= Lower limit for energy in numerov algorithm
	! epsil= Demanded precision for convergence
	! Rmin= Left hand side limit for the box in which we are solving SE
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
      implicit none
      REAL(kind=dp) :: Eup, Edown, epsil
      REAL(kind=dp) :: Rmin, Rmax, meshsize, a, Vvalue, rzero
      integer :: Nodemax, kpot, N_max, N_maxN, N_maxZ, N_act, Num_par
      integer :: Num_par_max
      REAL(kind=dp) :: Etrial, Eexp
      integer :: points, ifail, converg, i,j, charge, l, jj, nstat, jmin, jmax
      integer :: Nodecount, numN, numZ, strval, n_rad, st_pos, minN
      integer :: minZ, nstatN, nstatZ, num_stat, loop, HF_conv
      REAL(kind=dp) :: a1, a2, a3, pot, normal, deltae, pos 
      REAL(kind=dp) :: tot_E, t0, t3, alpha
      REAL(kind=dp) :: potV, rProt, valueR0, den_int, denN_int, denZ_int
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: trialwf, posiR
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: pott, density, denN, denZ
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:,:) :: states, states2
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: ordEneN, ordEneZ, occup
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: der_u, der_uN, der_uZ
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: kinden, kindenN, kindenZ
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: E_cou, E_cou1
      REAL(kind=dp) ::  Crho0,  Crho1, E_cou2
      CHARACTER(LEN=3):: char_nstat
      CHARACTER(LEN=3):: part
!==========     set initial values     ==============================
      Eup=100.E0_dp
      Edown=-100.E0_dp
      epsil=1.E-6_dp
      Rmin=0.1E0_dp
      Rmax=30.E0_dp !fermi
      meshsize=0.1E0_dp
      kpot= 5
	! kpot= 0 no potential      
	! kpot= 1 square well of size a and deepth Vvalue
	! kpot= 2 Wood-Saxon potential
	! kpot= 3 Coulomb
	! kpot= 4 spin-orbit W potential
	! kpot= 5 Radial potential Eq. 2.10
      numN= 126
      numZ= 82 
      a=0.67E0_dp
      Vvalue=-100.E0_dp 
      rzero=1.27E0_dp 

      rProt= rzero*(numZ**(1.0E0_dp/3.0E0_dp))
      N_maxN=8
      N_maxZ=8
!==================================================================     
!
! Calculates the number of points in a mesh and allocates the trial wavefunction
! the potential, position array, densities
! Gives error message if there are problems
      points=int((Rmax-Rmin)/meshsize)

         allocate(trialwf(0:points), stat = ifail)
         if (ifail .ne. 0) STOP

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
!
         do i=0,points
         posiR(i)=Rmin+i*meshsize
         end do

         num_stat=0
         nstatN=0
         nstatZ=0
         do charge=0,1
         if(charge .eq. 0) then
           N_max=N_maxN
         else if(charge .eq. 1) then
           N_max=N_maxZ
         end if
!
!Calculates the number of states which needed (for neutrons, protons, together)
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
!
!
!
         loop=0
         HF_conv=0
         E_before=0.0d0
         do while (HF_conv .ne. 1)
           if(loop .eq. 0) then
             kpot=5
           else
             kpot=5 !HF fields
           end if
         allocate(states(1:points+7,1:num_stat), stat = ifail)
         if (ifail .ne. 0) STOP
         allocate(states2(1:points+7,1:num_stat), stat = ifail)
         if (ifail .ne. 0) STOP
         allocate(ordEneN(1:nstatN), stat = ifail)
         if (ifail .ne. 0) STOP
         allocate(ordEneZ(1:nstatZ), stat = ifail)
         if (ifail .ne. 0) STOP
!states contain Energy,charge, nrad, l, jj, f(0), trialwf(0:points) then index of state

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
!         OPEN(UNIT=9, FILE='potential.dat', status='unknown')
         do i=0,points
         pott(i)=potV(i,meshsize,Rmax,Rmin,a,kpot, Vvalue, &
                  numN, numZ,rzero,rProt,charge, l, jj)
!         write(9,*) Rmin+i*meshsize, pott(i)
         end do
!         CLOSE(unit=9)
!

           call wavef(points,meshsize, n_rad, l, Eup, Edown, epsil, pott, &
                 trialwf, Etrial) 

           call parwf(points,meshsize, n_rad, l, numN, numZ, Etrial, &
                    epsil, pott, trialwf) 

! Calculates the integral over trial wavefunction and normalizes the wafefunction
         normal=0.0E0_dp

         do i=0,points-1
            normal = normal + trialwf(i)**2*meshsize !trialwf(r)= u_r(r)
         end do
         do i=0,points
         trialwf(i)=trialwf(i)/sqrt(normal)
         end do


! Calculates the output number of nodes excluded a possible node in the last point of the box
         Nodecount=0    
         do i=1,points
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
         states(6,nstat)=0.0E0_dp
         do i= 0,points
            states(7+i,nstat)=trialwf(i)  !trialwf=u(r)
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
!         OPEN(UNIT=7, FILE='u(r)_state'//char_nstat//'.dat', status='unknown')
!         
!         write(7,*) "#Radial number=", states2(2,i)
!         write(7,*) "#Energy= ", states2(1,i)
!         do j=7,points
!           write(7,*) Rmin+i*meshsize, states2(j,i)**2
!         end do
!         CLOSE(unit=7)




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

         allocate(occup(1:num_stat), stat = ifail)
         if (ifail .ne. 0) STOP
          
        do j=0,points-1
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
!         print*, int(states2(5,i)+1),occup, Num_par
!
            density(j)= density(j) + (occup(i))*states2(7+j,i)**2/4/pi &
                        /posiR(j)**2
!
            der_u(j)= (states2(7+j+1,i)-states2(7+j,i))/meshsize
!
            kinden(j)= kinden(j)+occup(i)/posiR(j)**2/4/pi*    &
                     ((der_u(j)-states2(7+j,i)/posiR(i))**2+ &
                     states2(4,i)*(states2(4,i)+1)/posiR(i)**2*states2(7+j,i)**2)

            denN(j)= denN(j) + (occup(i))*states2(7+j,i)**2/4/pi &
                        /posiR(j)**2
!
            der_uN(j)= (states2(7+j+1,i)-states2(7+j,i))/meshsize
!
            kindenN(j)= kindenN(j)+occup(i)/posiR(j)**2/4/pi*    &
                     ((der_u(j)-states2(7+j,i)/posiR(i))**2+ &
                     states2(4,i)*(states2(4,i)+1)/posiR(i)**2*states2(7+j,i)**2)

!            if((j .eq. 1 .or. j .eq. 2) .and. int(states2(4,i)) .eq. 0) then
!            density(0)= density(0) +(occup(i))*states2(10,i)/4/pi &
 !                       /posiR(j)**2
!            denN(0)= denN(0) +(occup(i))*states2(10,i)/4/pi &
!                        /posiR(j)**2
!            end if
         end do
         end do
         
! HERE SOME COMMENTS
!loop on proton density
         do j=0,points-1
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
            density(j)= density(j) + (occup(i))*states2(7+j,i)**2/4/pi &
                        /posiR(j)**2
!
            der_u(j)= (states2(7+j+1,i)-states2(7+j,i))/meshsize
!
            kinden(j)= kinden(j)+occup(i)/posiR(j)**2/4/pi*    &
                     ((der_u(j)-states2(7+j,i)/posiR(i))**2+ &
                     states2(4,i)*(states2(4,i)+1)/posiR(i)**2*states2(7+j,i)**2)

            denZ(j)= denZ(j) + (occup(i))*states2(7+j,i)**2/4/pi &
                        /posiR(j)**2
!
            der_uZ(j)= (states2(7+j+1,i)-states2(7+j,i))/meshsize
!
            kindenZ(j)= kindenZ(j)+occup(i)/posiR(j)**2/4/pi*    &
                     ((der_u(j)-states2(7+j,i)/posiR(i))**2+ &
                     states2(4,i)*(states2(4,i)+1)/posiR(i)**2*states2(7+j,i)**2)
!         if((j .eq. 1 .or. j .eq. 2) .and. int(states2(4,i)) .eq. 0) then
!            density(0)= density(0) +(occup(i))*states2(10,i)/4/pi &
!                        /posiR(j)**2
!            denZ(0)= denZ(0) +(occup(i))*states2(10,i)/4/pi &
!                        /posiR(j)**2
!            end if
         end do
         end do

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
           write(10,*) Rmin+i*meshsize, density(i), denN(i), denZ(i)
         end do
         CLOSE(unit=10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!continue from here!!!!!!!!
!
!        Coulomb density energy
         do i=0, points
              E_cou2=E_cou2+ meshsize*posiR(i)*denZ(i)
         end do
!
         allocate(E_cou(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
         allocate(E_cou1(0:points), stat = ifail)
         if (ifail .ne. 0) STOP
!
!
         E_cou(:)=0.0d0
         E_cou1(:)=0.0d0
         do j=1,points
            do i=0, j
              E_cou1(j)=E_cou1(j)+ meshsize*(posiR(i)**2*denZ(i)/posiR(j) &
                      -posiR(i)*2*denZ(i))
            end do
            E_cou(j)=4.0d0*pi*esquare*(E_cou1(j)+E_cou2)
            print*, E_cou(j)
         end do
         t0=-1132.400d0
         t3=23610.40d0
         alpha=1.0d0
         
!         allocate(C_rho0(0:points), stat = ifail)
!         if (ifail .ne. 0) STOP
!         allocate(C_rho1(0:points), stat = ifail)
!         if (ifail .ne. 0) STOP   

!         C_rho0(:)=0.0d0
!         C_rho1(:)=0.0d0  
    

         tot_E=0.0d0
         do j=1,points
            Crho0=3.0d0/8.0d0*(-1132.400d0) + 3.0d0/48.0d0*23610.40d0* &
              density(j)**alpha
            Crho1=-1.0d0/8.0d0*(-1132.4000d0)-1.0d0/48.0d0*23610.40d0*&
              density(j)**alpha
            tot_E= tot_E + meshsize*(posiR(j))**2 &
                *(hb2m* kinden(j)+ &
             (Crho0)*density(j)**2+ (Crho1)*(denN(j)-denZ(j))**2+ &
              E_cou(j)*denZ(j)/2.0d0)
!         print*, tot_E
         end do
         
         print*, tot_E
!
         if (abs(tot_E-E_before) .lt. epsil) then
            HF_conv=1
            print*,"Everything converge"
         end if
         end do
            

!
         
         deallocate(trialwf)
         deallocate(pott)
      end program numerov
! end
!
!
!
 
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
! trialwf = trial wavefunction
! Etrial = energy of trial wavefunction
         subroutine wavef(points,meshsize, n_rad,l, Eup0, Edown0, epsil, pott, &
                 trialwf, Etrial) 
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
         REAL(kind=dp), DIMENSION(0:points), INTENT(OUT) :: trialwf
         REAL(kind=dp) :: pot, a, Vvalue
         REAL(kind=dp) , DIMENSION(0:points), INTENT(IN) :: pott
!
         !loading the potential
           Eup=Eup0
           Edown=Edown0

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
              trialwf(i)=0.0E0_dp!meshsize**(l+1)
           end do

           trialwf(imin)=meshsize**(l+1)
           imax=points-1

!loop until convergence
         converg=0
         do while (converg .ne. 1)
         Etrial=(Eup+Edown)/2


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
         if(trialwf(points-1)*trialwf(points) .lt. 0.0E0_dp) &
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



!======================================================================
! CALCULATING THE WAVEFUNCTION IN TWO PARTS: LEFT AND RIGHT PART

         subroutine parwf(points,meshsize, n_rad, l, numN, numZ, Etrial, &
                    epsil, pott, trialwf) 
         use variables
         implicit none
         REAL(kind=dp), INTENT(IN) :: Etrial, epsil
         REAL(kind=dp), INTENT(IN)  :: meshsize
         integer, INTENT(IN)  :: n_rad, l, numN, numZ
         integer, INTENT(IN):: points
         integer :: Nodecount , ifail, i, idir, imin, imax, converg, i_sur
         REAL(kind=dp) :: a1, a2, a3, normal, normalr, normall, rcentr
         REAL(kind=dp), DIMENSION(0:points), INTENT(OUT) :: trialwf
         REAL(kind=dp), DIMENSION(0:points) ::  triall, trialr
         REAL(kind=dp) :: pot, surf
         REAL(kind=dp), DIMENSION(0:points), INTENT(IN) :: pott
         REAL(kind=dp) :: derivr, derivl, contir, contil, coeff
!
! half left wavefunction (positive direction = 1)
           imin=1
           idir=1
           triall(0)=0.0E0_dp

           triall(imin)=meshsize**(l+1)
           imax=points-1  !because f(x+2h) is necessary

         do i=imin,imax,idir
           pot=(Etrial-pott(i))/hb2m !insert V if V=0
           a1= 2.0E0_dp*(1.0E0_dp-5.0E0_dp/12.0E0_dp*pot*meshsize**2)
           a2=(1.0E0_dp+1.0E0_dp/12.0E0_dp*pot*meshsize**2)
           a3=(1.0E0_dp+1.0E0_dp/12.0E0_dp*pot*meshsize**2)
           triall(i+idir)= (a1*triall(i)-a2*triall(i-idir))/a3
         end do

!
! half right wavefunction (negative direction = -1)

           imin=1
           idir=-1
           trialr(0)=0.0E0_dp
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
            coeff= contil/contir
            trialr(:)= coeff*trialr(:)
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
! Calculates the potential function
	! kpot= 0 no potential      
	! kpot= 1 square well of size a and deepth Vvalue
	! kpot= 2 Wood-Saxon potential
	! kpot= 3 Coulomb
	! kpot= 4 spin-orbit W potential
	! kpot= 5 Radial potential Eq. 2.10
	! Variables:
	! vr= potential as a result
	! i= point in a mesh where potential is calculated
	! meshsize = distance between two meshpoints
	! Rmax = right bound for potential
	! Rmin = left bound for potential
	! Vvalue = depth of infinite potential 
	! numN = number of neutrons
	! numZ = number of protons
	! rzero = a constant related to a radius of nucleus
	! rProt = proton radius
	! charge = 0 (neutron), 1 (proton)
	! l = quantum number l
	! jj = quantum number j*2 (=>integer)
         function potV(i,meshsize,Rmax,Rmin,a,kpot, Vvalue, &
                  numN, numZ,rzero,rProt,charge, l, jj) result(vr)
         use variables
         implicit none
         real (kind = dp) :: meshsize, Rmax, Rmin,a,Vvalue, rzero, bigR, ferre
         real (kind=dp) :: vr, pos, leftlim, rightlim, pos0, rProt
         real (kind=dp) :: Uq_pot, Wq_pot1,Wq_pot2, V_Cou, V_cfug
         integer :: i, kpot, numZ, numN, points, l, jj, charge
         ! symmetric square well
         pos=Rmin+i*meshsize
         points=(Rmax-Rmin)/meshsize !redundant
         bigR= rzero*(numN+numZ)**(1.0E0_dp/3.0E0_dp)

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
           leftlim=Rmin
           rightlim=a
           if(pos .ge. leftlim .and. pos .le. rightlim) then
             vr= Vvalue
           else
             vr= 0
           end if
!-------------------------------------
! WOODS-SAXON

         else if(kpot .eq. 2) then
             ferre=1/(1+exp((abs(pos-pos0)-bigR)/a))
             vr=-44.0192307692*ferre
 !            vr=-(51-33*(numN-numZ)/(numN+numZ))*ferre
!---------------------------------------
! COULOMB

         else if(kpot .eq. 3) then
             if(abs(pos-pos0) .le. rProt) then
             vr=numZ*esquare/2/rProt*(3-(pos/rProt)**2)         
             else if(abs(pos-pos0) .gt. rProt) then
             vr=numZ*esquare/rProt 
             end if 
!-------------------------------------------
! SPIN-ORBIT

         else if(kpot .eq. 4) then
             ferre=-0.44*exp((abs(pos-pos0)-bigR)/a) * rzero**2 &
                   /a/pos/(1+exp((abs(pos-pos0)-bigR)/a))**2
             vr=-(51-33*(numN-numZ)/(numN+numZ))*ferre 
!--------------------------------------------  
! RADIAL POTENTIAL

         else if(kpot .eq. 5) then
             ! Uq potential is the Wood-Saxon
             if(i .eq. 0.1E0_dp)then
                pos = meshsize
             end if
             Uq_pot=- 44.0192307692 &
                    /(1+exp((abs(pos)-bigR)/a))
!                   -(51-33*(numN-numZ)/(numN+numZ)) &
!                   /(1+exp((abs(pos-pos0)-bigR)/a))
             !jj=2j so jj is integer
             ferre=-0.44*exp((abs(pos)-bigR)/a) * rzero**2 &
                   /a/(pos)/(1+exp((abs(pos)-bigR)/a))**2
!
             Wq_pot1=ferre*(dble(jj)/2*(dble(jj)/2+1)-l*(l+1)-3/4)/(pos) 
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
               if(abs(pos) .le. rProt) then
                V_Cou=numZ*esquare/2/rProt*(3-(pos/rProt)**2)         
               else if(abs(pos) .gt. rProt) then
                V_Cou=numZ*esquare/rProt 
               end if
             vr= vr+V_Cou
             end if
         end if
!
         end if
         
end function potV
