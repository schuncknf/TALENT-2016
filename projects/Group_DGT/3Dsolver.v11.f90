     program numerov
      use variables
      implicit none
!     program to solve 1D Schroedinger equation
!      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
      REAL(kind=dp) :: Eup, Edown, epsil
      REAL(kind=dp) :: Rmin, Rmax, meshsize, a, Vvalue, rzero
      integer :: Nodemax, kpot, N_max, N_maxN, N_maxZ, N_act, Num_par
      integer :: Num_par_max
      REAL(kind=dp) :: Etrial, Eexp
      integer :: points, ifail, converg, i,j, charge, l, jj, nstat, jmin, jmax
      integer :: Nodecount, numN, numZ, strval, n_rad, occup, st_pos, minN
      integer :: minZ, nstatN, nstatZ, num_stat
      REAL(kind=dp) :: a1, a2, a3, pot, normal, deltae, pos
      REAL(kind=dp) :: potV, rProt, valueR0, den_int, denN_int, denZ_int
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: trialwf, posiR
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: pott, density, denN, denZ
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:,:) :: states, states2
      REAL(kind=dp) , ALLOCATABLE, DIMENSION(:) :: ordEneN, ordEneZ
      CHARACTER(LEN=3):: char_nstat
      CHARACTER(LEN=5) :: string
      CHARACTER(LEN=3):: part
!==========     set initial values     ==============================
      Eup=100.E0_dp
      Edown=-100.E0_dp
      epsil=1.E-6_dp
!      Nodemax=5
      Rmin=0.1E0_dp
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
      N_maxN=8
      N_maxZ=8
!==================================================================     
!
      points=int((Rmax-Rmin)/meshsize)

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
 !               if(l .eq. 0) then
 !                jmin=1
 !                jmax=1
 !               else
 !                jmin=2*l-1
 !                jmax=2*l+1
 !               end if
               jmin=abs(2*l-1)
               do jj= jmin,2*l+1,2
               nstat=nstat+1

!               print*,nstat
               
                
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
!renormalization

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
! output number of nodes excluded a possible node in the last point of the box
! normalized wavefunction
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
            states(7+i,nstat)=trialwf(i)**2
         end do
                end do !loop j
31            continue
              end do ! loop l
            end do !  loop n_rad
         end do ! N_max
         end do ! charge             


         ordEneN(:)=states(1,1:nstatN)
         ordEneZ(:)=states(1,nstatN+1:num_stat)

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
!
!         OPEN(UNIT=7, FILE='u(r)_state'//char_nstat//'.dat', status='unknown')
!         
!         write(7,*) "#Radial number=", states2(2,i)
!         write(7,*) "#Energy= ", states2(1,i)
!         do j=7,points
!           write(7,*) Rmin+i*meshsize, states2(j,i)
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
          
        do j=2,points
!loop on neutron density
           charge=0
           N_max=N_maxN
           Num_par_max=numN
           Num_par=0
         do i=1, nstatN
         if(Num_par+int(states2(5,i)+1) .le. Num_par_max) then
           occup=(states2(5,i)+1)
         else if (Num_par .lt. Num_par_max) then
           occup=Num_par_max-Num_par
         else if (Num_par .gt. Num_par_max) then
           occup=0
         end if
         Num_par= Num_par+int(states2(5,i)+1)
!         print*, int(states2(5,i)+1),occup, Num_par
!
            density(j)= density(j) + (dble(occup))*states2(7+j,i)/4/pi &
                        /posiR(j)**2
!            print*, posiR(j), density(j), i

            denN(j)= denN(j) + (dble(occup))*states2(7+j,i)/4/pi  &
                        /posiR(j)**2
            if((j .eq. 1 .or. j .eq. 2) .and. int(states2(4,i)) .eq. 0) then
            density(0)= density(0) +(dble(occup))*states2(10,i)/4/pi &
                        /posiR(j)**2
            denN(0)= denN(0) +(dble(occup))*states2(10,i)/4/pi &
                        /posiR(j)**2
            end if
         end do
         end do
         
!loop on proton density
         do j=1,points
           charge=1
           N_max=N_maxZ
           Num_par_max=numZ
           Num_par=0
         do i=nstatN+1, num_stat
         if(Num_par+(states2(5,i)+1) .le. Num_par_max) then
           occup=(states2(5,i)+1)
         else if (Num_par .lt. Num_par_max) then
           occup=Num_par_max-Num_par
         else if (Num_par .gt. Num_par_max) then
           occup=0
         end if
         Num_par= Num_par+(states2(5,i)+1)
!         print*, int(states2(5,i)+1),occup
!
         density(j)= density(j) + (dble(occup))*states2(7+j,i)/4/pi &
                        /posiR(j)**2
         denZ(j)= denZ(j) +  (dble(occup))*states2(7+j,i)/4/pi &
                        /posiR(j)**2
         if((j .eq. 1 .or. j .eq. 2) .and. int(states2(4,i)) .eq. 0) then
            density(0)= density(0) +(dble(occup))*states2(10,i)/4/pi &
                        /posiR(j)**2
            denZ(0)= denZ(0) +(dble(occup))*states2(10,i)/4/pi &
                        /posiR(j)**2
            end if
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
! kinetic density
!         do i=0, points-1,
!           derdens(i)= (density(i+1)-density(i))/meshsize
!           derdenN(i)= (denN(i+1)-denN(i))/meshsize
!           derdenZ(i)= (denZ(i+1)-denZ(i))/meshsize
!         end do
!           derdens(points)= 0.0d0
!           derdenN(points)= 0.0d0
!           derdenZ(points)= 0.0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!continue from here!!!!!!!!
!         do i=0, points
!         end do

!
         
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
!
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
!           if(l .eq. 0) then
!           trialwf(0)=meshsize**(l+1)
!           end if
           trialwf(imin)=meshsize**(l+1)
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
 !          print '("Convergence obtained for wavefunction")'
!           print *, "with energy= ", Etrial
! expected energy from the formula Eq. 1.14 in manuale_hf.pdf
!          if (kpot .eq. 0) then 
!            Eexp= ((Node+1)*pi)**2*hb2m/(Rmax-Rmin)**2
!            print *, "Expected energy= ", Eexp
!          end if
           converg=1
         end if
41       continue
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
         REAL(kind=dp) :: a1, a2, a3, normal, normalr, normall, rcentr
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
!            if (l .gt. 0) then
!             rcentr=sqrt(Etrial/(l*(l+1))/hb2m)
!             imin=int(rcentr/meshsize)
!           else
             imin=1
!           end if
!
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
!         allocate(trialr(0:points), stat = ifail)
!         if (ifail .ne. 0) STOP
!
!           if (l .gt. 5) then
!             rcentr=sqrt(Etrial/(l*(l+1))/hb2m)
!             imin=int(rcentr/meshsize)
!           else
             imin=1
 !          end if
!
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
         bigR= rzero*(numN+numZ)**(1.0E0_dp/3.0E0_dp)
!
         if(i .eq. 0) then
           vr=0.0d0
         else if (i .ne. 0) then
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
