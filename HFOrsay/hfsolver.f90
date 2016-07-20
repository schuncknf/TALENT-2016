
      program HFsolver
!
      IMPLICIT NONE   
!
!      USE constants  
!
      INTEGER, PARAMETER :: maxit= 1000, lambda=1.d-8, lwmax=1000
      INTEGER :: it, N, Npart, LDA, LDVL, LDVR, INFO, i, j, k, l, LWORK
      DOUBLE PRECISION, ALLOCATABLE :: hf(:,:), eigvecR(:,:), eigvecL(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: eigvalR(:), eigvalL(:), eigvalOLD(:), WORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: rho(:,:), vpot(:,:,:,:), kin(:,:), gama(:,:)
      DOUBLE PRECISION esum, rhosum, gamasum, vnorm
!
      EXTERNAL dgeev
!
! --------------------------------------------------
!
      write(*,*) 'Type in dimension of basis:'
!
      read(*,*) N
!
      write(*,*) 'Type in number of particles:'
!
      read(*,*) Npart
!
      LDA = N
      LDVR = N
      LDVL = N
!
! ----- allocation of memory
! 
      ALLOCATE(hf(1:N,1:N),eigvecR(1:N,1:N),eigvecL(1:N,1:N))
      ALLOCATE(eigvalR(N),eigvalL(N),eigvalOLD(N),WORK(LWMAX))
!
      ALLOCATE(rho(1:N,1:N),vpot(1:N,1:N,1:N,1:N),kin(1:N,1:N),gama(1:N,1:N))
!
! --------- two-body matrix elements and kinetic energy (to be calculated from subroutines)
!
      do i = 1, N
 	 do j = 1, N
            do k = 1, N
 	       do l = 1, N
!
	          if(i .ne. j .and. j .ne. k .and. k .ne. l) then
                      vpot(i,j,k,l) = 0.0
                  else 
                      vpot(i,j,k,l) = real(0.5/i)
	          endif
!
               enddo
             enddo
          enddo
      enddo
!
      do i = 1, N ! Npart?
	 do j = 1, N
	    kin(i,j) = 0.0
	 enddo
	 kin(i,i) = 0.5*i
      enddo

! ---------- start of iteration loop
	
      do it = 1, maxit

	 if(it .eq. 1) then ! initializing eigenfunctions and eigenvalues
	    do i = 1, N
	       do j = 1, N
	          eigvecR(i,j) = 0.0
	       enddo
	       if(i .le. Npart) eigvecR(i,i) = 1.0 ! first Npart states occupied with Npart particles
	       eigvalOLD(i) = 0.0 
	    enddo
         endif

! --------- subroutines: (re)calculate rho and HF hamiltonian

         call compute_rho(rho,eigvecR,N)

         call compute_gamma(gama,vpot,rho,N)

         call compute_h(hf,kin,gama,N)

! --------- diagonalization of hamiltonian
	       
	 LWORK = -1
         call DGEEV('N','V',N, hf, LDA, eigvalR, eigvalL, eigvecL, LDVL, &
                  eigvecR, LDVR, WORK, LWORK, INFO ) 
         LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
         call DGEEV('N','V',N, hf, LDA, eigvalR, eigvalL, eigvecL, LDVL, &
                  eigvecR, LDVR, WORK, LWORK, INFO ) 

	 if(info .ne. 0 ) stop 'problem in diagonalization'

! --------- check for convergence

	 esum = 0.0 
	 do i = 1,N 
	    esum = esum + abs(eigvalR(i) - eigvalOLD(i))
	 enddo
	 esum = esum/N

	 if(esum .lt. lambda) exit ! calculation converged

! -------- save old eigenvalues

         do i = 1, N
	    eigvalOLD(i) = eigvalR(i)
	 enddo

      enddo !it

! ------- check out normalization of states

      do j = 1, N
	 vnorm= 0
	 do i = 1, N
	    vnorm = vnorm + eigvecR(i,j)*eigvecR(i,j)
         enddo 
	 if(abs(vnorm-1.0) .gt. 0.0001) stop 'problem in normalization'
      enddo

! -------- print out eigenfunctions and eigenvalues

      do j = 1, N
	 write(*,*) 'state number = ', j
         write(*,*) 'energy = ', eigvalR(j)
	 write(*,*) 'wave function:'
         do i = 1, N
	    write(*,*) eigvecR(i,j)
	 enddo
	 write(*,*) '------------------'
      enddo

! -------- deallocate memory

      DEALLOCATE(hf,eigvecR,eigvecL)
      DEALLOCATE(eigvalR,eigvalL,eigvalOLD,WORK)
      DEALLOCATE(rho,vpot,kin,gama)

      end program HFsolver		
!
