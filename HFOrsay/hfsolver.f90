
      program HFsolver

      IMPLICIT NONE   

!      USE constants  

      INTEGER, PARAMETER :: maxit= 1000, lambda=1.d-8, lwmax=1000
      INTEGER :: it, N, Npart, LDA, LDVL, LDVR, INFO, i, j, LWORK
      DOUBLE PRECISION, ALLOCATABLE :: hf(:,:),eigvecR(:,:),eigvecL(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: eigvalR(:),eigvalL(:),eigvalOLD(:),WORK(:)
      DOUBLE PRECISION :: esum

      EXTERNAL dgeev

      write(*,*) 'Type in dimension of HF matrix'

      read(*,*) N

      write(*,*) 'Type in number of particles'

      read(*,*) Npart

      LDA = N
      LDVR = N
      LDVL = N
 
      ALLOCATE(hf(1:N,1:N))
      ALLOCATE(eigvecR(1:N,1:N))
      ALLOCATE(eigvecL(1:N,1:N))
      ALLOCATE(eigvalR(N),eigvalL(N),eigvalOLD(N),WORK(LWMAX))

	
      do it = 1, maxit

	 if(it .eq. 1) then ! initializing eigenfunctions and eigenvalues
	    do i = 1, N
	       do j = 1, N
	          eigvecR(i,j) = 0.0
	       enddo
	       if(i .le. Npart) eigvecR(i,i) = 1.0
	       eigvalOLD(i) = 0.0
	    enddo
         endif

! subroutines: (re)calculate rho and HF hamiltonian

	LWORK = -1
        call DGEEV('N','V',N, hf, LDA, eigvalR, eigvalL, eigvecL, LDVL, &
                  eigvecR, LDVR, WORK, LWORK, INFO ) 
        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
        call DGEEV('N','V',N, hf, LDA, eigvalR, eigvalL, eigvecL, LDVL, &
                  eigvecR, LDVR, WORK, LWORK, INFO ) ! diagonalization of HF hamiltonian

	if(info .ne. 0 ) stop 'problem in diagonalization'

	esum = 0.0 ! check for convergence
	do i = 1,N 
	   esum = esum + abs(eigvalR(i) - eigvalOLD(i))
	enddo
	esum = esum/N

! check for convergence
	if(esum .lt. lambda) exit

! save old eigenvalues
        do i = 1, N
	   eigvalOLD(i) = eigvalR(i)
	enddo

      enddo !it


      DEALLOCATE(hf,eigvecR,eigvecL)
      DEALLOCATE(eigvalR,eigvalL,eigvalOLD,WORK)

      end program HFsolver		
