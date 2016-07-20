      program HFsolver
       use constants
      ! use lag
       use lag_pol
       use pot
       use maths
       use ho
!
      IMPLICIT NONE   
!
!      USE constants  
!
      INTEGER, PARAMETER :: maxit= 1000, lwmax=1000
      double precision,parameter::lambda=1.d-8
      INTEGER :: it, N, Npart, LDA, LDVL, LDVR, INFO, i, j, k, l, LWORK
      DOUBLE PRECISION, ALLOCATABLE :: hf(:,:), eigvecR(:,:), eigvecL(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: eigvalR(:), eigvalL(:), eigvalOLD(:), WORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: rho(:,:), vpot(:,:,:,:), kin(:,:), gama(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: vpotm(:,:,:,:),vpotp(:,:,:,:),vpotas(:,:,:,:)
      DOUBLE PRECISION esum, rhosum, gamasum, vnorm
      double precision::x,ri,rj
      double precision,allocatable::temp2(:,:)
!      external::compute_rho,compute_h,compute_gamma
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
      allocate(vpotm(n,n,n,n),vpotp(n,n,n,n),vpotas(n,n,n,n))
      allocate(temp2(1,n+1))
!----  Laguerre Mesh
     call lag_roots(n,0.5d0,.true.)
     !call gausslag(n,1,2,x)
! --------- two-body matrix elements and kinetic energy (to be calculated from subroutines)
!
!
      do i = 1, N ! Npart?
	 do j = 1, N
	    kin(i,j) = 0.0
	 enddo
	 kin(i,i) = (2.d0*i+1.5d0)*ama*2.d0/(bosc**2)
      enddo

do i=0,400
ri = 0.1*i
 do j=0,400
 rj = 0.1*j
write(11,*) ri,rj,potential(ri,rj,v0r,kr)
write(12,*) ri,rj,minnesota(ri,rj,v0r,kr)
enddo
write(11,*)
enddo


write(*,*) "Computing TBMEs"
     do i = 1,n
      do j = 1,n
       do k = 1,n
        do l = 1,n
           call tbme((n+1)/2,i,j,k,l,vpot(i,j,k,l))
           vpotp(i,j,k,l)=vpot(i,j,k,l)
           vpotm(i,j,k,l)=vpot(i,j,l,k)
           write(14,*) i,j,k,l,vpot(i,j,k,l)
        enddo !l
           write(14,*) 
       enddo !k
           write(14,*) 
      enddo !j
           write(14,*) 
     enddo !i
     vpotas = vpotp + vpotm
write(*,*) "TBMEs computed and antisymetrized"

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

         call compute_gamma(gama,vpotas,rho,N)

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
         write(*,*) "Iteration: ",it,"ediffi= ",esum

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
