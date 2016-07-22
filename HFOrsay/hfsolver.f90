      subroutine hfsolver()
       use constants
       use lag_pol
       use pot
       use maths
       use ho
!
      IMPLICIT NONE   
!
!      USE constants  
!
      INTEGER, PARAMETER :: lwmax=1000
      double precision,parameter::lambda=1.d-8
      INTEGER :: it, N , LDA, LDVL, LDVR, INFO, i, j, k, l, LWORK
      DOUBLE PRECISION, ALLOCATABLE :: hf(:,:), eigvecR(:,:), eigvecL(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: eigvalR(:), eigvalL(:), eigvalOLD(:), WORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: rho(:,:), vpot(:,:,:,:), kin(:,:), gama(:,:),trho(:,:),hrho(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: vpotm(:,:,:,:),vpotp(:,:,:,:),vpotas(:,:,:,:)
      DOUBLE PRECISION esum, rhosum, gamasum, vnorm
      double precision::x,ri,rj,hfenergy,hf2body,kin_energy
      double precision,allocatable::temp2(:,:)
      integer::n1,n2,n3,n4
      integer::i1,i2,i3,i4
!      external::compute_rho,compute_h,compute_gamma
!
      EXTERNAL dgeev
!
! --------------------------------------------------
!
      !write(*,*) 'Type in dimension of basis:'
!
      !read(*,*) N
      !call reader()
      !write(*,*) "Basis size",nbase
      n=nbase
      !write(*,*) "Npart",npart
      !write(*,*) "Iteration",maxit
      
!
!
!
      LDA = N
      LDVR = N
      LDVL = N
!
! ----- allocation of memory
! 
      ALLOCATE(hf(n,n),eigvecR(n,n),eigvecL(n,n))
      ALLOCATE(eigvalR(N),eigvalL(N),eigvalOLD(N),WORK(LWMAX))
!
      ALLOCATE(trho(n,n),hrho(n,n),rho(n,n),vpot(n,n,n,n),kin(n,n),gama(n,n))
      allocate(vpotm(n,n,n,n),vpotp(n,n,n,n),vpotas(n,n,n,n))
      allocate(temp2(1,n+1))
!----  Laguerre Mesh
    ! call lag_roots(n,0.5d0,.true.)
     !call gausslag(n,1,2,x)
! --------- two-body matrix elements and kinetic energy (to be calculated from subroutines)
!
!
        kin =0.d0
      do i = 1, N+1 ! Npart?
         kin(i,i) = (2.d0*(i-1)+1.5d0)*ama*2.d0/(bosc**2)
      enddo
vpot=0.d0
vpotp=0.d0
vpotas=0.d0
write(*,*) "Computing TBMEs"
     do i = 1,n
      n1=i-1
      do j = 1,n
       n2=j-1
       do k = 1,n
        n3=k-1
        do l = 1,n
         n4=l-1
           call tbme(n1,n2,n3,n4,vpot(i,j,k,l),.true.)
           call tbme(n1,n2,n4,n3,vpotp(i,j,k,l),.false.)
           !     vpotp(i,j,k,l)=vpot(i,j,k,l)
           !vpotm(i,j,k,l)=vpot(i,j,l,k)
           !vpotas(i,j,k,l) = (vpot(i,j,k,l) + vpotp(i,j,l,k))
           !vpotas(i,j,k,l) = vpot(i,j,k,l)+vpotp(i,j,k,l) 
           !vpotas(i,j,l,k) = vpot(i,j,k,l)+vpotp(i,j,k,l) 
           !vpotas(j,i,k,l) = vpot(i,j,k,l)+vpotp(i,j,k,l) 
           !vpotas(j,i,l,k) = vpot(i,j,k,l)+vpotp(i,j,k,l) 
           write(14,*) n1,n2,n3,n4,vpotas(i,j,k,l)
        enddo !l
           write(14,*) 
       enddo !k
           write(14,*) 
      enddo !j
           write(14,*) 
     enddo !i
           vpotas = vpot + vpotp 
           !vpotas = 0.d0
write(*,*) "TBMEs computed and antisymetrized"

! ---------- start of iteration loop

      do it = 1, maxit

         if(it .eq. 1) then ! initializing eigenfunctions and eigenvalues
            do i = 1, N
               do j = 1, N
                  eigvecR(i,j) = 0.d0
               enddo
               if(i .le. Npart) eigvecR(i,i) = 11.d0 ! first Npart states occupied with Npart particles
               eigvalOLD(i) = 0.d0 
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

         esum = 0.d0 
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
         vnorm= 0.d0
         do i = 1, N
            vnorm = vnorm + eigvecR(i,j)*eigvecR(i,j)
         enddo 
         if(abs(vnorm-1.d0) .gt. 0.0001d0) stop 'problem in normalization'
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

!-------- HF Energy


       !call compute_rho(rho,eigvecR,N)
       !call compute_h(hf,kin,gama,N)

       hrho=0.d0
       hrho = matmul(hf,rho)
       trho=0.d0
       trho = matmul(kin,rho)
       hfenergy = 0.d0
       do i=1,Npart
         hfenergy = hfenergy + half*trho(i,i) + half*hrho(i,i)
         !hfenergy = hfenergy + half*trho(i,i) + half*hrho(i,i)
       enddo
       write(*,*) 'Hartree-Fock Energy-1',hfenergy
       
       call sort(n,eigvalR)
       do i=1,n
         write(*,'(a,i3,a,f16.8)') 'E(',i,')= ',eigvalR(i)
       enddo
       hfenergy = 0.d0
       hf2body = 0.d0
       kin_energy = 0.d0
       do i=1,Npart
        hfenergy = hfenergy  + eigvalR(i)
        kin_energy = kin_energy + kin(i,i)
        do j=1,Npart
           hf2body = hf2body + vpotas(i,j,i,j)
        enddo
       enddo
       hfenergy = hfenergy  + kin_energy
       write(*,*) 'Hartree-Fock Energy',hfenergy*half
       write(*,*) "Kinetic energy",kin_energy


! -------- deallocate memory

      DEALLOCATE(hf,eigvecR,eigvecL)
      DEALLOCATE(eigvalR,eigvalL,eigvalOLD,WORK)
      DEALLOCATE(rho,vpot,kin,gama)

      end subroutine
!
