      subroutine hfsolver(pr)
       use constants
       use lag_pol
       use pot
       use maths
       use ho
      implicit none   
      integer, parameter :: lwmax=1000
      double precision,parameter::lambda=1.d-8
      integer :: it, n , lda, ldvl, ldvr, info, i, j, k, l, lwork
      double precision, allocatable :: hf(:,:), eigvecr(:,:), eigvecl(:,:)
      double precision, allocatable :: eigvalr(:), eigvall(:), eigvalold(:), work(:)
      double precision, allocatable :: rho(:,:), vpot(:,:,:,:), kin(:,:), gama(:,:),trho(:,:),hrho(:,:)
      double precision, allocatable :: vpotm(:,:,:,:),vpotp(:,:,:,:),vpotas(:,:,:,:)
      double precision esum, rhosum, gamasum, vnorm,tr
      double precision::x,ri,rj,hfenergy,hf2body,kin_energy,tbme1,tbme2
      integer::n1,n2,n3,n4
      integer::i1,i2,i3,i4
      external dgeev
      logical::pr
!
! --------------------------------------------------
!
      n=nbase
      lda = n
      ldvr = n
      ldvl = n
! ----- allocation of memory
! 
      allocate(hf(n,n),eigvecr(n,n),eigvecl(n,n))
      allocate(eigvalr(n),eigvall(n),eigvalold(n),work(lwmax))
!
      allocate(trho(n,n),hrho(n,n),rho(n,n),vpot(n,n,n,n),kin(n,n),gama(n,n))
      allocate(vpotm(n,n,n,n),vpotp(n,n,n,n),vpotas(n,n,n,n))
! --------- two-body matrix elements and kinetic energy (to be calculated from subroutines)


        call kinetic(kin,n)
vpot=0.d0
vpotas=0.d0
write(*,*) "computing TBME"
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(n,vpotas,vpotp) SCHEDULE(DYNAMIC)
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
           vpotas(i,j,k,l) = (vpot(i,j,k,l) + vpotp(i,j,k,l))
        enddo !l
       enddo !k
      enddo !j
     enddo !i
!$OMP END PARALLEL DO
!     do i=1,n
!      do j=1,n
!       do k=1,n
!        do l=1,n        
         !vpotas(i,j,k,l) = vpot(i,j,k,l) + vpotp(i,j,k,l)
         ! vpotas(i,j,l,k) = vpot(i,j,k,l)+vpotp(i,j,k,l) 
         ! vpotas(j,i,k,l) = vpot(i,j,k,l)+vpotp(i,j,k,l) 
         ! vpotas(j,i,l,k) = vpot(i,j,k,l)+vpotp(i,j,k,l) 
           !vpotas(i,j,l,k) = vpot(i,j,k,l)+vpotp(i,j,k,l) 
!           write(14,*) i-1,j-1,k-1,l-1,vpotas(i,j,k,l)
!        enddo
!           write(14,*) 
!       enddo
!           write(14,*) 
!      enddo
!           write(14,*) 
!     enddo
      !     vpotas = 0.d0
write(*,*) "TBME's computed and antisymetrized"

! ---------- start of iteration loop

      do it = 1, maxit

         if(it .eq. 1) then ! initializing eigenfunctions and eigenvalues
            do i = 1, n
               do j = 1, n
               if(i .le. n) eigvecr(i,i) = 1.d0 ! first npart states occupied with npart particles
                  eigvecr(i,j) = 0.d0
               enddo
               eigvalold(i) = 0.d0 
            enddo
         endif

! --------- subroutines: (re)calculate rho and hf hamiltonian
       

         call compute_rho(rho,eigvecr,n)
         call compute_gamma(gama,vpotas,rho,n)
         hf = kin + gama
! --------- diagonalization of hamiltonian
               
         lwork = -1
         call dgeev('n','v',n, hf, lda, eigvalr, eigvall, eigvecl, ldvl, &
                  eigvecr, ldvr, work, lwork, info ) 
         lwork = min( lwmax, int( work( 1 ) ) )
         call dgeev('n','v',n, hf, lda, eigvalr, eigvall, eigvecl, ldvl, &
                  eigvecr, ldvr, work, lwork, info ) 

         if(info .ne. 0 ) stop 'problem in diagonalization'

! --------- check for convergence

         esum = 0.d0 
         do i = 1,n 
            esum = esum + abs(eigvalr(i) - eigvalold(i))
         enddo
         esum = esum/n
         write(*,*) "iteration: ",it,"ediffi= ",esum

         if(esum .lt. lambda) exit ! calculation converged

! -------- save old eigenvalues

         do i = 1, n
            eigvalold(i) = eigvalr(i)
         enddo

      enddo !it

! ------- check out normalization of states

      do j = 1, n
         vnorm= 0.d0
         do i = 1, n
            vnorm = vnorm + eigvecr(i,j)*eigvecr(i,j)
         enddo 
         if(abs(vnorm-1.d0) .gt. 0.0001d0) stop 'problem in normalization'
      enddo

! -------- print out eigenfunctions and eigenvalues


!-------- hf energy


       call compute_rho(rho,eigvecr,n)
       call compute_gamma(gama,vpotas,rho,n)
       hf = kin + gama
       call sorteigv(n,eigvalr,eigvecr)
       hrho=0.d0
       hrho = matmul(hf,rho)
       trho=0.d0
       trho = matmul(kin,rho)
       hfenergy = 0.d0
       do i=1,n
         write(*,'(a,i3,a,f16.8)') 'e(',i,')= ',eigvalr(i)
       enddo
       hfenergy = 0.d0
       hf2body = 0.d0
       kin_energy = 0.d0
       do i=1,npart/2
        hfenergy = hfenergy  + 2.d0*eigvalr(i)
        kin_energy = kin_energy + kin(i,i)
       enddo
       do i=1,npart/2
        do j=1,npart/2
           hf2body = hf2body + 2.d0*vpotas(i,j,i,j)
         enddo
        enddo
       hfenergy = hfenergy  - hf2body*half
       write(*,*) 'hartree-fock energy',hfenergy
       write(*,*) "kinetic energy",kin_energy
       hfenergy = 0.d0
       do i=1,n
         hfenergy = hfenergy + trho(i,i) + hrho(i,i)
       enddo
       write(*,'(a,f16.9,a)') 'True Hartree-Fock Energy',hfenergy,'MeV'
! -------- Writing Outputs
open(22,file='hforsay.out')
write(22,*) "------ System ------"
write(22,'(a,i5)') 'Number of particles......',Npart
write(22,'(a,i5)') 'Oscillator Shell.........',n
write(22,*) "---- RESULTS ------"
write(22,'(a,f16.9,a)') 'Kinetic Energy..........',kin_energy,' MeV'
write(22,'(a,f16.9,a)') "Hartree-Fock Energy ....",hfenergy,' MeV'
write(22,*)
if (pr) then
write(22,*) "---- Spectrum and Wave-Function ----"
do j = 1, n
   write(22,*) 'state number = ', j
   write(22,*) 'energy = ', eigvalr(j),' MeV'
   write(22,*) 'Wave-function:'
do i = 1, n
      write(22,*) eigvecr(i,j)
   enddo
write(22,*) '------------------'
      enddo
endif

close(22)


! -------- deallocate memory

      deallocate(hf,eigvecr,eigvecl)
      deallocate(eigvalr,eigvall,eigvalold,work)
      deallocate(rho,vpot,kin,gama)

      end subroutine
!
