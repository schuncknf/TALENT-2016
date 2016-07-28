      subroutine hfsolver(pr)
       use constants
       use lag_pol
       use pot
       use maths
       use ho
       use basis
      implicit none   
      integer:: lwmax
      double precision,parameter::lambda=1.d-8
      integer :: it, n , lda, ldvl, ldz,ldvr, info, i, j, k, l, lwork,liwork,il,iu,n_val
      double precision, allocatable :: hf(:,:), eigvecr(:,:), eigvecl(:,:)
      double precision, allocatable :: eigvalr(:), eigvall(:), eigvalold(:), work(:)
      double precision, allocatable :: rho(:,:), vpot(:,:,:,:), kin(:,:), gama(:,:),trho(:,:),hrho(:,:),gammarho(:,:)
      double precision, allocatable :: vpotp(:,:,:,:),vpotas(:,:,:,:)
      double precision, allocatable :: vpotpr(:,:,:,:),vpotasr(:,:,:,:),vpotr(:,:,:,:)
      integer,allocatable::nl(:),nr(:),nj(:),isupz(:),iwork(:),nocc_diag(:)
      double precision esum, rhosum, gamasum, vnorm,tr
      double precision::x,ri,rj,hfenergy,hf2body,kin_energy,tbme1,tbme2
      double precision,allocatable::rhob(:,:,:),kinb(:,:,:),gamab(:,:,:),t_mat(:,:)
      integer::n1,n2,n3,n4
      integer::l1,l2,l3,l4
      integer::m1,m2,m3,m4
      integer::j1,j2,j3,j4
      integer::i1,i2,i3,i4
      integer::bn,bl,bj,nbloc,bloc_index,red_i,ii,ij
      !external dgeev
      external DSYEVR
      logical::pr
      double precision:: ABSTOL,vl,vu
!
! --------------------------------------------------
!
      !n=nbase
      ABSTOL = 0.0
      !n=red_size
      nbloc = maxval(n_red)
      n=nbloc
      !n=2
      ! n=size_full
      !n=ntx
      lda = n
      ldvr = n
      ldvl = n
      ldz = n
      lwmax = 26*n
! ----- allocation of memory
! 
      allocate(work(26*n),iwork(10*n),isupz(2*n))
      allocate(hf(n,n),eigvecr(n,n),eigvecl(n,n))
      allocate(eigvalr(n),eigvall(n),eigvalold(n))
      allocate(t_mat(n,n))
      allocate(trho(n,n),hrho(n,n),rho(n,n),vpot(n,n,n,n),kin(n,n),gama(n,n),gammarho(n,n))
      allocate(vpotp(n,n,n,n),vpotas(n,n,n,n))
      allocate(vpotr(n,n,n,n))
      allocate(nocc_diag(n))
      allocate(rhob(2*occ_states,n,n),kinb(2*occ_states,n,n),gamab(2*occ_states,n,n))
! --------- two-body matrix elements and kinetic energy (to be calculated from subroutines)


vpotas=0.d0
write(*,*) "computing TBME",n
!!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(n,vpotas) SCHEDULE(DYNAMIC)
     do i = 1,n
      n1=i!-1
      do j = 1,n
       n2=j!-1
       do k = 1,n
        n3=k!-1
        do l = 1,n
         n4=l!-1
!           call tbme(n1,n2,n3,n4,vpotr(i,j,k,l),.false.,1)
           !call tbme(n1,n2,n4,n3,vpotpr(i,j,k,l),.false.,1)
 !          vpotas(i,j,k,l) = (vpotr(i,j,k,l))! + vpotpr(i,j,k,l))
        enddo !l
       enddo !k
      enddo !j
     enddo !i
!!$OMP END PARALLEL DO
!vpotas(1:nbase,1:nbase,1:nbase,1:nbase) = vpotasr(1:nbase,1:nbase,1:nbase,1:nbase)
  !         vpotas = 0.d0
write(*,*) "TBME's computed and antisymetrized"

! ---------- start of iteration loop

        t_mat = 0 
        bloc_index = 0
        do l = minval(l_red),maxval(l_red) 
          call t_bloc(n,l,t_mat)
          do j= 0,1
          bloc_index = bloc_index + 1
          write(*,*) "# Of block",bloc_index
          if (bloc_index .gt. occ_states +1 ) exit
          if (l == 0) then 
          bj = 1
          elseif (j == 0) then 
          bj = 2*l - 1 
          else 
          bj = 2*l + 1
          endif
!           do bn1 = 0,maxval(n_red)
!            do bn2 = 0,maxval(n_red) 
!               bloc(bloc_index,bn1,bn2) =  

      eigvecr = 0.d0
      do it = 1, maxit
         if(it .eq. 1) then ! initializing eigenfunctions and eigenvalues
            nocc_diag = nocc
            eigvecr = 0.d0
            do i = 0,nbloc
               red_i = tag_hf(i,l,bj)
               if(nocc(red_i) .gt. 0) then 
               eigvecr(i,i) = 1.d0!*nocc(red_i) ! first npart states occupied with npart particles
               endif
            enddo
         call compute_rho(rho,eigvecr,nbloc,l,bj)
         call compute_gamma(gama,vpotas,rho,nbloc,l,bj)
         endif


! --------- subroutines: (re)calculate rho and hf hamiltonian
       



         hf = t_mat! + gama
! --------- diagonalization of hamiltonian

      il = 1
      iu = n
      vl = 0
      vu = 0
      n_val = 0 
      info = 0
      lwork = 26*n
      liwork =10*n

      call dsyevr('V','I','U',n,hf,n,vl,vu,1,n,abstol,n_val,eigvalr,eigvecr,ldz, isupz, work, lwork, iwork,liwork, info)
!      write(*,*) "héhé",eigvalr
         if(info .ne. 0 ) stop 'problem in diagonalization'
! --------- check for convergence
      nocc_diag = 0
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


! ------- check out normalization of states

      do ij = 1, n
         vnorm= 0.d0
         do ii = 1, n
            vnorm = vnorm + eigvecr(ii,ij)*eigvecr(ii,ij)
         enddo 
         if(abs(vnorm-1.d0) .gt. 0.0001d0) stop 'problem in normalization'
      enddo
         write(*,*) "Interm-energy= ",sum(eigvalr)
         call compute_rho(rho,eigvecr,nbloc,l,bj)
         call compute_gamma(gama,vpotas,rho,nbloc,l,bj)
         rhob(bloc_index,1:nbloc,1:nbloc) = rho 
         kinb(bloc_index,1:nbloc,1:nbloc) = t_mat 
         gamab(bloc_index,1:nbloc,1:nbloc) = gama 


         enddo !j
      enddo !l
      enddo !it

!-------- hf energy

       hfenergy = 0.d0
       tr = 0.d0
       do i=1,occ_states+5
          hfenergy = hfenergy + trace(matmul(kinb(i,1:nbloc,1:nbloc),rhob(i,1:nbloc,1:nbloc)),nbloc)& 
 &                      + half*trace(matmul(gamab(i,1:nbloc,1:nbloc),rhob(i,1:nbloc,1:nbloc)),nbloc)
          tr = tr + trace(rhob(i,1:nbloc,1:nbloc),nbloc)
       enddo
       write(*,'(a,f16.9,a)') 'Block True Hartree-Fock Energy',hfenergy,' MeV'
       write(*,'(a,f16.9)') 'Block Particles Number',tr


       hrho=0.d0
       hrho = matmul(hf,rho)
       trho=0.d0
       trho = matmul(kin,rho)
       hfenergy = 0.d0
       kin_energy= 0.d0
       gammarho=0.d0
       gammarho = matmul(gama,rho)
       do i=1,n
         kin_energy = kin_energy + kin(i,i)*nocc(i)
         write(*,'(a,i3,a,f16.8)') 'e(',i,')= ',eigvalr(i)
       enddo
       do i=1,n
        do j=1,n
           hf2body = hf2body + dsqrt(dble(nocc(i)*nocc(j)))*vpotas(i,j,i,j)
         enddo
        enddo
        write(*,*) "2-Body Interaction Energy",hf2body
       write(*,*) "Kinetic energy",kin_energy
!       hfenergy = trace(trho,n) + trace(hrho,n)
       write(*,*) "Part Num",trace(rho,n)
!       trho = matmul(kin,rho)
       hfenergy = 0.d0
       do i=1,n
         hfenergy = hfenergy + (trho(i,i) + half*gammarho(i,i))
       enddo
!
!       hfenergy=0.5*hfenergy (?)
!
! ---- Petar

       hfenergy = 0.d0
       do i=1,n
         hfenergy = hfenergy + nocc(i)*half*(kin(i,i) + eigvalr(i))
       enddo
       write(*,'(a,f16.9,a)') 'True Hartree-Fock Energy2',hfenergy,' MeV'
! -------- Writing Outputs
open(22,file='hforsay.out')
write(22,*) "------ System ------"
write(22,'(a,i5)') 'Number of particles......',Npart
write(22,'(a,i5)') 'Oscillator Shell.........',n
write(22,*) "---- RESULTS ------"
write(22,'(a,f16.9,a)') 'Kinetic Energy..........',kin_energy,' MeV'
write(22,'(a,f16.9,a)') "Hartree-Fock Energy ....",hfenergy,' MeV'
write(22,'(a,f16.9)') "Particles Number ....",tr
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
