      subroutine hfsolver(pr)
       use constants
       use lag_pol
       use pot
       use maths
       use ho
       use basis
      implicit none   
      integer:: lwmax,jf
      double precision,parameter::lambda=1.d-8
      integer :: it, n , lda, ldvl, ldz,ldvr, info, i, j, k, l, lwork,liwork,il,iu,n_val
      !double precision, allocatable :: hf(:,:), eigvecr(:,:), eigvecl(:,:)
      double precision, allocatable :: hf(:,:,:), eigvecr(:,:,:),eigvalold(:,:)
      !double precision, allocatable :: eigvalr(:), eigvall(:), eigvalold(:), work(:)
      double precision, allocatable :: eigvalr(:,:), work(:)
      double precision, allocatable :: rho(:,:,:), vpot(:,:,:,:,:), kin(:,:), gama(:,:,:),trho(:,:),hrho(:,:),gammarho(:,:)
      double precision, allocatable :: vpotas(:,:,:,:,:)
      double precision, allocatable :: vpotpr(:,:,:,:),vpotasr(:,:,:,:),vpotr(:,:,:,:)
      integer,allocatable::nl(:),nr(:),nj(:),isupz(:),iwork(:),nocc_diag(:)
      double precision esum, rhosum, gamasum, vnorm,tr
      double precision::x,ri,rj,hfenergy,hf2body,kin_energy,tbme1,tbme2
      double precision,allocatable::rhobo(:,:,:),kinbo(:,:,:),gamabo(:,:,:),t_mat(:,:,:),hfbo(:,:,:)
      double precision,allocatable::rhob(:,:,:),kinb(:,:,:),gamab(:,:,:),hfb(:,:,:)
      integer::n1,n2,n3,n4
      integer::l1,l2,l3,l4
      integer::m1,m2,m3,m4
      integer::j1,j2,j3,j4
      integer::i1,i2,i3,i4
      integer::bn,bl,bj,nbloc,bloc_index,red_i,ii,ij
      external DSYEVR
      logical::pr
      double precision:: ABSTOL,vl,vu
      double precision::hfold
!
! --------------------------------------------------
!
      !n=nbase
      ABSTOL = 0.0
      !n=red_size
      nbloc = maxval(n_red) + 1
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
      allocate(hf(0:2*occ_states,n,n),eigvecr(0:2*occ_states,n,n))
      allocate(eigvalr(0:2*occ_states,n))
      allocate(eigvalold(0:2*occ_states,n))
      allocate(t_mat(0:2*occ_states,n,n))
      allocate(trho(n,n),hrho(n,n),rho(0:2*occ_states,n,n),vpot(0:2*occ_states,n,n,n,n),kin(n,n))
      allocate(gama(0:2*occ_states,n,n),gammarho(n,n))
      allocate(vpotas(0:2*occ_states,n,n,n,n))
      allocate(vpotr(n,n,n,n))
      allocate(nocc_diag(n))
      allocate(rhobo(0:2*occ_states,n,n),kinbo(0:2*occ_states,n,n),gamabo(0:2*occ_states,n,n),hfbo(0:2*occ_states,n,n))
      allocate(rhob(0:2*occ_states,n,n),kinb(0:2*occ_states,n,n),gamab(0:2*occ_states,n,n),hfb(0:2*occ_states,n,n))
! --------- two-body matrix elements and kinetic energy (to be calculated from subroutines)



! ---------- start of iteration loop

  write(*,*) "Occ states",occ_states,"It. Number",maxit
  do it = 1,maxit
!        t_mat = 0 
! ------------------ Constructing Blocks        
        bloc_index = 0
        do l = minval(l_red),maxval(l_red) 
          if (l == 0) then 
          jf = 0
          else
          jf=1
          endif
          do j= 0,jf
          if (bloc_index .gt. occ_states - 1) then
          hf = 0.d0 
          exit
          endif 
          if (l == 0) then 
          bj = 1
          elseif (j == 0) then 
          bj = 2*l - 1 
          else 
          bj = 2*l + 1
          endif
        
!        eigvecr = 0.d0
         if(it .eq. 1) then ! initializing eigenfunctions and eigenvalues
            nocc_diag = nocc
            eigvecr = 0.d0
            eigvalold = 0.d0
            do i = 1,nbloc
               red_i = tag_hf(i-1,l,bj)
               if(nocc(red_i) .gt. 0) then 
               eigvecr(bloc_index,i,i) = 1.d0!*nocc(red_i) ! first npart states occupied with npart particles
               endif
            enddo
         endif
         call t_bloc(n,l,t_mat(bloc_index,:,:))
         call compute_rho(rho(bloc_index,:,:),eigvecr(bloc_index,:,:),nbloc,l,bj)
         call compute_gamma(gama(bloc_index,:,:),rho(bloc_index,:,:),nbloc,l,bj)
         hf(bloc_index,:,:) = t_mat(bloc_index,:,:)  + 0.d0*gama(bloc_index,:,:)
         !write(*,*) "Gama",bloc_index 
         !write(*,*) gama(bloc_index,:,:)
         !write(*,*)
         

! --------- subroutines: (re)calculate rho and hf hamiltonian
        rhobo(bloc_index,1:nbloc,1:nbloc) = rho(bloc_index,:,:)
        gamabo(bloc_index,1:nbloc,1:nbloc) = gama(bloc_index,:,:)
        kinbo(bloc_index,1:nbloc,1:nbloc) = t_mat(bloc_index,:,:)
        hfbo(bloc_index,1:n,1:n) = hf(bloc_index,:,:)




! --------- diagonalization of hamiltonian

      il = 1
      iu = n
      vl = 0
      vu = 0
      n_val = 0 
      info = 0
      lwork = 26*n
      liwork =10*n
      call dsyevr('V','I','U',n,hf(bloc_index,:,:),n,vl,vu,1,n,abstol,n_val,eigvalr(bloc_index,:),eigvecr(bloc_index,:,:),ldz, isupz, work, lwork, iwork,liwork, info)
         if(info .ne. 0 ) stop 'problem in diagonalization'
      call compute_rho(rho(bloc_index,:,:),eigvecr(bloc_index,:,:),nbloc,l,bj)
      call compute_gamma(gama(bloc_index,:,:),rho(bloc_index,:,:),nbloc,l,bj)
      call t_bloc(n,l,t_mat(bloc_index,:,:))
        rhob(bloc_index,1:nbloc,1:nbloc) = rho(bloc_index,:,:)
        gamab(bloc_index,1:nbloc,1:nbloc) = gama(bloc_index,:,:)
        kinb(bloc_index,1:nbloc,1:nbloc) = t_mat(bloc_index,:,:)
        hfb(bloc_index,1:n,1:n) = hf(bloc_index,:,:)
!         write(*,*) "Bloc,gammarho",bloc_index,matmul(gama(bloc_index,:,:),rho(bloc_index,:,:))
!         write(*,*) "Bloc,Rho",bloc_index,rho(bloc_index,:,:)
!         write(*,*) "Bloc,Gamma",bloc_index,gama(bloc_index,:,:)
!         write(*,*) "Bloc,RhoGamma",bloc_index,matmul(gama(bloc_index,:,:),rho(bloc_index,:,:))
         do i = 1, n
            eigvalold(bloc_index,i) = eigvalr(bloc_index,i)
         enddo

      
      bloc_index = bloc_index + 1
         enddo !j

      enddo !l
! ------------ End Diagonalization for all blocks



! --------- check for convergence


       hfenergy = 0.d0
       hfold    = 0.d0
       esum = 0.d0
       do i=0,occ_states-1
!      write(*,*) "Before"
!      write(*,*) hfbo(i,:,:)
!      write(*,*) "After"
!      write(*,*) hfb(i,:,:)
          hfenergy = hfenergy + trace(matmul(t_mat(i,1:nbloc,1:nbloc),rhob(i,1:nbloc,1:nbloc)),nbloc)& 
 &                      + half*trace(matmul(gamab(i,1:nbloc,1:nbloc),rhob(i,1:nbloc,1:nbloc)),nbloc)
          hfold = hfold + trace(matmul(kinbo(i,1:nbloc,1:nbloc),rhobo(i,1:nbloc,1:nbloc)),nbloc)& 
 &                      + half*trace(matmul(gamabo(i,1:nbloc,1:nbloc),rhobo(i,1:nbloc,1:nbloc)),nbloc)
!          !tr = tr + trace(rho(i,1:nbloc,1:nbloc),nbloc)
       enddo
       esum = abs(hfenergy - hfold) 






!      nocc_diag = 0
!         esum = 0.d0 
!      do j=0,occ_states
!         do i = 1,n 
!            esum = esum + abs(eigvalr(j,i) - eigvalold(j,i))
!         enddo
!       enddo
         esum = esum/n
         write(*,*) "iteration: ",it,"ediffi= ",esum
         if(esum .lt. lambda) exit ! calculation converged

! -------- save old eigenvalues


! ------- check out normalization of states

!      do ij = 1, n
!         vnorm= 0.d0
!         do ii = 1, n
!            vnorm = vnorm + eigvecr(ii,ij)*eigvecr(ii,ij)
!         enddo 
!         if(abs(vnorm-1.d0) .gt. 0.0001d0) stop 'problem in normalization'
!      enddo
!         write(*,*) "Interm-energy= ",sum(eigvalr)
      enddo !it

!-------- hf energy

       hfenergy = 0.d0
       tr = 0.d0
       do i=0,occ_states-1 
          hfenergy = hfenergy + trace(matmul(t_mat(i,1:nbloc,1:nbloc),rhob(i,1:nbloc,1:nbloc)),nbloc)& 
 &                      + half*trace(matmul(gamab(i,1:nbloc,1:nbloc),rhob(i,1:nbloc,1:nbloc)),nbloc)
          tr = tr + trace(rho(i,1:nbloc,1:nbloc),nbloc)
       enddo
       write(*,'(a,f16.9,a)') 'Block True Hartree-Fock Energy',hfenergy,' MeV'
       write(*,'(a,f16.9)') 'Block Particles Number',tr


!       hrho=0.d0
!       hrho = matmul(hf,rho)
!       trho=0.d0
!       trho = matmul(kin,rho)
!       hfenergy = 0.d0
!       kin_energy= 0.d0
!       gammarho=0.d0
!       gammarho = matmul(gama,rho)
!!       do i=1,n
!         kin_energy = kin_energy + kin(i,i)*nocc(i)
!         write(*,'(a,i3,a,f16.8)') 'e(',i,')= ',eigvalr(i)
!       enddo
!       do i=1,n
!        do j=1,n
!           hf2body = hf2body + dsqrt(dble(nocc(i)*nocc(j)))*vpotas(i,j,i,j)
!         enddo
!        enddo
!        write(*,*) "2-Body Interaction Energy",hf2body
!       write(*,*) "Kinetic energy",kin_energy
!!       hfenergy = trace(trho,n) + trace(hrho,n)
!       write(*,*) "Part Num",trace(rho,n)
!!       trho = matmul(kin,rho)
!       hfenergy = 0.d0
!       do i=1,n
!         hfenergy = hfenergy + (trho(i,i) + half*gammarho(i,i))
!       enddo
!
!       hfenergy=0.5*hfenergy (?)
!
! ---- Petar

!       hfenergy = 0.d0
!       do i=1,n
!         hfenergy = hfenergy + nocc(i)*half*(kin(i,i) + eigvalr(i))
!       enddo
!       write(*,'(a,f16.9,a)') 'True Hartree-Fock Energy2',hfenergy,' MeV'
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
!   write(22,*) 'energy = ', eigvalr(j),' MeV'
!   write(22,*) 'Wave-function:'
do i = 1, n
!      write(22,*) eigvecr(i,j)
   enddo
write(22,*) '------------------'
      enddo
endif

close(22)


! -------- deallocate memory

      deallocate(hf,eigvecr)
      deallocate(eigvalr,work)
      deallocate(rho,vpot,kin,gama)

      end subroutine
!
