subroutine solver(pr)
use basis
use maths
use constants
use pot
use bcs
implicit none
logical::pr
integer::i,l,j,k,jf,q
integer::sizebloc
integer::info,lwork
integer::it
double precision,parameter::eps=1.d-11 
double precision,allocatable::rho(:,:,:,:),eigvec(:,:,:,:),gama(:,:,:,:)
double precision,allocatable::t(:,:,:),h(:,:,:,:),eigval(:,:,:),work(:)
double precision,allocatable::conv(:,:,:),v2(:),esp(:),qesp(:)
double precision::diff,hfenergy,partnum,kinenergy,hfenergybcs

sizebloc=maxval(n_red) + 1
allocate(rho(lmin:lmax,jmin:jmax,sizebloc,sizebloc))
allocate(gama(lmin:lmax,jmin:jmax,sizebloc,sizebloc))
allocate(eigvec(lmin:lmax,jmin:jmax,sizebloc,sizebloc))
allocate(h(lmin:lmax,jmin:jmax,sizebloc,sizebloc))
allocate(t(lmin:lmax,sizebloc,sizebloc))
allocate(work(26*sizebloc))
allocate(conv(lmin:lmax,jmin:jmax,maxit))
allocate(eigval(lmin:lmax,jmin:jmax,sizebloc))
allocate(v2(red_size),esp(red_size),qesp(red_size))


lwork = 26*sizebloc
conv = 0.d0
do it = 1,maxit
 do l = lmin,lmax
  if (l == 0) then
     jf = 0
  else
     jf = 1
  endif
  do k=0,jf
  if (l == 0) then 
    j = 1
   elseif (k==0) then
    j = 2*l - 1
   else
    j = 2*l+1
  endif
  call compute_rho(rho(l,j,:,:),eigvec(l,j,:,:),sizebloc,l,j,it)
    if (it .eq. 1) then  ! Initialisation Kinetic Matrix
      call t_bloc(sizebloc,l,t(l,:,:))
    endif !it
  enddo !j
 enddo !l
 do l = lmin,lmax
  if (l == 0) then
     jf = 0
  else
     jf = 1
  endif
  do k=0,jf
  if (l == 0) then 
    j = 1
   elseif (k==0) then
    j = 2*l - 1
   else
    j = 2*l+1
  endif

  call compute_gamma(gama(l,j,:,:),rho,sizebloc,l,j)
  h(l,j,:,:) = t(l,:,:) + gama(l,j,:,:)
  eigvec(l,j,:,:) = h(l,j,:,:)
  info = 0
  call dsyev('V','U',sizebloc,eigvec(l,j,:,:),sizebloc,eigval(l,j,:),work,lwork,info)
  if (info .ne. 0) stop "Diagonalization Failed"
  do q = 1,sizebloc
  conv(l,j,it) = conv(l,j,it) + eigval(l,j,q)
  enddo
  conv(l,j,it) = conv(l,j,it)/sizebloc
  enddo !j
 enddo !l

hfenergy = 0.d0
kinenergy = 0.d0
partnum = 0.d0
diff = 0.d0
 do l = lmin,lmax
  if (l == 0) then
     jf = 0
  else
     jf = 1
  endif
  do k=0,jf
  if (l == 0) then 
    j = 1
   elseif (k==0) then
    j = 2*l - 1
   else
    j = 2*l+1
  endif
  if (it .gt. 1) then
  diff = diff + abs(conv(l,j,it) - conv(l,j,it-1))
endif
  kinenergy = kinenergy + trace(matmul(t(l,:,:),rho(l,j,:,:)),sizebloc)
  hfenergy = hfenergy + trace(matmul(t(l,:,:),rho(l,j,:,:)),sizebloc) + half*trace(matmul(gama(l,j,:,:),rho(l,j,:,:)),sizebloc)
  partnum = partnum + trace(rho(l,j,:,:),sizebloc)
  do q = 1,sizebloc
        esp(tag_hf(q-1,l,j)) = eigval(l,j,q)
        if (nocc(tag_hf(q-1,l,j)) .ne. 0.d0) then
        v2(tag_hf(q-1,l,j)) = nocc(tag_hf(q-1,l,j))
        else
        v2(tag_hf(q-1,l,j)) = 0.d0
        endif
  enddo !q
  enddo !k
 enddo !l
diff = diff/(lmax*2)
write(*,*) "Iteration. ",it,"diff= ",diff
 if (it .gt. 1 .and. diff .lt. eps) exit
enddo !it

if (flagbcs .eq. 1) call pairing(red_size,esp,v2,.false.)
write(*,*) "Hartree-Fock Energy",hfenergy
if (flagbcs .eq. 1) then
partnum = two*sum(v2,red_size)
hfenergybcs = 0.d0
do k=1,red_size
 qesp(k) = esp(k) - x2 
 hfenergybcs = hfenergybcs + (two*qesp(k)*v2(k))
enddo
 hfenergybcs = hfenergybcs - gap**2/g_pair
endif
write(*,*) "Particles Number",partnum
open(22,file='hforsay.out')
write(22,*) "------ System ------"
write(22,'(a,i5)') 'Number of particles......',Npart
write(22,*) "---- RESULTS ------"
write(22,'(a,f16.9,a)') "Hartree-Fock Energy ....",hfenergy,' MeV'
write(22,'(a,f16.9)') "Particles Number......",partnum
if (flagbcs .eq. 1) then
write(22,'(a,f16.9)') "Pairing Gap......",gap
write(22,'(a,f16.9,a)') "Pairing Energy......",hfenergy - hfenergybcs,'Mev'
write(22,'(a,f16.9,a)') "Total Energy HF-BCS......",hfenergybcs,'MeV' 
endif
write(22,*)
if (pr) then
write(22,*) "--------- Single Particle Energies --------"
  write(22,'(a)') "   n  l  j          Energy(MeV)     Occ. Num"
  do k=1,red_size
  write(22,'(a,3i3,f18.6,a,f12.4)') ' ',n_red(k),l_red(k),j_red(l),esp(k),' MeV',v2(k)
  enddo
 write(22,*) "-------------------------------------------"
endif
close(22)
if(iplot .eq. 1) then
open(33,file='occ.dat')
call sorteigv(red_size,esp,v2)
do k=1,red_size
  write(33,'(i3,2f18.12)') k,esp(k),v2(k)
enddo
close(33)
endif


deallocate(rho,h,t,gama,work,conv)
end subroutine














