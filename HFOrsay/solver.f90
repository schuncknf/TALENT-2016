subroutine solver()
use basis
use maths
use constants
use pot
implicit none
integer::i,l,j,k,jf,q
integer::sizebloc
integer::info,lwork
integer::it
double precision,parameter::eps=1.d-10 
double precision,allocatable::rho(:,:,:,:),eigvec(:,:,:,:),gama(:,:,:,:)
double precision,allocatable::t(:,:,:),h(:,:,:,:),eigval(:,:,:),work(:)
double precision,allocatable::conv(:,:,:)
double precision::diff,hfenergy,partnum

sizebloc=maxval(n_red) + 1
allocate(rho(lmin:lmax,jmin:jmax,sizebloc,sizebloc))
allocate(gama(lmin:lmax,jmin:jmax,sizebloc,sizebloc))
allocate(eigvec(lmin:lmax,jmin:jmax,sizebloc,sizebloc))
allocate(h(lmin:lmax,jmin:jmax,sizebloc,sizebloc))
allocate(t(lmin:lmax,sizebloc,sizebloc))
allocate(work(26*sizebloc))
allocate(conv(lmin:lmax,jmin:jmax,maxit))
allocate(eigval(lmin:lmax,jmin:jmax,sizebloc))


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
  hfenergy = hfenergy + trace(matmul(t(l,:,:),rho(l,j,:,:)),sizebloc) + half*trace(matmul(gama(l,j,:,:),rho(l,j,:,:)),sizebloc)
  partnum = partnum + trace(rho(l,j,:,:),sizebloc)
  enddo !k
 enddo !l
diff = diff/(lmax*2)
write(*,*) "Iteration. ",it,"diff= ",diff
 if (it .gt. 1 .and. diff .lt. eps) exit
enddo !it

write(*,*) "Hartree-Fock Energy",hfenergy
write(*,*) "Particles Number",partnum
end subroutine














