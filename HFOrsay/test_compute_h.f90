program test_compute_rho

implicit none
external compute_rho, compute_h, compute_gamma
! use constants

complex, dimension (:,:), allocatable :: test_D
double precision, dimension (:,:), allocatable :: test_rho,test_gamma,test_h,test_t
double precision, dimension (:,:,:,:), allocatable :: test_TBME
integer :: test_dim
integer :: i,j,flag

write(*,*), "please select matrix dimension"
read(*,*), test_dim

allocate(test_D(test_dim,test_dim))
allocate(test_rho(test_dim,test_dim))
allocate(test_gamma(test_dim,test_dim))
allocate(test_h(test_dim,test_dim))
allocate(test_t(test_dim,test_dim))

allocate(test_TBME(test_dim,test_dim,test_dim,test_dim))

test_D(:,:) = 0.d0
test_t(:,:) = 0.d0
test_TBME(:,:,:,:) = 1.d0

write(*,*), "real : 1 , Complex : 2, Null : 3"
read(*,*) flag

if (flag .eq. 1) then
   do i=1,test_dim
    do j=1,test_dim
       if (i == j) test_D(i,j) = (1.d0,0.d0)
    enddo
    enddo
endif
if (flag .eq. 2) then
   do i=1,test_dim
    do j=1,test_dim
       if (i == j) test_D(i,j) = (0.d0,1.d0)
    enddo
    enddo
endif
if (flag .eq. 3) then
   do i=1,test_dim
    do j=1,test_dim
       if (i == j) test_D(i,j) = (0.d0,0.d0)
    enddo
    enddo
endif

call compute_rho(test_rho,test_D,test_dim)
call compute_gamma(test_gamma,test_TBME,test_rho,test_dim)
call compute_h(test_h,test_t,test_gamma,test_dim)


write(*,*) "rho values"
do i=1,test_dim
  do j=1,test_dim
      write (*,*) test_rho(i,j)
  enddo
write(*,*) "    "
enddo

write(*,*) "hamiltonian values"
do i=1,test_dim
  do j=1,test_dim
      write (*,*) test_h(i,j)
  enddo
write(*,*) "    "
enddo

deallocate(test_D)
deallocate(test_rho)
deallocate(test_h)
deallocate(test_t)
deallocate(test_TBME)
deallocate(test_gamma)

end program test_compute_rho
