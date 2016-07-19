program test_compute_rho

implicit none
external compute_rho
! use constants

complex, dimension (:,:), allocatable :: test_D
double precision, dimension (:,:), allocatable :: test_rho
integer :: test_dim
integer :: i,j,flag

write(*,*), "please select matrix dimension"
read(*,*), test_dim

allocate(test_D(test_dim,test_dim))
allocate(test_rho(test_dim,test_dim))

test_D(:,:) = 0.d0

write(*,*), "real : 1 , Complex : 2"
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

call compute_rho(test_D,test_rho,test_dim)

write(*,*) "rho values"
do i=1,test_dim
  do j=1,test_dim
      write (*,*) test_rho(i,j)
  enddo
  write(*,*) "    "
enddo

deallocate(test_D)
deallocate(test_rho)

end program test_compute_rho
