program hforsay
!use omp_lib
use constants
use lag_pol
use maths
use pot
use ho
use basis

implicit none
double precision::test
double precision::test11(4)
double precision::start,finish
integer::j


call reader()
!write(*,*) "tbme",tbme_ext(17,18,89,18)
!call cpu_time(start)
call external_basis()
!call external_tbme()
call filled_number()
write(*,*) "maxval",maxval(n_red)+1
write(*,*) "# Of occ states",occ_states
read(*,*)
write(*,'(a)') "State, n , l , j , nocc"
do j=1,red_size
write(*,'(a,5i4)') "State: ",j,n_red(j),l_red(j),j_red(j),nocc(j)
write(*,*) 
enddo
!start = omp_get_wtime()
call lag_roots(ngauss,0.5d0,.false.)
call hfsolver(.true.)
!finish = omp_get_wtime()
!call cpu_time(finish)
!print '("Real cpu-time = ",f6.3," Seconds")',finish-start

end program
