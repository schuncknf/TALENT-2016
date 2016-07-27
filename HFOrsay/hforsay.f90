program hforsay
use omp_lib
use constants
use lag_pol
use maths
use pot
use ho
use basis

implicit none
double precision::ml,x
double precision::test
double precision::test11(4)
double precision::start,finish


call reader()
!write(*,*) "tbme",tbme_ext(17,18,89,18)
!call cpu_time(start)
call external_basis()
write(*,*) n_ext
read(*,*)
start = omp_get_wtime()
call lag_roots(ngauss,0.5d0,.true.)
call hfsolver(.true.)
finish = omp_get_wtime()
!call cpu_time(finish)
print '("Real cpu-time = ",f6.3," Seconds")',finish-start

end program
