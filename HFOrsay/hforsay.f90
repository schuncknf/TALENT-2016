program hforsay
use constants
use lag_pol
use maths
use pot
use ho
use basis
implicit none
double precision::start,finish
integer::j


call fancy(.true.)
call reader()
call cpu_time(start)
call external_basis()
call external_tbme(.true.)
call filled_number()
write(*,*) "There is",occ_states,"occupied states"
!call lag_roots(ngauss,0.5d0,.false.)
call solver(.true.)
call cpu_time(finish)
print '("Real cpu-time = ",f6.3," Seconds")',finish-start
end program
