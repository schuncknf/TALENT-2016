program hforsay
use constants
use lag_pol
use maths
use pot
use ho

implicit none
double precision::ml,x
double precision::test
integer::i,nl,n1,n2
double precision::test11(4)
double precision::start,finish


call reader()
call cpu_time(start)
call lag_roots(ngauss,0.5d0,.true.)
call hfsolver(.true.)
call cpu_time(finish)
print '("Cpu-time = ",f6.3," Seconds")',finish-start

end program
