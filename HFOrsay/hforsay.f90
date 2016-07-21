program hforsay
use constants
use lag_pol
use maths
use pot
use ho

implicit none


call reader()
call lag_roots(nbase,0.5d0,.true.)
call hfsolver()





end program
