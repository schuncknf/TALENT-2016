module maths
contains
recursive function fac(n) result (fak)
implicit none
integer::n,fak
if (n == 0) then
    fak = 1
else
    fak = n*fac(n-1)
endif
end function fac
recursive function ffac(n) result (ffak)
implicit none
integer::n,ffak
if (n .gt. -2 .and. n .le. 0) then
    ffak = 1
else
    ffak = n*ffac(n-2)
endif
end function ffac
end module maths





