!***************************************************
!* this subroutine is simpson's formula *
!* N: the number of division *
!* dx: axis step size *
!* f: array of function *
!* res: integrated value *
!***************************************************
subroutine simpson(N, dx, f, res)
implicit none
real(8) :: res, even, odd, dx
integer :: N
real(8) :: f(0:N)
integer :: i
res = 0d0
odd = 0d0
even = 0d0
do i=1,N-1, 2
odd = odd + f(i)
end do
odd = 4.*odd
do i=2,N-2, 2
even = even + f(i)
end do
even = 2.*even
res = f(0) + f(N) + odd + even
res = dx/3.*res
end subroutine simpson

