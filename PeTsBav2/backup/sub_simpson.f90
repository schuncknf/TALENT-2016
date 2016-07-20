!***************************************************
!*      this subroutine is simpson's formula       *
!*      N: the number of division                  *
!*      dx: axis step size                         *
!*      f: array of function                       *
!*      int: integrated value                      *
!***************************************************

subroutine simpson(N, dx, f, int)
implicit none
real(8) :: int, int1, int2, x, dx
integer :: N
real(8) :: f(0:N)
integer :: i
int = 0d0
int1 = 0d0
int2 = 0d0
do i=1,int(N/2)-1
   int1 = int1 + f(2*i)
end do
int1 = 2*int1
do i=1,int(N/2)
   int2 = int2 + f(2*i-1)
end do
int2 = 4*int2
int = int + f(0) + f(N) + int1 + int2
int = dx/3*int
end subroutine simpson

