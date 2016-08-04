subroutine LU(N, A, b, x, err)
implicit none

complex(8),parameter :: zI=(0d0,1d0) ! unit of imaginary part
complex(8) :: A(1:N,1:N), b(1:N), x(1:N) 
complex(8) :: y(1:N), L_U(1:N,1:N), tmp(1:N,1:N), tmp2(1:N)
integer :: N, i, j, kk, err

err = 0

do j=1,N
   L_U(1,j) = A(1,j)
!   write(*,'(a6,i2,a2,2e16.6)') 'LU( 1,' , j, ')=', L_U(1,j)
end do
if (L_U(1,1)==0d0) then
   err = -1
   write(*,*) 'Ax=b cannot be solved.'
   goto 99
end if

do i=2,N
   do j=1,i-1
      tmp = 0d0
      if (j>=2) then
         do kk=1,j-1
            tmp(i,j) = tmp(i,j) + L_U(i,kk)*L_U(kk,j)
!            write(*,*) tmp(i,j)
         end do
      end if
      L_U(i,j) = (A(i,j) - tmp(i,j))/L_U(j,j)
!      write(*,'(a3,i2,a1,i2,a2,2e16.6)') 'LU(', i, ',', j, ')=', L_U(i,j)
   end do
   do j=i,N
      tmp = 0d0
      if (i>=2) then
         do kk=1,i-1
            tmp(i,j) = tmp(i,j) + L_U(i,kk)*L_U(kk,j)
!            write(*,*) tmp(i,j)
         end do
      end if
      L_U(i,j) = A(i,j) - tmp(i,j)
!      write(*,'(a3,i2,a1,i2,a2,2e16.6)') 'LU(', i, ',', j, ')=', L_U(i,j)
   end do
   if (L_U(i,i)==0d0) then
      err = -1
      write(*,*) 'Ax=b cannot be solved.'
      goto 99
   end if
end do

y(1) = b(1)
!write(*,'(a6,2e16.6)') 'y( 1)=', y(1)
do i=2,N
   tmp2 = 0d0
   do kk=1,i-1
      tmp2(i) = tmp2(i) + L_U(i,kk)*y(kk)
   end do
   y(i) = b(i) - tmp2(i)
!   write(*,'(a2,i2,a2,2e16.6)') 'y(', i, ')=', y(i)
end do

x(N) = y(N)/L_U(N,N)
do i=2,N
   tmp2 = 0d0
   do kk=1,i-1
      tmp2(N-i+1) = tmp2(N-i+1) + L_U(N-i+1,N-i+1+kk)*x(N-i+1+kk)
   end do
   x(N-i+1) = (y(N-i+1) - tmp2(N-i+1))/L_U(N-i+1,N-i+1)
end do

99 continue
end subroutine
