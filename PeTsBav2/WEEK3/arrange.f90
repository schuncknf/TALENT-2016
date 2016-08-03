subroutine arrange(E, dim, func, N)
implicit none
real(8) :: E(dim,4), box(4)
real(8) :: func(0:N,1:dim), func_tmp(0:N)
real(8) :: number(dim)
integer :: i, j, dim, N

do i=1,dim-1
   do j=1,dim-i
      number(j) = E(j,1)
      number(j+1) = E(j+1,1)
      if (number(j) > number(j+1)) then
         func_tmp(:) = func(:,j)
         box(:) = E(j,:)
         func(:,j) = func(:,j+1)
         E(j,:) = E(j+1,:)
         func(:,j+1) = func_tmp(:)
         E(j+1,:) = box(:)
      end if
   end do
end do
end subroutine arrange

