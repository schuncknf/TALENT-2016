module Minnesota
  use :: types
  use :: variables
  implicit none
contains

  subroutine Initialize_Minnseota
    implicit none
    integer :: i
    allocate(t_mat(1:Nsize,1:Nsize))
    allocate(v_mat(1:Nsize,1:Nsize,1:Nsize,1:Nsize))
    allocate(n_ho(1:Nsize))
    allocate(l_ho(1:Nsize))
    t_mat = 0
    v_mat = 0
    l_ho = 0
    do i = 1,nsize
       n_ho(i) = i-1
       t_mat(i,i) = hw*(2*n_ho(i) + l_ho(i) + 1.5_dp)
    enddo
  end subroutine Initialize_Minnseota

end module Minnesota
