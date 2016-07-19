module HartreeFock
  use :: types
  use :: variables
  implicit none
contains

  subroutine Initialize_HF
    implicit none
    integer :: i
!    if(allocated(D_mat)) deallocate(D_mat)
    allocate(D_mat(1:Nsize,1:Nsize))
    allocate(rho_mat(1:Nsize,1:Nsize))
    allocate(h_mat(1:Nsize,1:Nsize))
    allocate(Gamma_mat(1:Nsize,1:Nsize))
    allocate(E_values(1:Nsize))
    allocate(E_prev(1:Nsize))
    E_prev = 0
    D_mat = 0
    do i = 1,nsize
       D_mat(i,i) = 1
    enddo
  end subroutine Initialize_HF

  subroutine Construct_rho
    implicit none
    integer :: i,j,k
    real(dp) :: D
    do i = 1,Nsize
       do j = i,Nsize
          D = 0
          do k = 1,Nsize
             D = D + D_mat(i,k)*D_mat(j,k)
          enddo
          rho_mat(i,j) = D
          rho_mat(j,i) = D
       enddo
    enddo
  end subroutine Construct_rho

  subroutine Construct_gamma
    implicit none
    integer :: i,j,k,l
    real(dp) :: gamma 
    do i = 1,Nsize
       do j = i,Nsize
          gamma = 0
          do k = 1,Nsize
             do l = 1,Nsize
                gamma = gamma + v_mat(i,l,j,k)*rho_mat(k,l)
             enddo
          enddo
          gamma_mat(i,j) = gamma
          gamma_mat(j,i) = gamma
       enddo
    enddo
  end subroutine Construct_gamma

  subroutine Diagonalize_h 
    implicit none
    real(dp),  dimension(1:3*nsize-1) :: Work
    integer :: lwork,info
    lwork = 3*nsize-1
    D_mat = h_mat
    call dsyev('V','U',Nsize,D_mat,Nsize,E_Values,Work,lwork,info)
  end subroutine Diagonalize_h


end module HartreeFock
