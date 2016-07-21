module HartreeFock
  use :: types
  use :: variables
  implicit none
contains

  subroutine Initialize_HF
    implicit none
    integer :: i
    if(allocated(D_mat)) then
       deallocate(D_mat,rho_mat,h_mat,Gamma_mat,E_values,E_prev)
    endif
    allocate(D_mat(1:Nsize,1:Nsize))
    allocate(rho_mat(1:Nsize,1:Nsize))
    allocate(h_mat(1:Nsize,1:Nsize))
    allocate(Gamma_mat(1:Nsize,1:Nsize))
    allocate(E_values(1:Nsize))
    allocate(E_prev(1:Nsize))
    E_prev = 0
    E_values = 1
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
          do k = 1,Nparticles/2
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
    lwork = 30*nsize-1
    D_mat = h_mat
!    write(*,*) D_mat
    call dsyev('V','U',Nsize,D_mat,Nsize,E_Values,Work,lwork,info)
!    write(*,*) E_values
  end subroutine Diagonalize_h

  function Trace_product(A,B) result(TrAB)
    implicit none
    real(dp), dimension(:,:) :: A,B
    real(dp) :: TrAB
    integer :: N,i,j
    N = size(A,1)
    TrAB = 0
    do i = 1,N
       do j = 1,N
          TrAB = TrAB + A(i,j)*B(j,i)
       enddo
    enddo
  end function Trace_product


end module HartreeFock
