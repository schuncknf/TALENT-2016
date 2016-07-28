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
    allocate(D_prev(1:Nsize,1:Nsize))
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
    D_prev = D_mat
  end subroutine Initialize_HF

  subroutine Construct_rho
    implicit none
    integer :: i,j,k
    real(dp) :: D
    rho_mat = 0._dp
    do i = 1,Nsize
       do j = i,Nsize
          if(l_hf(i).ne.l_hf(j).or.j_hf(i).ne.j_hf(j)) cycle
          D = 0
          do k = 1,3 !3 is the number of occupied states for 8 particles
             D = D + (j_hf(k)+1)*D_mat(i,k)*D_mat(j,k)
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
    gamma_mat = 0
    do i = 1,Nsize
       do j = i,Nsize
          if(l_hf(i).ne.l_hf(j).or.j_hf(i).ne.j_hf(j)) cycle
          gamma = 0
          do k = 1,Nsize
             do l = 1,Nsize
                if(l_hf(k).ne.l_hf(l).or.j_hf(k).ne.j_hf(l)) cycle
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
    call dsyev('V','U',Nsize,D_mat,Nsize,E_Values,Work,lwork,info)
  end subroutine Diagonalize_h

  function Trace_product(A,B) result(TrAB)
    implicit none
    real(dp), dimension(:,:) :: A,B
    real(dp) :: TrAB
    integer :: N,i,j
    N = size(A,1)
    TrAB = 0
    do i = 1,nsize
       do j = 1,nsize
          TrAB = TrAB + A(i,j)*B(j,i)
       enddo
    enddo
  end function Trace_product

  function Trace(A) result(TrA)
    implicit none
    real(dp), dimension(:,:) :: A
    real(dp) :: TrA
    integer :: N,i
    N = size(A,1)
    TrA = 0
    do i = 1,nsize
       TrA = TrA + A(i,i)
    enddo
  end function Trace


end module HartreeFock
