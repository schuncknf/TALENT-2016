!---------------------------------------------------------------------
!Module: HartreeFock
!---------------------------------------------------------------------
!> Contains functions and subroutines for the Hartree-Fock self-
!! consistent part of the calculation, include building the density
!! matrix \f$\rho_{\mu\nu}\f$ and the single particle potential
!! \f$\Gamma_{\alpha\beta}\f$, and diagonalizing the Hamiltonian to
!! extract its eigenvectors and eigenvalues (see pg. 2 of
!! HF_truncated_v2.pdf).
!---------------------------------------------------------------------
module HartreeFock
  use :: types
  use :: variables
  implicit none
contains

  subroutine read_input()
    implicit none
    open(100,file='HF_input.dat')
    read(100,*)
    read(100,*) Nparticles
    read(100,*)
    read(100,*)
    read(100,*) Maxit
    read(100,*)
    read(100,*)
    read(100,*) k_fermi
    read(100,*)
    read(100,*)
    read(100,*) orbitals_file
    read(100,*)
    read(100,*)
    read(100,*) elements_file
    read(100,*)
    read(100,*)
    read(100,*) output_file
    read(100,*)
    read(100,*)
    read(100,*) density_file
    read(100,*)
    read(100,*)
    read(100,*) type_of_calculation
    if(trim(type_of_calculation).eq.'truncated') then
       read(100,*)
       read(100,*)
       read(100,*) Nsize
       truncated = .true.
       nsize = nsize + 1
    else
       truncated = .false.
    endif
    close(100)
    if(trim(type_of_calculation).ne.'truncated'.and.&
         trim(type_of_calculation).ne.'spherical'.and.&
         trim(type_of_calculation).ne.'LDA'.and.&
         trim(type_of_calculation).ne.'DME') then
       write(*,*) 'Unrecognized type of calculation ಠ_ಠ'
       write(*,*) type_of_calculation
       write(*,*) 'please check HF_input.dat file'
       write(*,*) 'stoping now'
       stop
    endif
    if(trim(type_of_calculation).eq.'LDA'.or.&
         trim(type_of_calculation).eq.'DME') then
       approximated_rho = .true.
    else
       approximated_rho = .false.
    endif
  end subroutine read_input

!> Allocates the arrays which will be used in the Hartree-Fock
!! calculation, and initializes the matrix \f$D_{\mu i}=\delta{\mu i}\f$
!! (see eqns. 1, 5 of HF_truncated_v2.pdf).
  subroutine Initialize_HF
    implicit none
    integer :: i
    if(allocated(D_mat)) then
       deallocate(D_mat,rho_mat,h_mat,Gamma_mat,E_values,E_prev)
    endif
    allocate(D_mat(1:Nsize,1:Nsize))
!    allocate(D_prev(1:Nsize,1:Nsize))
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
!    D_prev = D_mat
  end subroutine Initialize_HF

!> Constructs the density matrix 
!! \f$\rho_{\mu\nu}=\sum_{i=1}^ND_{\mu i}D^*_{\nu i}\f$.
  subroutine Construct_rho
    implicit none
    integer :: i,j,k
    real(dp) :: D
    rho_mat = 0._dp
    do i = 1,Nsize
       do j = i,Nsize
          if(l_hf(i).ne.l_hf(j).or.j_hf(i).ne.j_hf(j)) cycle
          D = 0
          do k = 1,Noccupied
             D = D + (j_hf(k)+1)*D_mat(i,k)*D_mat(j,k)
          enddo
          rho_mat(i,j) = D
          rho_mat(j,i) = D
       enddo
    enddo
  end subroutine Construct_rho

!> Constructs the single-particle potential
!! \f$\Gamma_{\alpha\beta}=\sum_{\mu\nu}v_{\alpha\nu\beta\mu}\rho{\mu\nu}\f$
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

!> Diagonalizes a matrix, returning its eigenvalues and eigenvectors, by
!! calling the LAPACK routine dsyev. Unless I'm mistaken, eigenvalues are
!! returned in the array \f$\texttt{E\_Values}\f$, and eigenvectors in
!! the array \f$\texttt{D\_mat}\f$.
  subroutine Diagonalize_h 
    implicit none
    real(dp),  dimension(1:3*nsize-1) :: Work
    integer :: lwork,info
    lwork = 30*nsize-1
    D_mat = h_mat
    call dsyev('V','U',Nsize,D_mat,Nsize,E_Values,Work,lwork,info)
  end subroutine Diagonalize_h

!> Computes the trace of the matrix produce \f$Tr(AB)\f$.
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

!> Computes the trace of a square matrix \f$Tr(A)\f$
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
