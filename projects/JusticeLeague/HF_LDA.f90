!---------------------------------------------------------------------
!Module: LDA
!---------------------------------------------------------------------
!> Contains functions and subroutines which allow you to approximate the
!! potential using the Local-Density Approximation or the Density
!! Matrix Expansion.
!---------------------------------------------------------------------
module LDA
  use :: types
  use :: variables
  use :: HO_Quadrature
  implicit none
  real(dp), private, parameter :: V0R=200.00_dp, muR=1.487_dp !< Parameters for the
  real(dp), private, parameter :: V0s= 91.85_dp, mus=0.465_dp !< Minnesota potential
  real(dp), private, parameter :: hw=10._dp, hbar=197.327_dp, M_n=939.57_dp
  real(dp), private, parameter :: h2_2m = 20.736209412_dp
  real(dp), private, parameter :: b_ho=sqrt(2*h2_2m/hw)
contains

!> Computes the density in the coordinate basis in the Local Density
!! Approximation [see eqn. 60 of HF_extensions.pdf].
  function rho_LDA(r) result(rho)
    implicit none
    real(dp), intent(in) :: r
    real(dp) :: rho
    real(dp) :: phi,xi, Anl,Anpl
    real(dp), dimension(0:5,0:2) :: L_nl
    integer :: i,nj,lj,j,ji,nk,lk,jj,jk,k
    xi = r/b_ho
    call LaguerreALL(5,2,xi**2,L_nl)
    rho = 0
    do i = 1,Noccupied!3
       phi = 0
       ji = j_hf(i)
       do j = 1,Nsize
          nj = n_hf(j)
          lj = l_hf(j)
          jj = j_hf(j)
          do k = 1,Nsize
             nk = n_hf(k)
             lk = l_hf(k)
             jk = j_hf(k)
             if(lk.ne.lj.or.jj.ne.jk) cycle
             Anl =  HO_Normalization(nj,lj)
             Anpl = HO_Normalization(nk,lk)
             phi = phi + (D_mat( i,j)*D_mat( i,k)*0.1_dp + &
                          D_prev(i,j)*D_prev(i,k)*0.9_dp)*&
                          xi**(lj+lk)*exp(-xi**2)*&
                          L_nl(nj,lj)*L_nl(nk,lk)*Anl*Anpl/(b_ho**3)
          enddo
       enddo
       rho = rho + (ji+1)*phi
    enddo
    rho = rho/(4*pi)
  end function rho_LDA

!> Calls the function \f$\texttt{gamma\_{LDA}(n,n',l)}\f$ to calculate
!! the matrix elements \f$\Gamma_{ij}, \Gamma_{ji}$ (see
!! HF_Extensions.pdf, eqn. 58).
  subroutine calculate_gamma_LDA 
    implicit none
    integer :: i,j, ni,nj,li,lj,ji,jj
    real(dp) :: gamma
    gamma_mat = 0
    do i = 1,nsize
       ni = n_hf(i)
       li = l_hf(i)
       ji = j_hf(i)
       do j = i,nsize
          nj = n_hf(j)
          lj = l_hf(j)
          jj = j_hf(j)
          if(li.ne.lj.or.ji.ne.jj) cycle
          gamma = gamma_LDA(ni,nj,li)
          gamma_mat(i,j) = gamma
          gamma_mat(j,i) = gamma
       enddo
    enddo
    
  end subroutine calculate_gamma_LDA

!> Performs the integral in HF_Extensions.pdf, eqn. 58)
  function gamma_LDA(n,np,l) result(gamma)
    implicit none
    integer, intent(in) :: n,np,l
    real(dp) :: gamma
    real(dp) :: alpha,xi,wi,An,Anp,Ln,Lnp,rho,Ivc
    integer :: i
    save Ivc
    if(calc_Ivc) then
       Ivc = IntegralVc()
       calc_ivc = .false.
    endif
    An  = HO_Normalization(n ,l)
    Anp = HO_Normalization(np,l)
    gamma = 0
    do i = 1,N_quad
       call LaguerreL( n,l+0.5_dp,x_quad(i),Ln)
       call LaguerreL(np,l+0.5_dp,x_quad(i),Lnp)
       rho = rho_quad(i)
       gamma = gamma + w_quad(i)*x_quad(i)**l*Ln*Lnp*rho
    enddo
    gamma = Ivc*gamma*An*Anp/2._dp
  end function gamma_LDA

!> Evaluates the integral ______the whole thing or just part?_____ in
!! HF_Extensions.pdf eqn. 54
  function IntegralVc() result(Ivc)
    real(dp) :: Ivc
    real(dp) :: alpha,xi,wi,Ir,Is
    integer, parameter :: Ngauss = 77
    real(dp), dimension(1:Ngauss) :: w,x
    logical :: first_call = .true.
    integer :: i
    save w,x,first_call
    if(first_call) then
       alpha = -0.5_dp
       call GaussLaguerreWX(alpha,w,x)
       first_call = .false.
    endif
    Ivc = 0
    Ir = 0
    Is = 0
    do i = 1,Ngauss
       Ir = Ir + w(i)*SphericalBesselJ1(k_Fermi*sqrt(x(i)/muR))**2
       Is = Is + w(i)*SphericalBesselJ1(k_Fermi*sqrt(x(i)/mus))**2
    enddo
    Ir = Ir*V0R/sqrt(muR)
    Is = Is*V0S/sqrt(mus)
    IVc = (Ir-Is)*pi*9/(2*k_Fermi**2)
    Ivc = Ivc + pi**(1.5_dp)*(V0R/(muR**1.5_dp)-V0s/(mus**1.5_dp))/4._dp
    
  end function IntegralVc

!> Spherical Bessel function \f$J_1(x)\f$
  function SphericalBesselJ1(x) result(j1)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: j1
    j1 = sin(x)/x**2 - cos(x)/x
  end function SphericalBesselJ1

!> Writes \f$r\f$ and \f$\rho_{LDA}(r)\f$ (or rather,
!! \f$r^2\rho_{LDA}(r)\f$) to a file \f$\texttt{mixed\_rho\_plot**.dat}\f$
!! for plotting.
  subroutine plot_rho_LDA(j)
    implicit none
    integer, intent(in) :: j
    real(dp) :: Trho
    real(dp), parameter :: dr=0.001_dp
    real(dp) :: ri
    integer :: i
    character(2) :: index
    Trho = 0
    write(index,'(i2.2)') j
    open(100,file='mixed_rho_plot'//index//'.dat')
    do i = 1,10000
       ri = i*dr
       Trho = Trho  + rho_LDA(ri)*ri**2*dr
       write(100,*) ri, rho_LDA(ri)*ri**2
    enddo
    close(100)
  end subroutine Plot_rho_LDA

!> Computes the trace of the density matrix in LDA by integrating
!! the density over coordinate space. If everything is working
!! properly, this integral should return the number of particles.
  function Trace_rho_LDA() result(Trho)
    implicit none
    real(dp) :: Trho
    real(dp), parameter :: dr=0.001_dp
    real(dp) :: ri
    integer :: i
    character(2) :: index
    Trho = 0
    do i = 1,10000
       ri = i*dr
       Trho = Trho  + rho_LDA(ri)*ri**2*dr
    enddo
    Trho = Trho/pi
  end function Trace_rho_LDA

!> I'm not sure what this subroutine was used for.
  subroutine sample_rho_LDA
    implicit none
    real(dp) :: alpha
    logical :: first_call = .true.
    integer :: i
    save first_call
    if(first_call) then
       alpha = 0.5_dp
       call GaussLaguerreWX(alpha,w_quad,x_quad)
       first_call = .false.
    endif
    do i = 1,N_quad
       rho_quad(i) = rho_LDA(b_ho*x_quad(i)**0.5_dp)
    enddo
  end subroutine sample_rho_LDA

!> I don't know what this one does, either.
!! I commented it out for now since I don't have a function
!! DME_fields to call, giving me a compilation error.
!  subroutine sample_DME_fields
!    implicit none
!    real(dp) :: alpha
!    logical :: first_call = .true.
!    integer :: i
!    real :: rho, tau, del_rho
!    save first_call
!    if(first_call) then
!       alpha = 0.5_dp
!       call GaussLaguerreWX(alpha,w_quad,x_quad)
!       first_call = .false.
!    endif
!    do i = 1,N_quad
!       call DME_fields(b_ho*x_quad(i)**0.5_dp,rho,tau,del_rho)
!       rho_quad(i) = rho_LDA(b_ho*x_quad(i)**0.5_dp)
!    enddo
!  end subroutine sample_DME_fields

end module LDA
