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
!! Approximation (see eqn. 60 of HF_extensions.pdf).
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
             phi = phi + D_mat( j,i)*D_mat( k,i)* &
!                          D_prev(j,i)*D_prev(k,i)*0.0_dp)*&
                          xi**(lj+lk)*exp(-xi**2)*&
                          L_nl(nj,lj)*L_nl(nk,lk)*Anl*Anpl/(b_ho**3)
          enddo
       enddo
       rho = rho + (ji+1)*phi
    enddo
    rho = rho/(4*pi)
  end function rho_LDA

!> Calls the function \f$\texttt{gamma\_LDA(n,n',l)}\f$ to calculate
!! the matrix elements \f$\Gamma_{ij}, \Gamma_{ji}\f$ (see
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
          if(type_of_calculation.eq.'DME') then
             gamma = gamma_DME(ni,nj,li)
          else
             gamma = gamma_LDA(ni,nj,li)
          endif
          gamma_mat(i,j) = gamma
          gamma_mat(j,i) = gamma
       enddo
    enddo
    
  end subroutine calculate_gamma_LDA

!> Performs the integral in HF_Extensions.pdf, eqn. 58.
!! \f$\texttt{(n, np, l)}=n,n',l\f$
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

!> Evaluates the integral in HF_Extensions.pdf eqn. 54
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

  function SphericalBesselJ2(x) result(j2)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: j2
    j2 = (-3*Cos(x))/x**2 + ((3 - x**2)*Sin(x))/x**3
  end function SphericalBesselJ2

!> Spherical Bessel function \f$J_3(x)\f$
  function SphericalBesselJ3(x) result(j3)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: j3
    j3 = (15/x**3-6/x)*sin(x)/x - (15/x**2-1)*cos(x)/x
  end function SphericalBesselJ3

  function SphericalBesselJ4(x) result(j4)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: j4
    j4 = (5*(-21 + 2*x**2)*Cos(x))/x**4 + &
         ((105 - 45*x**2 + x**4)*Sin(x))/x**5
  end function SphericalBesselJ4


!> Writes \f$r\f$ and \f$\rho_{LDA}(r)\f$ (or rather,
!! \f$r^2\rho_{LDA}(r)\f$) to a file \f$\texttt{mixed\_rho\_plot**.dat}\f$
!! for plotting.
  subroutine plot_rho_LDA
    implicit none
!    integer, intent(in) :: j
    real(dp) :: Trho
    real(dp), parameter :: dr=0.001_dp
    real(dp) :: ri
    integer :: i
    character(2) :: index
!    Trho = 0
!    write(index,'(i2.2)') j
    open(100,file=trim(density_file))
    do i = 1,10000
       ri = i*dr
!       Trho = Trho  + rho_LDA(ri)*ri**2*dr
       write(100,*) ri, rho_LDA(ri)!*ri**2
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
    Trho = Trho*4*pi
  end function Trace_rho_LDA

!> Calls upon \f$\texttt{rho\_LDA}\f$ at an array of points defined
!! by the Gauss-
!! Laguerre quadrature mesh, and stores the result in a global array
!! \f$\texttt{rho\_LDA}\f$ to be used by \f$\texttt{gamma\_LDA}\f$
!! to calculate integral 58 of HF_extensions.pdf.
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

!> Calculates, at an arbitrary point \f$R\f$, the fields \f$\rho(R)\f$,
!! \f$\tau(R)=\nabla_1\cdot\nabla_2\rho(r_1r_2)|_{r_1=r_2=R}\f$
!! and \f$\nabla\rho(r)\f$.
!! The subroutine \f$\texttt{sample\_DME\_fields}\f$ will call this
!! subroutine for each point in the integration/quadrature mesh.
  subroutine DME_fields(r,rho,tau,del_rho)
    implicit none
    real(dp), intent(in) :: r
    real(dp), intent(out) :: rho,tau,del_rho
    real(dp), dimension(0:2,0:5,0:2) :: Rnl
    real(dp) :: xi
    integer :: i,ji,j,nj,lj,jj,k,nk,lk,jk
    xi = r/b_ho
    call RadialHOALL(5,2,xi,b_ho,Rnl)
    rho = 0
    tau = 0
    del_rho = 0
    do i = 1,Noccupied
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
             rho = rho + (ji+1)*D_mat( j,i)*D_mat( k,i)*&
                  Rnl(0,nj,lj)*Rnl(0,nk,lk)
             tau = tau + (ji+1)*D_mat( j,i)*D_mat( k,i)*&
                  (Rnl(1,nj,lj)*Rnl(1,nk,lk)+&
                  lj*(lj+1)/r**2*Rnl(0,nj,lj)*Rnl(0,nk,lk))
             del_rho = del_rho + (ji+1)*D_mat( j,i)*D_mat( k,i)*&
                  (2*(Rnl(0,nj,lj)*Rnl(1,nk,lk)/r + &
                      Rnl(0,nk,lk)*Rnl(1,nj,lj)/r + &
                      Rnl(1,nj,lj)*Rnl(1,nk,lk)) + & 
                   Rnl(0,nj,lj)*Rnl(2,nk,lk) + Rnl(0,nk,lk)*Rnl(2,nj,lj))

          enddo
       enddo
    enddo
    rho = rho/(4*pi)
    tau = tau/(4*pi)
    del_rho = del_rho/(4*pi)
  end subroutine DME_fields


!> Performs the integration (using Guass-Laguerre quadrature with weight
!! function w) in HF_extensions.pdf eqn. 58, except using
!! DME instead of LDA. \f$\texttt{(n, np, l)}=n,n',l\f$
  function gamma_DME(n,np,l) result(gamma)
    implicit none
    integer, intent(in) :: n,np,l
    real(dp) :: gamma
    real(dp) :: alpha,xi,wi,An,Anp,Ln,Lnp,rho,tau,delrho
    real(dp) :: crho,ctau,cdelrho,dcrho,dctau,dcdelrho
    integer :: i
!    if(calc_couplings) then
!       call calculte_couplings
!       calc_couplings = .false.
!    endif
    An  = HO_Normalization(n ,l)
    Anp = HO_Normalization(np,l)
    gamma = 0
    do i = 1,N_quad
       call LaguerreL( n,l+0.5_dp,x_quad(i),Ln)
       call LaguerreL(np,l+0.5_dp,x_quad(i),Lnp)
       rho = rho_quad(i)
       tau = tau_quad(i)
       delrho = delrho_quad(i)
       crho = c_rhorho_quad(i)
       ctau = c_rhotau_quad(i)
       cdelrho = c_rhodelrho_quad(i)
       dcrho = dc_rhorho_quad(i)
       dctau = dc_rhotau_quad(i)
       dcdelrho = dc_rhodelrho_quad(i)
       gamma = gamma + w_quad(i)*x_quad(i)**l*Ln*Lnp &
!            *(-18.25_dp*rho + 4.57_dp*tau - 1.8_dp*delrho)
            *(dcrho*rho**2 + 2*Crho*rho + dctau*rho*tau + Ctau*tau + &
            dcdelrho*rho*delrho + Cdelrho*delrho)
    enddo
    gamma = gamma*An*Anp/2._dp
  end function gamma_DME

!> Calls upon \f$\texttt{DME\_fields}\f$ at an array of points defined
!! by the Gauss-
!! Laguerre quadrature mesh, and stores the result in the global arrays
!! \f$\texttt{rho\_quad}\f$, \f$\texttt{tau\_quad}\f$, and
!! \f$\texttt{delrho\_quad}\f$ to be used by \f$\texttt{gamma\_DME}\f$
!! to calculate integral 58 of HF_extensions.pdf.
  subroutine sample_DME_fields
    implicit none
    real(dp) :: alpha
    logical :: first_call = .true.
    integer :: i
    real(dp) :: rho,tau,del_rho
    real(dp) :: xi
    real(dp), dimension(0:3,0:5,0:5) :: Rnl
    save first_call
   if(first_call) then
      alpha = 0.5_dp
      call GaussLaguerreWX(alpha,w_quad,x_quad)
      first_call = .false.
   endif
    do i = 1,N_quad
       call DME_fields(b_ho*x_quad(i)**0.5_dp,rho,tau,del_rho)
       rho_quad(i) = rho
       tau_quad(i) = tau
       delrho_quad(i) = del_rho
    enddo
  end subroutine sample_DME_fields

  subroutine sample_DME_couplings
    implicit none
    integer :: i
    real(dp) :: crho,ctau,cdelrho,dcrho,dctau,dcdelrho
    do i = 1,N_quad
       call calculate_couplings(i,crho,ctau,cdelrho,dcrho,dctau,dcdelrho)
       C_rhorho_quad(i) = crho
       C_rhotau_quad(i) = ctau
       C_rhodelrho_quad(i) = cdelrho
       dC_rhorho_quad(i) = dcrho
       dC_rhotau_quad(i) = dctau
       dC_rhodelrho_quad(i) = dcdelrho
    enddo
  end subroutine sample_DME_couplings


!> Here the coupling constants \f$C^{\rho\rho}, C^{\rho\tau}\f$, and
!! \f$C^{\rho\nabla^2\rho}\f$ are calculated using the kernels defined in
!! the functions \f$\texttt{C\_rhorho\_kernel(r)}\f$ and
!! \f$\texttt{C\_rhotau\_kernel(r)}\f$ (the kernel of the
!! \f$C^{\rho\nabla^2\rho}\f$ term is equal to
!! \f$-\frac{C^{\rho\tau}}{4}\f$; see  HF_extensions.pdf eqn. 63).
  subroutine calculate_couplings(i_quad,crho,ctau,cdelrho,dcrho,dctau,dcdelrho)
    implicit none
    integer, intent(in) :: i_quad
    real(dp), intent(out) :: crho,ctau,cdelrho,dcrho,dctau,dcdelrho
    real(dp) :: alpha,xi,wi,KF
    integer, parameter :: Ngauss = 77
    real(dp), dimension(1:Ngauss) :: w,x
    logical :: first_call = .true.
    integer :: i
    real(dp) :: I_R, I_S, J_R, J_S, C_Hdelrho, dI_R, dI_S, dJ_R, dJ_S
    save w,x,first_call
    if(first_call) then
       alpha = -0.5_dp
       call GaussLaguerreWX(alpha,w,x)
       first_call = .false.
    endif
    C_hartree =  pi**(1.5_dp)*(V0R/(muR**1.5_dp)-V0s/(mus**1.5_dp))/8._dp
    C_Hdelrho =3*pi**(1.5_dp)*(V0R/(muR**2.5_dp)-V0s/(mus**2.5_dp))/8._dp
    I_R = 0
    I_S = 0
    J_R = 0
    J_S = 0
    dI_R = 0
    dI_S = 0
    dJ_R = 0
    dJ_S = 0
    do i = 1,N_quad
       I_R = I_R + w_quad(i)*C_rhorho_kernel((x_quad(i)/muR)**0.5_dp,i_quad)
       I_S = I_S + w_quad(i)*C_rhorho_kernel((x_quad(i)/muS)**0.5_dp,i_quad)
       J_R = J_R + w_quad(i)*C_rhotau_kernel((x_quad(i)/muR)**0.5_dp,i_quad)
       J_S = J_S + w_quad(i)*C_rhotau_kernel((x_quad(i)/muS)**0.5_dp,i_quad)
       dI_R=dI_R + w_quad(i)*DC_rhorho_kernel((x_quad(i)/muR)**0.5_dp,i_quad)
       dI_S=dI_S + w_quad(i)*DC_rhorho_kernel((x_quad(i)/muS)**0.5_dp,i_quad)
       dJ_R=dJ_R + w_quad(i)*DC_rhotau_kernel((x_quad(i)/muR)**0.5_dp,i_quad)
       dJ_S=dJ_S + w_quad(i)*DC_rhotau_kernel((x_quad(i)/muS)**0.5_dp,i_quad)
    enddo
    Crho=pi*(V0R*I_R/muR**1.5_dp-V0S*I_S/muS**1.5_dp)*0.5_dp
    Ctau=pi*(V0R*J_R/muR**1.5_dp-V0S*J_S/muS**1.5_dp)*0.5_dp
    Cdelrho = (-Ctau/4._dp+C_Hdelrho/8._dp)
    Crho = Crho + C_hartree
    dCrho=pi*(V0R*dI_R/muR**1.5_dp-V0S*dI_S/muS**1.5_dp)*0.5_dp
    dCtau=pi*(V0R*dJ_R/muR**1.5_dp-V0S*dJ_S/muS**1.5_dp)*0.5_dp
    dCdelrho = -dCtau/4._dp
  end subroutine calculate_couplings

  function pi0(x) result(p0)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: p0
    if(x.lt.0.01_dp) then
       p0 = 1-x**2/10._dp + x**4/280._dp
    else
       p0 = 3*SphericalBesselJ1(x)/x
    endif
  end function pi0

  function pi2(x) result(p2)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: p2
    if(x.lt.0.01_dp) then
       p2 = 1-x**2/18._dp + x**4/792._dp
    else
       p2 = 105*SphericalBesselJ3(x)/x**3
    endif
  end function pi2

  function dpi0(x) result(p0)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: p0
    if(x.lt.0.01_dp) then
       p0 = -x/5._dp + x**3/70._dp
    else
       p0 = 3*(sin(x)-3*SphericalBesselJ1(x))/x**2
    endif
  end function dpi0

  function dpi2(x) result(p2)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: p2
    if(x.lt.0.01_dp) then
       p2 = -x/9._dp + x**3/198._dp
    else
       p2 = 105*(x*SphericalBesselJ2(x)-7*SphericalBesselJ3(x))/x**4
    endif
  end function dpi2
  

  function DC_rhorho_kernel(r,i) result(DCk)
    implicit none
    integer, intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: DCk, Kf,dKf,p0,p2,dp0,dp2
    Kf = (3*pi**2*rho_quad(i))**(1/3._dp)
    dKF = (3*pi**2)**(1/3._dp)*rho_quad(i)**(-2/3._dp)/3._dp
    p0 = pi0(kf*r)
    p2 = pi2(kf*r)
    dp0 = dpi0(kf*r)
    dp2 = dpi2(kf*r)
    DCk = r*dKF*(2*r*Kf*p0*p2 + 10*p0*dp0 + &
         (r*kf)**2*(p2*dp0+p0*dp2))*0.1_dp
  end function DC_rhorho_kernel

  function DC_rhotau_kernel(r,i) result(DCk)
    implicit none
    integer, intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: DCk, Kf,dKf,p0,p2,dp0,dp2
    Kf = (3*pi**2*rho_quad(i))**(1/3._dp)
    dKF = (3*pi**2)**(1/3._dp)*rho_quad(i)**(-2/3._dp)/3._dp
    p0 = pi0(kf*r)
    p2 = pi2(kf*r)
    dp0 = dpi0(kf*r)
    dp2 = dpi2(kf*r)
    DCk = -r**3*dkf*(p2*dp0+p0*dp2)/6._dp
  end function DC_rhotau_kernel


!> Contains the kernel of the integral used to compute \f$C^{\rho\rho}\f$.
!! Essentially, it is everything inside the integral in eqn. 64 of
!! HF_extensions.pdf, except the integral has been transformed into a
!! form that permits it to be evaluated using our Gauss-Laguerre
!! quadrature scheme.
  function C_rhorho_kernel(r,i) result(Ck)
    implicit none
    integer, intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: Ck
    real(dp) :: Kf
    Kf = (3*pi**2*rho_quad(i))**(1/3._dp)
    Ck = (pi0(kf*r)**2+(kf*r)**2*pi0(kf*r)*pi2(kf*r)/5._dp)*0.5_dp
!    Ck =  9*SphericalBesselJ1(Kf*r)**2 + &
!         63*SphericalBesselJ1(Kf*r)*SphericalBesselJ3(kf*r)
  end function C_rhorho_kernel

!> Contains the kernel of the integral used to compute \f$C^{\rho\tau}\f$
!! (see  eqn. 63 of HF_extensions.pdf) As in the case of
!! \f$C^{\rho\rho}\f$, the integral has been transformed into a
!! form that permits it to be evaluated using our Gauss-Laguerre
!! quadrature scheme.
  function C_rhotau_kernel(r,i) result(Ck)
    implicit none
    integer, intent(in) :: i
    real(dp), intent(in) :: r
    real(dp) :: Ck
    real(dp) :: Kf
    Kf = (3*pi**2*rho_quad(i))**(1/3._dp)
    Ck = -pi0(kf*r)*pi2(kf*r)*r**2/6._dp
  end function C_rhotau_kernel


!> Writes \f$r\f$ and \f$\rho_{DME}(r)\f$ (or rather,
!! \f$r^2\rho_{DME}(r)\f$) to a file
!! \f$\texttt{dmefields\_surface\_sat.dat}\f$ for plotting.
  subroutine plot_DME_fields
    implicit none
    integer :: i
    real(dp) :: rho,tau,del_rho
    real(dp) :: xi
    real(dp) :: Tr, dr
    dr = 0.001
    Tr = 0._dp
    open(100,file='dmefields_surface_sat.dat')
    do i = 1,10000
       xi = i*dr
       call DME_fields(xi,rho,tau,del_rho)
       Tr = Tr + rho*xi**2*dr
       write(100,*) xi,rho,tau,del_rho
    enddo
    close(100)
    Tr = Tr*4*pi
    write(*,*) 'Trace of DME rho', Tr
  end subroutine plot_DME_fields

!> Computes the RMS radius of the neutron drop in the DME approximation.
  function rms_DME() result(rms)
    implicit none
    real(dp) :: rms
    integer :: i
    real(dp) :: dr, xi, rho,tau,del_rho
    dr = 0.01
    rms = 0._dp
    do i = 1,1000
       xi = i*dr
       call DME_fields(xi,rho,tau,del_rho)
       rms = rms + rho*xi**4*dr
    enddo
    rms = sqrt(rms*4*pi/real(Nparticles,kind=dp))
  end function rms_DME

end module LDA
