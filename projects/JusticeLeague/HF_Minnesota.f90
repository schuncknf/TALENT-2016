module Minnesota
  use :: types
  use :: variables
  use :: HO_Quadrature
  implicit none
  real(dp), private, parameter :: V0R=200.00_dp, muR=1.487_dp
  real(dp), private, parameter :: V0s= 91.85_dp, mus=0.465_dp
  real(dp), private, parameter :: hw=10._dp, hbar=197.327_dp, M_n=939.57_dp
  real(dp), private, parameter :: h2_2m = 20.736209412_dp
  real(dp), private, parameter :: b_ho=sqrt(2*h2_2m/hw)
contains

  subroutine read_orbitals
    implicit none
    integer :: i,j
    character(20) :: text
    integer :: ni,li,ji,mi,tzi
    allocate(HO_inverse(1:n_orbitals))
    allocate(ho_flag(1:n_orbitals))
    allocate(n_ho(1:n_orbitals))
    allocate(m_ho(1:n_orbitals))
    allocate(l_ho(1:n_orbitals))
    allocate(j_ho(1:n_orbitals))
    allocate(tz_ho(1:n_orbitals))
    allocate(n_hf(1:n_orbitals))
    allocate(l_hf(1:n_orbitals))
    allocate(j_hf(1:n_orbitals))
    open(100,file='spM_n5l2.dat')
    do i = 1,1
       read(100,*)
    enddo
    n_hf = 0
    l_hf = 0
    j_hf = 0
    HO_inverse = 0
    ho_flag = 0
    j = 0
    do i = 1,n_orbitals
       read(100,'(a20,6i4)') text, ni,li,ji,mi,tzi
       n_ho(i) = ni
       l_ho(i) = li
       j_ho(i) = ji
       m_ho(i) = mi
       tz_ho(i) = tzi
       if(tzi.eq.1) then
          ho_flag(i) = 1
          if(ni.ne.n_ho(i-1).or.li.ne.l_ho(i-1).or.ji.ne.j_ho(i-1)&
               .or.tzi.ne.tz_ho(i-1)) then
             j = j + 1
             n_hf(j) = ni
             l_hf(j) = li
             j_hf(j) = ji
           endif
           HO_inverse(i) = j
      endif
    enddo
    close(100)
    nsize = 30
  end subroutine read_orbitals

  subroutine read_TBME
    implicit none
    integer :: i1,i2,i3,i4,tz,P,J
    integer :: j1,j2,j3,j4
    integer :: k1,k2,k3,k4
    real(dp) :: TBME
    allocate(v_mat(1:Nsize,1:Nsize,1:Nsize,1:Nsize))
    v_mat = 0
    open(100,file='VM-scheme_n5l2.dat')
    read(100,*)
    read(100,*)
    do
       read(100,*) i1,i2,i3,i4,TBME
       if(i1.eq.0) exit 
       if(ho_flag(i1).eq.1.and.ho_flag(i2).eq.1.and.&
            ho_flag(i3).eq.1.and.ho_flag(i4).eq.1) then
          if(m_ho(i1).ne.m_ho(i3).or.m_ho(i2).ne.m_ho(i4)) cycle
          j1 = HO_inverse(i1)
          j2 = HO_inverse(i2)
          j3 = HO_inverse(i3)
          j4 = HO_inverse(i4)
          v_mat(j1,j2,j3,j4) = v_mat(j1,j2,j3,j4) + &
               1/((j_ho(i1)+1._dp)*(j_ho(i2)+1._dp))*TBME
       endif
    enddo
    close(100)
  end subroutine read_TBME

  subroutine Initialize_Minnseota
    implicit none
    integer :: i,n,np,l
    real(dp) :: alpha,xi,xj,wi,wj,A1,A2,A3,A4,s,L1,L2
    integer, parameter :: Ngauss = 95
    real(dp), dimension(1:Ngauss) :: w,x
    allocate(t_mat(1:Nsize,1:Nsize))
    t_mat = 0
    do i = 1,nsize
       t_mat(i,i) = hw*(2*n_hf(i) + l_hf(i) + 1.5_dp)
    enddo
  end subroutine Initialize_Minnseota

  function fermi_level() result(Noccupied)
    integer :: i,j,k,Noccupied
    k=1
    do i = 1,n_orbitals
        j = HO_inverse(i)

        if (j.ne.0.and.k.le.Nparticles) then
         Noccupied = j
         k=k+1
     elseif (k.gt.Nparticles) then
         write (*,*) "nOccupied = ", Noccupied
         exit
     endif
    end do
  end function fermi_level

  subroutine calculate_TBME
    implicit none
    integer :: n1,n2,n3,n4,i1,i2,i3,i4
    real(dp) :: M34, M43
    do i1 = 1,Nsize
       n1 = n_ho(i1)
       do i2 = i1,Nsize
          n2 = n_ho(i2)
          do i3 = 1,Nsize
             n3 = n_ho(i3)
             do i4 = i3,Nsize
                n4 = n_ho(i4)
                M34 = Minnesota_TBME(n1,n2,n3,n4)
                if(n3.eq.n4) then
                   M43 = M34
                else
                   M43 = Minnesota_TBME(n1,n2,n4,n3)
                endif
                v_mat(i1,i2,i3,i4) = M34+M43
                v_mat(i1,i2,i4,i3) = M34+M43
                v_mat(i2,i1,i3,i4) = M34+M43
                v_mat(i2,i1,i4,i3) = M34+M43
             enddo
          enddo
       enddo
    enddo
  end subroutine calculate_TBME

  function Minnesota_TBME(n1,n2,n3,n4) result(TBME)
    implicit none
    integer, intent(in) :: n1,n2,n3,n4
    integer :: i,j
    integer, parameter :: Ngauss = 95
    real(dp) :: TBME
    real(dp) :: alpha,xi,xj,wi,wj,A1,A2,A3,A4
    real(dp), dimension(1:Ngauss) :: w,x
    logical :: first_call = .true.
    save w,x
    if(first_call) then
       alpha = 0.5_dp
       call GaussLaguerreWX(alpha,w,x)
       first_call = .false.
    endif
    A1 = HO_Normalization(n1,0)
    A2 = HO_Normalization(n2,0)
    A3 = HO_Normalization(n3,0)
    A4 = HO_Normalization(n4,0)
    TBME = 0._dp
    do i = 1,Ngauss
       xi = x(i)
       wi = w(i)
       do j = 1, Ngauss
          xj = x(j)
          wj = w(j)
          TBME = TBME + wi*wj*Minnesota_Kernel(n1,n2,n3,n4,xi,xj)
       enddo
    enddo

    TBME = TBME*A1*A2*A3*A4/(4._dp)
  end function Minnesota_TBME

  function Minnesota_Kernel(n1,n2,n3,n4,xi,xj) result(MK)
    implicit none
    integer, intent(in) :: n1,n2,n3,n4
    real(dp), intent(in) :: xi,xj
    real(dp) :: MK
    real(dp) :: L1,L2,L3,L4,ri,rj
    call LaguerreL(n1,0.5_dp,xi,L1)
    call LaguerreL(n2,0.5_dp,xj,L2)
    call LaguerreL(n3,0.5_dp,xi,L3)
    call LaguerreL(n4,0.5_dp,xj,L4)
    ri = b_ho*sqrt(xi)
    rj = b_ho*sqrt(xj)
    MK = (V0R/mur*(Exp(-mur*(ri**2+rj**2-2*ri*rj))&
                  -Exp(-mur*(ri**2+rj**2+2*ri*rj))) &
         -V0s/mus*(Exp(-mus*(ri**2+rj**2-2*ri*rj))&
                  -Exp(-mus*(ri**2+rj**2+2*ri*rj))))&
                  *L1*L2*L3*L4/(4*ri*rj)*0.5_dp
  end function Minnesota_Kernel

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
    do i = 1,3
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
             phi = phi + (D_mat( j,i)*D_mat( k,i)*1.0_dp + &
                          D_prev(j,i)*D_prev(k,i)*0.0_dp)*&
                          xi**(lj+lk)*exp(-xi**2)*&
                          L_nl(nj,lj)*L_nl(nk,lk)*Anl*Anpl/(b_ho**3)
          enddo
       enddo
       rho = rho + (ji+1)*phi
    enddo
    rho = rho/(4*pi)
  end function rho_LDA

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
!          if(i.eq.1.and.i.eq.j) write(*,*) 'Gamma 001', gamma
       enddo
    enddo
    
  end subroutine calculate_gamma_LDA

  function gamma_LDA(n,np,l) result(gamma)
    implicit none
    integer, intent(in) :: n,np,l
    real(dp) :: gamma
    integer, parameter :: Ngauss = 95
    real(dp) :: alpha,xi,wi,An,Anp,Ln,Lnp,rho,Ivc
    real(dp), dimension(1:Ngauss) :: w,x
    logical :: first_call = .true.
    integer :: i
    save w,x,first_call,Ivc
    if(first_call) then
       Ivc = IntegralVc()
       alpha = 0.5_dp
       call GaussLaguerreWX(alpha,w,x)
       first_call = .false.
    endif
    An  = HO_Normalization(n ,l)
    Anp = HO_Normalization(np,l)
    gamma = 0
    do i = 1,Ngauss
       call LaguerreL( n,l+0.5_dp,x(i),Ln)
       call LaguerreL(np,l+0.5_dp,x(i),Lnp)
       rho = rho_LDA(b_ho*x(i)**0.5_dp)
       gamma = gamma + w(i)*x(i)**l*Ln*Lnp*rho
    enddo
    gamma = Ivc*gamma*An*Anp/2._dp
  end function gamma_LDA

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

  function SphericalBesselJ1(x) result(j1)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: j1
    j1 = sin(x)/x**2 - cos(x)/x
  end function SphericalBesselJ1

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
!    stop
  end subroutine Plot_rho_LDA

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

end module Minnesota
