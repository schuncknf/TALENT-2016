module Minnesota
  use :: types
  use :: variables
  use :: HO_Quadrature
  implicit none
  real(dp), private, parameter :: V0R=200.00_dp, muR=1.487_dp
  real(dp), private, parameter :: V0s= 91.85_dp, mus=0.465_dp
  real(dp), private, parameter :: hw=10._dp, hbar=197.327_dp, M_n=939.57_dp
  real(dp), private, parameter :: b_ho=hbar/sqrt(M_n*hw)
contains

  subroutine Initialize_Minnseota
    implicit none
    integer :: i
    if(allocated(t_mat)) then
       deallocate(t_mat,v_mat,n_ho,l_ho)
    endif
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
    call calculate_TBME

  end subroutine Initialize_Minnseota

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
                M43 = Minnesota_TBME(n1,n2,n4,n3)
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
    integer, parameter :: Ngauss = 90
    real(dp) :: TBME
    real(dp) :: alpha,xi,xj,wi,wj,A1,A2,A3,A4
    real(dp), dimension(1:Ngauss) :: w,x
    alpha = 0.5_dp
    call GaussLaguerreWX(alpha,w,x)
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

    TBME = TBME*A1*A2*A3*A4/4._dp
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
                  -Exp(-mur*(ri**2+rj**2))) &
         -V0s/mus*(Exp(-mus*(ri**2+rj**2-2*ri*rj))&
                  -Exp(-mus*(ri**2+rj**2))))*L1*L2*L3*L4/(4*ri*rj)*0.5_dp
  end function Minnesota_Kernel

end module Minnesota
