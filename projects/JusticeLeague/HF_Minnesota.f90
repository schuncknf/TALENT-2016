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
!    character(11) :: text
    character(18) :: text
    integer :: ni,li,ji,mi,tzi
    allocate(HO_inverse(1:n_orbitals))
    allocate(ho_flag(1:n_orbitals))
    allocate(n_ho(1:n_orbitals))
    allocate(m_ho(1:n_orbitals))
    allocate(l_ho(1:n_orbitals))
    allocate(j_ho(1:n_orbitals))
!    open(100,file='spJ.dat')
    open(100,file='spM.dat')
!    do i = 1,10
    do i = 1,1
       read(100,*)
    enddo
    HO_inverse = 0
    ho_flag = 0
    j = 0
    do i = 1,n_orbitals
!       read(100,'(a11,4i6)') text, ni,li,ji,tzi
       read(100,'(a18,6i4)') text, ni,li,ji,mi,tzi
       if(li.eq.0.and.tzi.eq.1) then
          j = j + 1
          n_ho(j) = ni
          l_ho(j) = li
          j_ho(j) = ji
          m_ho(j) = mi
          ho_flag(i) = 1
          HO_inverse(i) = j
       endif
    enddo
    close(100)
!    nsize = j
    nsize = j/2
  end subroutine read_orbitals

  subroutine read_TBME
    implicit none
    integer :: i1,i2,i3,i4,tz,P,J
    integer :: j1,j2,j3,j4
    integer :: k1,k2,k3,k4
    real(dp) :: TBME
    allocate(v_mat(1:Nsize,1:Nsize,1:Nsize,1:Nsize))
    v_mat = 0
!    open(100,file='VJ-scheme.dat')
    open(100,file='VM-scheme.dat')
    read(100,*)
    read(100,*)
    do
!       read(100,*) tz,P,J,i1,i2,i3,i4,TBME
       read(100,*) i1,i2,i3,i4,TBME
       
!       if(tz.eq.9) exit
       if(i1.eq.0) exit 
       if(ho_flag(i1).eq.1.and.ho_flag(i2).eq.1.and.&
            ho_flag(i3).eq.1.and.ho_flag(i4).eq.1) then
          k1 = HO_inverse(i1)
          k2 = HO_inverse(i2)
          k3 = HO_inverse(i3)
          k4 = HO_inverse(i4)
          if(m_ho(k1).ne.m_ho(k3).or.m_ho(k2).ne.m_ho(k4)) cycle
          j1 = (HO_inverse(i1)+1)/2
          j2 = (HO_inverse(i2)+1)/2
          j3 = (HO_inverse(i3)+1)/2
          j4 = (HO_inverse(i4)+1)/2
          ! if(j1.eq.1.and.j2.eq.1.and.j3.eq.1.and.j4.eq.2) &
          !      write(*,*) m_ho(k1),m_ho(k2),m_ho(k3),m_ho(k4),TBME
!          v_mat(j1,j2,j3,j4) = v_mat(j1,j2,j3,j4) + &
!               (J+1)/((j_ho(i1)+1._dp)*(j_ho(i2)+1._dp))*TBME
!          if(j1.eq.4.and.j2.eq.4.and.j3.eq.4.and.j4.eq.4) &
!               write(*,*) TBME
          v_mat(j1,j2,j3,j4) = v_mat(j1,j2,j3,j4) + &
               1/((j_ho(k1)+1._dp)*(j_ho(k2)+1._dp))*TBME
       endif
    enddo
    close(100)
    v_mat = 2*v_mat
    ! do j1 = 1,nsize
    !    do j2 = j1,nsize
    !       do j3 = 1,nsize
    !          do j4 = j3,nsize
    !             v_mat(j1,j2,j4,j3) = v_mat(j1,j2,j3,j4)
    !             v_mat(j2,j1,j3,j4) = v_mat(j1,j2,j3,j4)
    !             v_mat(j2,j1,j4,j3) = v_mat(j1,j2,j3,j4)
    !          enddo
    !       enddo
    !    enddo
    ! enddo
   !  do j1 = 1,nsize
   !     do j2 = 1,nsize
   !        do j3 = 1,nsize
   !           do j4 = 1,nsize
   !              write(*,'(4i4,1f15.8)') j1,j2,j3,j4,v_mat(j1,j2,j3,j4)
   !           enddo
   !        enddo
   !     enddo
   !  enddo
   ! stop
  end subroutine read_TBME

  subroutine Initialize_Minnseota
    implicit none
    integer :: i,n,np,l
    real(dp) :: alpha,xi,xj,wi,wj,A1,A2,A3,A4,s,L1,L2
    integer, parameter :: Ngauss = 95
    real(dp), dimension(1:Ngauss) :: w,x
    allocate(t_mat(1:Nsize,1:Nsize))
!    allocate(v_mat(1:Nsize,1:Nsize,1:Nsize,1:Nsize))
!    allocate(n_ho(1:Nsize))
!    allocate(l_ho(1:Nsize))
    t_mat = 0
!    v_mat = 0
    n_ho = 0
    l_ho = 0
    do i = 1,nsize
       n_ho(i) = i-1
       t_mat(i,i) = hw*(2*n_ho(i) + l_ho(i) + 1.5_dp)
!       t_mat(i,i) = hw*(2*n_ho(2*i) + l_ho(2*i) + 1.5_dp)
    enddo
!    call calculate_TBME
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
                if(n3.eq.n4) then
                   M43 = M34
                else
                   M43 = Minnesota_TBME(n1,n2,n4,n3)
                endif
!                write(*,*) i1,i2,i3,i4,M34+M43
                v_mat(i1,i2,i3,i4) = M34+M43
                v_mat(i1,i2,i4,i3) = M34+M43
                v_mat(i2,i1,i3,i4) = M34+M43
                v_mat(i2,i1,i4,i3) = M34+M43
             enddo
          enddo
       enddo
    enddo
    do i1 = 1,Nsize
       do i2 = 1,Nsize
          do i3 = 1,Nsize
             do i4 = 1,Nsize
                write(*,'(4i4,f15.8)') i1,i2,i3,i4,v_mat(i1,i2,i3,i4)
             enddo
          enddo
       enddo
    enddo
    stop
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

end module Minnesota
