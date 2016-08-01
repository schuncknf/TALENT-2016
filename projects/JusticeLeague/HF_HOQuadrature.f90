!---------------------------------------------------------------------
!Module:HO_Quadrature 
!---------------------------------------------------------------------
!> Contains the functions and subroutines needed for defining the
!! Laguerre polynomials, which describe the radial part of the
!! harmonic oscillator wave function (which is our starting basis).
!> Also contains functions and subroutines for evaluating integrals
!! of these functions, using Gauss-Laguerre quadrature.
!---------------------------------------------------------------------
module HO_Quadrature
  use :: types
  implicit none
contains
 
  subroutine GaussLaguerreWX(alfa,w,x)
    implicit none
    real(dp), intent(in) :: alfa
    real(dp), intent(out), dimension(:) :: x
    real(dp), intent(out), dimension(:) :: w
    integer, parameter :: maxit = 100
    integer i,its,n,ai
    real(dp), parameter :: eps = 2.e-14_dp
    real(dp) :: z,z1,Ln,Lnp,Lnm1
    n = size(w)
    if(size(w).ne.size(x)) then !Check that both arrays have the same size
       write(*,*) 'Arrays w and x must have the same size in GaussLaguerreWX'
       stop
    end if
    do i = 1,n
       !intial guess for every root of the Laguerre Polynomial
       if(i.eq.1)then
          z=(1+alfa)*(3+0.92_dp*alfa)/(1+2.4_dp*n+1.8_dp*alfa)
       else if(i.eq.2)then
          z=z+(15+6.25_DP*alfa)/(1+0.9_dp*alfa+2.5_dp*n)
       else
          ai=i-2
          z=z+((1+2.55_dp*ai)/(1.9_dp*ai)+1.26_dp*ai*alfa/(1+3.5_dp*ai))&
               *(z-x(i-2))/(1+0.3_dp*alfa)
       endif
       ! Newton method to refine the roots
       do its = 1, maxit
          call LaguerreL(n,alfa,z,Ln,Lnp,Lnm1)
          z1 = z
          z = z1 - Ln/Lnp
!          if(abs(z-z1).le.eps) exit
          if(abs(z-z1).le.eps*z) exit
       enddo
!       if(its==maxit+1) write(*,*) 'maxit exceeded in GaussLaguerreWX',i,abs(z-z1),z
       if(its==maxit+1) write(*,*) 'maxit exceeded in GaussLaguerreWX',its,abs(z-z1)/z,z
       ! Save root and calculate weight
       x(i) = z
       if(alfa.eq.0._dp) then
          w(i) = -1/(n*Lnm1*Lnp)
       else
          w(i) = -gamma(n+alfa)/(factrl(n)*Lnm1*Lnp)
       endif
    end do
  end subroutine GaussLaguerreWX


  subroutine LaguerreL(n,alpha,x,Ln,Lnp,Lnm1)
    implicit none
    integer, intent(in) :: n    
    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: x   
    real(dp), intent(out) :: Ln 
    real(dp), optional, intent(out) :: Lnp  
    real(dp), optional, intent(out) :: Lnm1 
    real(dp) :: Ljm2, Ljm1, Lj
    integer :: j
    Lj = 1._dp
    Ljm1 = 0._dp
    do j = 1,n
       Ljm2 = Ljm1
       Ljm1 = Lj
       Lj = ((-x+2*j-1+alpha)*Ljm1-(j-1+alpha)*Ljm2)/real(j,dp)
    end do
    Ln = Lj
    if(present(Lnp)) Lnp = (n*Ln-(n+alpha)*Ljm1)/x
    if(present(Lnm1)) Lnm1 = Ljm1
    return
  end subroutine LaguerreL

  subroutine LaguerreALL(n_max,l_max,x,Ln)
    implicit none
    integer, intent(in) :: n_max
    integer, intent(in) :: l_max
    real(dp), intent(in) :: x   
    real(dp), intent(out), dimension(0:n_max,0:l_max) :: Ln 
    real(dp) :: alpha
    real(dp) :: Ljm2, Ljm1, Lj
    integer :: j,l
    do l = 0,l_max
       alpha = l + 0.5_dp
       Lj = 1._dp
       Ljm1 = 0._dp
       Ln(0,l) = Lj
       do j = 1,n_max
          Ljm2 = Ljm1
          Ljm1 = Lj
          Lj = ((-x+2*j-1+alpha)*Ljm1-(j-1+alpha)*Ljm2)/real(j,dp)
          Ln(j,l) = Lj
       end do
    enddo
    return
  end subroutine LaguerreALL


!> Numerically evaluates the value of the factorial function,
!! \f$n!=n(n-1)(n-2)...\f$.
  function factrl(n) result(fact)
    implicit none
    integer, intent(in) :: n 
    real(dp) :: fact
    integer, save :: ntop = 0
    integer, parameter :: nmax = 170
    integer :: i
    real(dp), dimension(0:nmax), save :: a=0._dp 
    if(n.lt.0) then
       write(*,*) 'negative integer in factrl'
       stop
    endif
    if(ntop.eq.0) a(0)=1._dp
    if(n.le.ntop) then
       fact = a(n)
    elseif(n.le.nmax) then
       do i = ntop+1,n
          a(i) = real(i,kind=dp)*a(i-1)
       enddo
       ntop = n
       fact = a(n)
    else
       fact = exp(log_gamma(n+1._dp))
    endif
  end function factrl

!> Numerically evaluates the value of the double-factorial function,
!! \f$n!!=n(n-2)(n-4)...\f$.
  function doublefactrl(n) result(fact)
    implicit none
    integer, intent(in) :: n !< An integer
    real(dp) :: fact
    integer, save :: ntop = 1
    integer, parameter :: nmax = 300
    integer :: i
    real(dp), dimension(0:nmax), save :: a=0._dp 
    if(n.lt.0) then
       write(*,*) 'negative integer in doublefactrl'
       stop
    endif
    if(ntop.eq.1) then
       a(0)=1._dp
       a(1)=1._dp
    endif
    if(n.le.ntop) then
       fact = a(n)
    elseif(n.le.nmax) then
       do i = ntop+1,n
          a(i) = real(i,kind=dp)*a(i-2)
       enddo
       ntop = n
       fact = a(n)
    else
       fact = exp(log_gamma(1+n*0.5_dp)+(n+1)*0.5_dp*log(2._dp))/sqrt(pi)
    endif
  end function doublefactrl

!> Evaluates the normalization prefactor for the 3D harmonic oscillator,
!! as written in eqn. 5 and 6 of HO_basis.pdf, for a given value of
!! \f$n\f$ and \f$l\f$.
  function HO_Normalization(n,l) result(Nnl)
    implicit none
    integer, intent(in) :: n,l
    real(dp) :: Nnl
    Nnl = sqrt((2**(n+l+2)*factrl(n))/(sqrt(pi)*doublefactrl(2*n+2*l+1)))
  end function HO_Normalization

end module HO_Quadrature
