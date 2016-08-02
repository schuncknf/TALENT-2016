!---------------------------------------------------------------------
!Module:HO_Quadrature 
!---------------------------------------------------------------------
!> Contains the functions and subroutines needed for defining the
!! Laguerre polynomials, which describe the radial part of the
!! harmonic oscillator wave function (which is our starting basis).
!! Also contains functions and subroutines for evaluating integrals
!! of these and other functions, using Gauss-Laguerre quadrature.
!---------------------------------------------------------------------
module HO_Quadrature
  use :: types
  implicit none
contains

!> Gauss-Laguerre quadrature routine 
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

!> Evaluates the Laguerre polynomial \f$L_n^\alpha(x)\f$ using the
!! recurrence relation (see section 1.2 of HO_basis.pdf). Returns
!! \f$\texttt{Ln}\f$. Optionally, the user may include a fifth and sixth
!! argument, in which case the subroutine will also return the derivative
!! \f$L_n'\rightarrow\texttt{Lnp}\f$ and the previous Laguerre polynomial
!! \f$L_{n-1}\rightarrow\texttt{Lnm1}\f$.
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

!> Similar to \f$\texttt{LaguerreAll}\f$, this subroutine computes and
!! stores an array \f$\texttt{Rnl}\f$ in which \f$\texttt{Rnl(0,n,l)}\f$
!! is the radial part of the 3D harmonic oscillator wave function
!! \f$R_{nl}\f$, \f$\texttt{Rnl(1,n,l)}\f$ is its first derivative, and
!! \f$\texttt{Rnl(0,n,l)}\f$ is its second derivative. See eqn. 5 of
!! HO_basis.pdf.
  subroutine RadialHOALL(n_max,l_max,xi,b,Rnl)
    implicit none
    integer, intent(in) :: n_max
    integer, intent(in) :: l_max
    real(dp), intent(in) :: xi
    real(dp), intent(in) :: b
    real(dp), intent(out), dimension(0:2,0:n_max,0:l_max) :: Rnl
    real(dp), dimension(0:n_max,0:l_max+2) :: Lnl
    real(dp), dimension(0:n_max,0:l_max+2) :: Anl
    real(dp) :: Rnm1lp1, Rnm2lp2
    integer :: n,l
    Rnl = 0
    call NormalizationAll(n_max,l_max+2,Anl)
    call LaguerreAll(n_max,l_max+2,xi**2,Lnl)
    do n = 0,n_max
       do l = l_max,0,-1
          Rnl(0,n,l) = Anl(n,l)*xi**l*exp(-xi**2/2._dp)&
               *Lnl(n,l)/(b**(1.5_dp))
          if(n.gt.0) then
             if(l.lt.l_max) then
                Rnm1lp1 = Rnl(0,n-1,l+1)/Anl(n-1,l+1)
             else
                Rnm1lp1 = xi**(l+1)*exp(-xi**2/2._dp)*&
                     Lnl(n-1,l+1)/(b**(1.5_dp))
             endif
             if(n.gt.1) then
                if(l.lt.l_max-1) then
                   Rnm2lp2 = Rnl(0,n-2,l+2)/Anl(n-2,l+2)
                else
                   Rnm2lp2 = xi**(l+2)*exp(-xi**2/2._dp)*&
                        Lnl(n-2,l+2)/(b**(1.5_dp))
                endif
             else
                Rnm2lp2 = 0
             endif
          else
             Rnm1lp1 = 0
             Rnm2lp2 = 0
          endif
          Rnl(1,n,l) = ( (l/xi-xi)*Rnl(0,n,l) - 2*Anl(n,l)*Rnm1lp1)/b
          Rnl(2,n,l) = ( ((l/xi-xi)**2 - (l/xi**2+1))*Rnl(0,n,l) -&
               2*Anl(n,l)*((2*l+1)/xi-2*xi)*Rnm1lp1 + &
               4*Anl(n,l)*Rnm2lp2)/b**2
       enddo
    enddo
  end subroutine RadialHOALL

!> Calls the function \f$\texttt{HO\_Normalization(n,l)}\f$ for all values
!! of \f$n,l\f$ up to \f$n_{max},l_{max}\f$, and stores the result to
!! an array \f$\texttt{Anl(n,l)}\f$.
  subroutine NormalizationAll(n_max,l_max,Anl)
    implicit none
    integer, intent(in) :: n_max
    integer, intent(in) :: l_max
    real(dp), intent(out), dimension(0:n_max,0:l_max) :: Anl
    integer :: n,l
    do n = 0,n_max
       do l = 0,l_max
          Anl(n,l) = HO_Normalization(n,l)
       enddo
    enddo
  end subroutine NormalizationAll

!> Similar to \f$\texttt{LaguerreL}\f$, except that this computes ALL
!! Laguerre polynomials up to \f$\texttt{n\_max, l\_max}\f$ using the
!! recurrence relation (see section 1.2 of HO_basis.pdf).
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
