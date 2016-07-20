!----------------------------------------------------------------------
!Module: types
!---------------------------------------------------------------------
!> Defines the integer parameters <B>sp</B>, <B>dp</B> and <B>qp</B> to
!! be used as kinds to define real variables as single precision 
!! (32-bit), double precision (64-bit) and quadruple precision 
!! (128-bit) (depending on the machine and compiler, this last one may
!! not always be available).
!!
!!The intrinsic module iso_fortran_env (Fortran 2008 and later) is
!!used.
!!
!!More types may be defined later (i.e. larger integers)
!!
!!Fundamental constants (i.e. &pi;, e, ...) may also be defined here if
!!desired. This is a good place to define this type of constants as
!!all modules and the main program will (in principle) use this module
!!
!! @param sp single precision kind 
!! @param dp double precision kind 
!! @param qp quadruple precision kind (not available in every machine)
!! @param pi \f$ \pi = 3.141592 \ldots \f$ 
!! @author Rodrigo Navarro Perez
!---------------------------------------------------------------------
module types
  use, intrinsic :: iso_fortran_env
  implicit none
  integer, parameter :: sp = REAL32 !< single precision kind
  integer, parameter :: dp = REAL64 !< double precision kind
  integer, parameter :: qp = REAL128!< quadruple precision kind
  real(dp), parameter :: pi=acos(-1._dp)!< &pi; = 3.141592...
end module types
