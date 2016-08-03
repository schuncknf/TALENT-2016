!***********************************************************************************************************************************
!
!                                                               C G
!
!  Program:      CG
!
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!
!  Date:         April 20, 2005
!
!  Language:     Fortran-90
!
!  Version:      1.00a
!
!  Description:  Calculates Clebsch-Gordan cg_coefficients.
!
!  Files:        Source files:
!
!                   cg.f90                   Main program
!
!  Notes:
!
!***********************************************************************************************************************************

module cg
use maths

contains
    function cgordan(J1,J2,J3,M1,M2) result(cg_coeff)
    IMPLICIT NONE
    INTEGER :: I, K
    DOUBLE PRECISION :: J1, J2, J3, M1, M2, M3, C, SUMK, TERM, cg_coeff
    
!-----------------------------------------------------------------------------------------------------------------------------------

      M3 = M1 + M2

!
!     Check for invalid input.
!

      IF (ISFRAC(J1+J2+J3) .OR. ISFRAC(J1+M1)     .OR. ISFRAC(J2+M2) .OR.  &
          ISFRAC(J3+M3)    .OR. ISFRAC(-J1+J3-M2) .OR. ISFRAC(-J2+J3+M1)) THEN
         WRITE (UNIT=*, FMT='(/A)') ' Invalid input.'
         STOP
      END IF

!
!     Check for conditions that give C = 0.
!

      IF ( (J3 .LT. ABS(J1-J2)) .OR.  &
           (J3 .GT. (J1+J2))    .OR.  &
           (ABS(M1) .GT. J1)    .OR.  &
           (ABS(M2) .GT. J2)    .OR.  &
           (ABS(M3) .GT. J3)) THEN
         C = 0.0D0
      ELSE

!
!     Compute Clebsch-Gordan coefficient.
!

         cg_coeff = SQRT((J3+J3+1)/fac(NINT(J1+J2+J3+1)))
         cg_coeff = cg_coeff * SQRT(fac(NINT(J1+J2-J3))*FAC(NINT(J2+J3-J1))*FAC(NINT(J3+J1-J2)))
         cg_coeff = cg_coeff * SQRT(fac(NINT(J1+M1))*FAC(NINT(J1-M1))*fac(NINT(J2+M2))*fac(NINT(J2-M2))* &
         fac(NINT(J3+M3))*fac(NINT(J3-M3)))
         SUMK = 0.0D0
         DO K = 0, 99
            IF (J1+J2-J3-K .LT. 0.0D0) CYCLE
            IF (J3-J1-M2+K .LT. 0.0D0) CYCLE
            IF (J3-J2+M1+K .LT. 0.0D0) CYCLE
            IF (J1-M1-K    .LT. 0.0D0) CYCLE
            IF (J2+M2-K    .LT. 0.0D0) CYCLE
            TERM = fac(NINT(J1+J2-J3-K))*fac(NINT(J3-J1-M2+K))*fac(NINT(J3-J2+M1+K))*fac(NINT(J1-M1-K))*  &
               fac(NINT(J2+M2-K))*fac(K)
            IF (MOD(K,2) .EQ. 1) TERM = -TERM
            SUMK = SUMK + 1.0D0/TERM
         END DO
         cg_coeff = cg_coeff * SUMK
      END IF

    END function cgordan

!***********************************************************************************************************************************
!  ISFRAC
!
!  Return .TRUE. if the argument has a fractional part, and .FALSE. if the argument is an integer.
!***********************************************************************************************************************************

      FUNCTION ISFRAC (X) RESULT (Y)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: X
      LOGICAL :: Y
      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-8


      IF ((ABS(X)-INT(ABS(X))) .GT. EPS) THEN
         Y = .TRUE.
      ELSE
         Y = .FALSE.
      END IF

      RETURN
      END FUNCTION ISFRAC

end module cg
