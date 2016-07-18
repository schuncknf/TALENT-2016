SUBROUTINE PRINT_EIGENVALUES( DESC, N, WR, WI )
      CHARACTER*(*)    DESC
      INTEGER          N,NIN
      DOUBLE PRECISION WR( * ), WI( * )
      DOUBLE PRECISION ZERO
      PARAMETER        ( ZERO = 0.0 )
      INTEGER          J
      NIN = N 
      if (N > 10) then
      WRITE(*,*) '10 First Eigenvalues'
      NIN = 10
      endif
      WRITE(*,*) DESC
      DO J = 1, NIN
         IF( WI( J ).EQ.ZERO ) THEN
            WRITE(*,9998,ADVANCE='NO') WR( J )
         ELSE
            WRITE(*,9999,ADVANCE='NO') WR( J ), WI( J )
         END IF
      END DO
      WRITE(*,*)
 9998 FORMAT( 11(:,1X,F6.2) )
 9999 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
      RETURN
      END
