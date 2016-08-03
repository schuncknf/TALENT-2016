      program integrate
      real*8 x, rho, suma
      open(unit=1, file='density_T.dat', status='unknown')
      suma=0.d0
      do i=1, 2999
      read(1,*), x, rho
      suma= suma+rho*0.01d0
      enddo
      print*, suma
      end program integrate
