module solver
  use grid
  implicit none

contains

  subroutine solve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine is closely based on the notes provided by the organizers
    ! of the 2016 Density Functional Theory TALENT Course.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: i, ir, nnodes
    real(wp) :: Etrial, Eupper, Elower, a1, a2, a3, norm
    real(wp), allocatable :: potential(:), test(:)
    logical :: sign

    allocate(potential(0:nbox),test(0:nbox))

    Eupper = 100_wp
    Elower = 0_wp
    do i=1,10000
      Etrial = (Eupper+Elower)/2.0
      ! Trying very large values right now, may change this if unstable
      potential(0) = 1E20_wp
      potential(nbox) = 1E20_wp
      ! Attempting to set the potential before hand, if this does not work, we
      ! can do it "on the fly"
      do ir=1,nbox-1
        potential(ir) = Etrial/hbar22m
      end do

      wavefunctions(0) = 0.0
      wavefunctions(1) = 1.0
      do ir=1,nbox-1
        a1 = 2.0 * (1.0 - (5.0/12.0) * h**2 * potential(ir))
        !a1 = 2.0 * (1.0- 5.0*h**2/12.0 * potential(ir))
        a2 = (1.0 + (1.0/12.0) * h**2 * potential(ir-1))
        a3 = (1.0 + (1.0/12.0) * h**2 * potential(ir+1))

        wavefunctions(ir+1) = (a1*wavefunctions(ir) - a2*wavefunctions(ir-1))/a3
      end do
      nnodes = 0
      ! I believe this is appropriate logic for checking the sign change of the
      ! wavefunction
      sign = wavefunctions(1) > 0
      do ir =1,nbox
        if ((wavefunctions(ir) > 0) .NEQV. sign) then
          sign = wavefunctions(ir) > 0
          nnodes = nnodes + 1
        end if
      end do

      if (nnodes > nodes) then
        Eupper = Etrial
      else if (nnodes < nodes) then
        Elower = Etrial
      else
        ! This is a variation on the lower energy in order to "squeeze" the
        ! energies together. by moving it a small amount (arbitrarily here)
        ! we can force the solution to converge.
        Elower = Elower + 0.01
      end if

      if (abs(Eupper - Elower) < conv) then
        write (6,*) "Converged!"
        write (6,*) "|Eupper - Elower| =", abs(Eupper - Elower)
        write (6,*) "Energy =", Etrial
        write (6,*) "Exact Energy =",(nodes+1)**2 *pi**2 *hbar22m/(nbox*h)**2
        write (6,*) "Difference between calculated and exact =",Etrial-(nodes+1)**2 *pi**2 *hbar22m/(nbox*h)**2
        exit
      end if
    end do

    if (abs(Eupper - Elower) > conv) then
      write (6,*) "Program did not converge!"
      write (*,*) abs(Eupper - Elower)
    end if
    ! Normalisation
    norm = sum(abs(wavefunctions(:)))
    do ir = 0,nbox
      wavefunctions(ir) = wavefunctions(ir)/norm
    end do
    norm = sum(abs(wavefunctions(:)))
    write (*,*) norm
    ! Printing points for plotting. I run $ xmgrace plt

    do ir=0,nbox
      test(ir) = sqrt(2/(nbox*h)) * sin((nodes+1) * pi * meshpoints(ir)/(nbox*h))
      write(13,*) ir*h, wavefunctions(ir), test(ir)
    end do

  end subroutine solve

end module solver
