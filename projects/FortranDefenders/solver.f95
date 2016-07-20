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

    Eupper = 10000_wp
    Elower = -v0
    do i=1,100000
      Etrial = (Eupper+Elower)/2.0
      ! Trying very large values right now, may change this if unstable

      ! Attempting to set the potential before hand, if this does not work, we
      ! can do it "on the fly"
      do ir=0,nbox
        potential(ir) = finitepot(ir,Etrial)
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
        write (6,*) "Converged1!"
        write (6,*) "|Eupper - Elower| =", abs(Eupper - Elower)
        write (6,*) "Energy =", Etrial
        write (6,*) "Exact Energy =",finite_exact()
        write (6,*) "Difference between calculated and exact =",Etrial-finite_exact()
        exit
      end if
    end do

    if (abs(Eupper - Elower) > conv) then
      write (6,*) "Program did not converge!"
      write (*,*) abs(Eupper - Elower)
    end if
    ! Normalisation
    norm = sqrt(sum(h*wavefunctions(:)*wavefunctions(:)))
    wavefunctions(:) = wavefunctions(:)/norm

    ! Printing points for plotting. I run $ xmgrace plt

    do ir=0,nbox
      test(ir) = sqrt(2/(nbox*h)) * sin((nodes+1) * pi * meshpoints(ir)/(nbox*h))
      write(13,*) ir*h, wavefunctions(ir)
    end do

  end subroutine solve


  subroutine solvelr
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine is closely based on the notes provided by the organizers
    ! of the 2016 Density Functional Theory TALENT Course.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: i, j, ir, nnodes
    real(wp) :: Etrial, Eupper, Elower, a1, a2, a3, norm
    real(wp), allocatable :: potential(:), test(:), test2(:)
    logical :: sign

    allocate(potential(0:nbox),test(0:nbox),test2(0:nbox))

    Eupper = 100_wp
    Elower = -v0
    do i=1,10000
      Etrial = (Eupper+Elower)/2.0
      ! Trying very large values right now, may change this if unstable

      ! Attempting to set the potential before hand, if this does not work, we
      ! can do it "on the fly"
      do ir=0,nbox
        potential(ir) = finitepot(ir,Etrial)
      end do

      wfl(0) = 0.0
      wfr(nbox) = 0.0
      wfl(1) = 1.0
      wfr(nbox-1) = 1.0
      do ir=1,nbox-1
        a1 = 2.0 * (1.0 - (5.0/12.0) * h**2 * potential(ir))
        !a1 = 2.0 * (1.0- 5.0*h**2/12.0 * potential(ir))
        a2 = (1.0 + (1.0/12.0) * h**2 * potential(ir-1))
        a3 = (1.0 + (1.0/12.0) * h**2 * potential(ir+1))

        wfl(ir+1) = (a1*wfl(ir) - a2*wfl(ir-1))/a3
      end do
      a1 = 0
      a2 = 0
      a3 = 0
      do ir=nbox-1,1,-1
        a1 = 2.0 * (1.0 - (5.0/12.0) * h**2 * potential(ir))
        !a1 = 2.0 * (1.0- 5.0*h**2/12.0 * potential(ir))
        a2 = (1.0 + (1.0/12.0) * h**2 * potential(ir+1))
        a3 = (1.0 + (1.0/12.0) * h**2 * potential(ir-1))

        wfr(ir-1) = (a1*wfr(ir) - a2*wfr(ir+1))/a3

      end do

      if (.NOT. (mod(nodes,2) .EQ. 0)) then
        wfr(:) = -wfr(:)
      end if

      nnodes = 0
      ! I believe this is appropriate logic for checking the sign change of the
      ! wavefunction
      sign = wfr(nbox-1) > 0
      do ir =nbox-1,1,-1
        if ((wfr(ir) > 0) .NEQV. sign) then
          sign = wfr(ir) > 0
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
        Elower = Elower + .01

      end if


      if (abs(Eupper - Elower) < conv) then
        !if (abs(wfr(nbox/2) - wfl(nbox/2)) < conv) then
          !if (abs((wfr(nbox/2+1)-wfr(nbox/2))/h - (wfl(nbox/2)-wfl(nbox/2-1))/h) < conv) then
            write (6,*) abs(wfr(nbox/2) - wfl(nbox/2)), abs((wfr(nbox/2+1)-wfr(nbox/2))/h - (wfl(nbox/2)-wfl(nbox/2-1))/h)
            write (6,*) "Converged!"
            write (6,*) "|Eupper - Elower| =", abs(Eupper - Elower)
            write (6,*) "Energy =", Etrial
            write (6,*) "Exact Energy =",finite_exact()
            write (6,*) "Difference between calculated and exact =",Etrial-finite_exact()
            exit
          !end if
        !end if
      end if
    end do

    if (abs(Eupper - Elower) > conv) then
      write (6,*) "Program did not converge!"
      write (*,*) abs(Eupper - Elower)
    end if

    wavefunctions(0:nbox/2-1) = wfl(0:(nbox/2-1))
    wavefunctions(nbox/2:nbox) = wfr(nbox/2:nbox)
    ! Normalization
    norm = sqrt(sum(h*wavefunctions(:)*wavefunctions(:)))
    wavefunctions(:) = wavefunctions(:)/norm

    ! Printing points for plotting. I run $ xmgrace plt

    do ir=0,nbox
      test(ir) = sqrt(2/(nbox*h)) * sin((nodes+1) * pi * meshpoints(ir)/(nbox*h))
      test2(ir) = woodsaxon(ir)
      write(13,*) ir*h, wavefunctions(ir),test2(ir)
    end do

  end subroutine solvelr

  function infwell_exact() result(energy)
    real(wp) :: energy

    energy = (nodes+1)**2 *pi**2 *hbar22m/(nbox*h)**2

  end function

  function finite_exact() result(energy)
    real(wp) :: energy

    energy = (nodes+1)**2 *pi**2 *hbar22m/(nbox*h)**2 -v0

  end function

  function infwell(ir, Etrial) result(pot)
    real(wp) :: pot
    real(wp), intent(in) :: Etrial
    integer, intent(in) :: ir

    if ((ir .EQ. 0) .OR. (ir .EQ. nbox)) then
      pot = 1E20_wp
    else
      pot = Etrial/hbar22m
    end if

  end function

  function finitepot(ir, Etrial) result(pot)
    real(wp) :: pot
    real(wp), intent(in) :: Etrial
    integer, intent(in) :: ir

    if ((ir <= nbox/4) .OR. (ir >= nbox*3/4)) then
      pot = Etrial/hbar22m
    else
      pot = (Etrial+v0)/hbar22m
    end if

  end function

 !function woodsaxon(ir) result(pot)
  !  real(wp) :: pot
  !  integer, intent(in) :: ir

!    if (ir < nbox/2) then
!      pot = v0 * 1 / (1 * exp((-meshpoints(nbox/2-ir)+meshpoints!!(radius))/0.67))
!    else
!      pot = v0 * 1 / (1 * exp((-meshpoints(ir-nbox/2)+meshpoints(radius))/0.67))
!    end if


!  end function

 function woodsaxons( r0, r ) result(res)
    implicit none
    real (wp), intent(in) :: r0, r
    real (wp) :: res
    !
    res =  1 / ( 1 + exp( (   r - r0 ) / 0.67 ) )    &
         + 1 / ( 1 + exp( ( - r - r0 ) / 0.67 ) ) - 1
    !
  end function woodsaxons

  function woodsaxon(ir) result(pot)
    real(wp) :: pot

    integer, intent(in) :: ir

  !pot = -v0 * 1 / (1 + exp((meshpoints(ir)-300*h)/0.67))
    if (ir < nbox/2) then
      pot =-v0 * 1 / (1 + exp((+meshpoints(nbox/2-ir)-300*h)/0.67))
    else
      pot = -v0 * 1 / (1 + exp((+meshpoints(ir-nbox/2)-300*h)/0.67))
  end if


  end function

end module solver
