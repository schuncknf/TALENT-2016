program main
    use grid
    use solver
    integer :: iq,i,npr,ir,n,l
    real(wp) :: j
    open(unit=5,file='in',status='old',form='formatted')
    open(unit=6,file='out',form='formatted')
    open(unit=13,file='plt',form='formatted')
    open(unit=14,file='densities',form='formatted')

    call init_params
    call init_grids
    call init_wavefunctions

    call solve_r

    call energy_sort
    call build_densities
    write(6,*) "Single Particle States | "
    do iq =1,2

      if (iq == 1) then
      write(6,*) "Neutrons | "
      npr = nn
      else
      write(6,*) "Protons  | "
      npr = np
      end if
      do i = 1, npr
            if (sortenergies(i,iq) < - small) then
              j = sortstates(i,2,iq) + spin(sortstates(i,3,iq))
              if (sortstates(i,2,iq) == 0) j = 0
              write (6,*) "n=", sortstates(i,1,iq), &
                          &"l=",sortstates(i,2,iq),&
                          &"is=",sortstates(i,3,iq),&
                          &"Energy=", sortenergies(i,iq)

            end if
      end do
    end do

    write(6,*) "Total Neutrons =", sum(4*pi*h*meshpoints(:)**2*rho(:,1)), &
    "Total Protons =", sum(4*pi*h*meshpoints(:)**2 *rho(:,2))

    close(13)
    close(6)
    close(5)
    close(14)

end program main
