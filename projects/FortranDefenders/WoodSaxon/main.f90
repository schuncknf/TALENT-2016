program main
    use grid
    use solver
    integer :: n,l,is,iq
    open(unit=5,file='in',status='old',form='formatted')
    open(unit=6,file='out',form='formatted')
    open(unit=13,file='plt',form='formatted')

    call init_params
    call init_grids
    call init_wavefunctions

    call solve_r
    write(6,*) "Single Particle States | "
    do iq =1,2
      if (iq == 1) write(6,*) "Neutrons | "
      if (iq == 2) write(6,*) "Protons  | "
      do n=1,lmax-2
        do l=0,lmax
          do is=1,2
            if (vocc(n,l,is,iq) > small) then
              write (6,*) "n=",n,"l=",l,"is=",is,"Energy=", energies(n,l,is,iq)
            end if
          end do
        end do
      end do
    end do
    write(6,*) "Total Bound States =", sum(vocc(1:lmax-2,:,:,:)), "Vpb(1) =", vpb(1)

    close(13)
    close(6)
    close(5)

end program main
