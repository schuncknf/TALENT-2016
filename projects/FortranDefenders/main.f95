program main
    use grid
    use solver

    open(unit=5,file='in',status='old',form='formatted')
    open(unit=6,file='out',form='formatted')
    open(unit=13,file='plt',form='formatted')

    call init_params
    call init_grids
    call init_wavefunctions

    call solve_r


    close(6)
    close(13)
    close(5)

end program main
