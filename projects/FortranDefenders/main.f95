program main
        use grid

        open(unit=5,file='in',status='old',form='formatted')
        open(unit=6,file='out',form='formatted')

        call init_params
        call init_grids

        write (6,*) meshpoints(:)

        close(6)
        close(5)
end program main
