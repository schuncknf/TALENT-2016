program main
    use grid
    use solver
    integer :: iq,i,npr,j
    open(unit=5,file='in',status='old',form='formatted')
    open(unit=6,file='out',form='formatted')
    open(unit=13,file='plt',form='formatted')

    call init_params
    call init_grids
    call init_wavefunctions

    call solve_r
    call energy_sort
   ! do i = 1,nmax
   !  print *, sortenergies(i,1), sortstates(i,1,1), &
   !  &       sortstates(i,2,1),sortstates(i,3,1)
   ! end do
!    write(6,*) "Single Particle States | "
!    do iq =1,2
!      if (iq == 1) write(6,*) "Neutrons | "
!      if (iq == 2) write(6,*) "Protons  | "
!      do n=1,lmax-2
!        do l=0,lmax
!          do is=1,2
!            if (vocc(n,l,is,iq) > small) then
!              write (6,*) "n=",n,"l=",l,"is=",is,"Energy=", energies(n,l,is,iq)
!            end if
!          end do
!        end do
!      end do
!    end do
!    write(6,*) "Total Bound States =", sum(vocc(1:lmax-2,:,:,:)), "Vpb(1) =", vpb(1)
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
              j = (2*sortstates(i,2,iq)+1)
              if (j >= npr) exit
              write (6,*) "n=", sortstates(i,1,iq), &
                          &"l=",sortstates(i,2,iq),&
                          &"is=",sortstates(i,3,iq),&
                          &"Energy=", sortenergies(i,iq)
            end if
      end do
    end do
  write(6,*) "Total Bound States =", sum(vocc(1:lmax-2,:,:,:)), "Vpb(1) =", vpb(1)

    close(13)
    close(6)
    close(5)

end program main
