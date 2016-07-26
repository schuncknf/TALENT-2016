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
  do iq =1,2
    density(:,iq) = 0.
    do n=1,lmax-2
      do l=0,lmax
        do is=1,2
          if(energies(n,l,is,iq) < 0 .AND. energies(n,l,is,iq) > vpb(iq)+.01) then

            j = l + spin(is)
            if (l == 0) j = 0
            do ir = 0, nbox
              density(ir,iq) = density(ir,iq) + (2*j+1)*wfr(ir,n,l,is,iq) &
              *wfr(ir,n,l,is,iq) / (4*pi*meshpoints(ir)**2)
            end do
          end if
        end do
      end do
    end do
  end do
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
      density(:,iq)=0.
      do i = 1, npr
            if (sortenergies(i,iq) < - small) then
              j = sortstates(i,2,iq) + spin(sortstates(i,3,iq))
              if (sortstates(i,2,iq) == 0) j = 0
              do ir=0,nbox
                !density(ir,iq) = density(ir,iq) + (2*j+1)*wfr(ir,sortstates(i,1,iq),sortstates(i,2,iq),sortstates(i,3,iq),iq)&
                !*wfr(ir,sortstates(i,1,iq),sortstates(i,2,iq),sortstates(i,3,iq),iq) / (4*pi*meshpoints(ir)**2)
              end do
              write (6,*) "n=", sortstates(i,1,iq), &
                          &"l=",sortstates(i,2,iq),&
                          &"is=",sortstates(i,3,iq),&
                          &"Energy=", sortenergies(i,iq)

            end if
      end do
    end do

  write(6,*) "Total Neutrons =", sum(4*pi*h*meshpoints(:)**2 * density(:,1)), &
  "Total Protons =", sum(4*pi*h*meshpoints(:)**2 *density(:,2))
    do ir=0,nbox
      write (13,*) ir*h, wfr(ir,3,0,1,1), wfr(ir,3,1,1,1), wfr(ir,3,2,1,1)
      !wfr(ir,1,3,1,1), wfr(ir,1,4,1,1), wfr(ir,1,5,1,1), wfr(ir,1,6,1,1)
    end do
    close(13)
    close(6)
    close(5)
    close(14)

end program main
