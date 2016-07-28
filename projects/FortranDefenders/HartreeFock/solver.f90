module solver
  use grid
  implicit none

  real(wp), allocatable :: vocc(:,:,:,:), energies(:,:,:,:), &
                         &  sortenergies(:,:)
  integer, allocatable :: sortstates(:,:,:)

contains

  subroutine solve_r
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine is closely based on the notes provided by the organizers
    ! of the 2016 Density Functional Theory TALENT Course.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: i, ir, nnodes,l, is, iq, n
    real(wp) :: Etrial, Eupper, Elower, a1, a2, a3, b1, b2, b3, norm, j, diff
    real(wp), allocatable :: potential(:), woodsaxon(:), woodsaxond(:), spinorbitmat(:), coulombmat(:)

    allocate(potential(0:nbox),vocc(lmax,0:lmax,2,2),energies(lmax,0:lmax,2,2),&
    woodsaxon(0:nbox),woodsaxond(0:nbox),spinorbitmat(0:nbox),coulombmat(0:nbox))
    ! Storing the Woodsaxon and derivative for future usage
    do ir = 0,nbox
      woodsaxon(ir) = fullwoodsaxon(ir)
      spinorbitmat(ir) = r0**2 * dfullwoodsaxon(ir)*woodsaxon(ir)/meshpoints(ir)
      coulombmat(ir) = coulomb(ir)
    end do
    spinorbitmat(0) = 0.0
    wfr(:,:,:,:,:) = 0.0
    !Main loop begins here; be careful
    do iq =1,2
      do n =1,lmax-2
        do l =0,lmax
            do is = 1,2
                j = l + spin(is)
                if (l==0) j=0.5
                ! Bound States only
                Eupper = 1_wp
                Elower = vpb(iq)
                do i=1,1000000
                  ! Trial Energy for Numerov Algorithm
                  Etrial = (Eupper+Elower)/2.0
                  do ir=0,nbox
                    ! Isospin dependent potential using matrices from before
                    if (iq .EQ. 1) then
                       potential(ir) = (-vpb(iq)*woodsaxon(ir)-vso*spinorbitmat(ir)*0.5*(j*(j+1) - l*(l+1) - 0.75) &
                      -hbar22m*l*(l+1)/meshpoints(ir)**2+Etrial)/hbar22m
                    else
                       potential(ir) = (-vpb(iq)*woodsaxon(ir)-vso*spinorbitmat(ir)*0.5*(j*(j+1) - l*(l+1) - 0.75) &
                      -hbar22m*l*(l+1)/meshpoints(ir)**2+Etrial-coulombmat(ir))/hbar22m
                    end if
                  end do

                  ! Setting initial values for left and right W.F.
                  wfr(nbox,n,l,is,iq) = 0.
                  wfr(nbox-1,n,l,is,iq) = 1.
                  wfl(0,n,l,is,iq) = 0.
                  wfl(1,n,l,is,iq) = meshpoints(1)**(l+1)

                  nnodes = 0
                  ! Our Friend, Numerov
                  do ir=nbox-1,1,-1
                    ! These are for the right propagating W.F.'s
                    a1 = 2.0 * (1.0 - (5.0/12.0) * h**2 * potential(ir))
                    a2 = (1.0 + (1.0/12.0) * h**2 * potential(ir+1))
                    a3 = (1.0 + (1.0/12.0) * h**2 * potential(ir-1))
                    ! These are for the left propagating W.F.'s
                    b1 = 2.0 * (1.0 - (5.0/12.0) * h**2 * potential(nbox-ir))
                    b2 = (1.0 + (1.0/12.0) * h**2 * potential(nbox - ir-1))
                    b3 = (1.0 + (1.0/12.0) * h**2 * potential(nbox - ir+1))

                    wfr(ir-1,n,l,is,iq) = (a1*wfr(ir,n,l,is,iq) - a2*wfr(ir+1,n,l,is,iq))/a3
                    wfl(nbox - ir+1,n,l,is,iq) = (b1*wfl(nbox - ir,n,l,is,iq) - b2*wfl(nbox - ir-1,n,l,is,iq))/b3
                    ! Setting coefficient for left W.F. scaling W.R.T. right
                    if (ir == njoin) diff = wfr(ir,n,l,is,iq) / wfl(ir,n,l,is,iq)
                    ! Scale all the left W.F. from njoin -> 0
                    if (ir <= njoin) then
                      wfl(ir,n,l,is,iq) = diff * wfl(ir,n,l,is,iq)
                    end if
                    ! Count the nodes (Except those spurious ones near the boundary)
                    if((wfl(nbox-ir,n,l,is,iq)*wfl(nbox-ir+1,n,l,is,iq) < 0)) nnodes = nnodes+1
                  end do
                  ! Unite the left and right
                  wfr(0:njoin,n,l,is,iq) = wfl(0:njoin,n,l,is,iq)

                  if (nnodes > n-1) then
                    Eupper = Etrial
                  else if (nnodes <= n-1) then
                    Elower = Etrial
                  end if
                  ! If the energies converge on something that is not vpb, select that solution
                  if (abs(Eupper - Elower) < conv) then
                    if (Etrial < 0 .AND. Etrial > vpb(iq)+.01) then
                          vocc(n,l,is,iq) = 2*l+1
                          energies(n,l,is,iq) = etrial
                          norm = sqrt(sum(h*wfr(:,n,l,is,iq)*wfr(:,n,l,is,iq)))
                          wfr(:,n,l,is,iq) = wfr(:,n,l,is,iq)/norm
                       end if
                    ! Leave this place
                    exit
                  end if
                end do
              end do
          end do
        end do
      end do
  end subroutine solve_r

  subroutine ddensities
    integer :: ir,iq
    do iq=1,2
      do ir=0,nbox
        if(ir < 1) then
          djsc(ir,iq) = (-jsc(ir+2,iq) + 8*jsc(ir+1,iq) &
                 +8*jsc(ir+1,iq) - jsc(ir+2,iq))/(12*h)

          drho(ir,iq) = (-rho(ir+2,iq) + 8*rho(ir+1,iq) &
                 -8*rho(ir+1,iq) + rho(ir+2,iq))/(12*h)

          ddrho(ir,iq) = (-rho(ir+2,iq)+16*rho(ir+1,iq)-30*rho(ir,iq)&
                  +16*rho(ir+1,iq)-rho(ir+2,iq))/(12*h**2)
        else if (ir < 2) then
          djsc(ir,iq) = (-jsc(ir+2,iq) + 8*jsc(ir+1,iq) &
                 +8*jsc(ir-1,iq) - jsc(ir,iq))/(12*h)

          drho(ir,iq) = (-rho(ir+2,iq) + 8*rho(ir+1,iq) &
                 -8*rho(ir-1,iq) + rho(ir,iq))/(12*h)

          ddrho(ir,iq) = (-rho(ir+2,iq)+16*rho(ir+1,iq)-30*rho(ir,iq)&
                  +16*rho(ir-1,iq)-rho(ir,iq))/(12*h**2)
        else if ((ir >= 2) .AND. (ir <= nbox-2)) then
          drho(ir,iq) = (-rho(ir+2,iq) + 8*rho(ir+1,iq) &
                 -8*rho(ir-1,iq) + rho(ir-2,iq))/(12*h)

          djsc(ir,iq) = (-jsc(ir+2,iq) + 8*jsc(ir+1,iq) &
                 -8*jsc(ir-1,iq) + jsc(ir-2,iq))/(12*h)

          ddrho(ir,iq) = (-rho(ir+2,iq)+16*rho(ir+1,iq)-30*rho(ir,iq)&
                  +16*rho(ir-1,iq)-rho(ir-2,iq))/(12*h**2)


        else if ((ir > nbox-2) .AND. (ir/=nbox)) then
          drho(ir,iq) = 0.
          ddrho(ir,iq) = 0.
          djsc(ir,iq) = 0.
        else
          drho(ir,iq) = 0.
          ddrho(ir,iq) = 0.
          djsc(ir,iq) = 0.
        end if
      end do
      ddrho(0:3,iq)=ddrho(4,iq)
    end do

    drho(:,3) = drho(:,1) + drho(:,2)
    drho(:,4) = drho(:,1) - drho(:,2)
    djsc(:,3) = djsc(:,1) + djsc(:,2)
    djsc(:,4) = djsc(:,1) - djsc(:,2)
    ddrho(:,3) = ddrho(:,1) + ddrho(:,2)
    ddrho(:,4) = ddrho(:,1) - ddrho(:,2)

  end subroutine ddensities

  function dwavefunction(ir,n,l,is,iq) result(derv)
    integer, intent(in) :: ir,n,l,is,iq
    real(wp) :: derv

    if(ir == 0) then
      if(mod(l,2)==1) then
        derv = (-wfr(ir+2,n,l,is,iq) + 8*wfr(ir+1,n,l,is,iq) &
               -8*wfr(ir+1,n,l,is,iq) + wfr(ir+2,n,l,is,iq))/(12*h)
      else
        derv = (-wfr(ir+2,n,l,is,iq) + 8*wfr(ir+1,n,l,is,iq) &
               +8*wfr(ir+1,n,l,is,iq) - wfr(ir+2,n,l,is,iq))/(12*h)
      end if
    else if (ir == 1) then
      if(mod(l,2)==1) then
        derv = (-wfr(ir+2,n,l,is,iq) + 8*wfr(ir+1,n,l,is,iq) &
               -8*wfr(ir-1,n,l,is,iq) + wfr(ir,n,l,is,iq))/(12*h)
      else
        derv = (-wfr(ir+2,n,l,is,iq) + 8*wfr(ir+1,n,l,is,iq) &
               +8*wfr(ir-1,n,l,is,iq) - wfr(ir,n,l,is,iq))/(12*h)
      end if
    else if ((ir >= 2) .AND. (ir <= nbox-2)) then
      derv = (-wfr(ir+2,n,l,is,iq) + 8*wfr(ir+1,n,l,is,iq) &
             -8*wfr(ir-1,n,l,is,iq) + wfr(ir-2,n,l,is,iq))/(12*h)
    else if ((ir > nbox-2) .AND. (ir/=nbox)) then
      derv = 0.!(2*wfr(ir+1,n,l,is,iq) + 3*wfr(ir,n,l,is,iq) &
             !-6*wfr(ir-1,n,l,is,iq) + wfr(ir-2,n,l,is,iq))/(6*h)
    else
      derv = 0.
    end if

  end function


  subroutine energy_sort
    integer :: n, l, iq, is,k,i,nfill,nfull
    real(wp) :: temp,j
    integer, dimension(1:3) :: state
    allocate(sortenergies(1:nmax,2),sortstates(1:nmax,1:3,2))
    ! This sorts the energies
    sortenergies = small
    sortstates = small
    do iq =1,2
     nfull = nn
     if (iq == 2) nfull = np
     k = 1
     nfill = 0
     do while(nfill.lt.nfull)
        temp = 0._wp
        state = 1
       do n = 1,lmax
        do l = 0,lmax
          do is=1,2
            if(temp > energies(n,l,is,iq)) then
            temp = energies(n,l,is,iq)
            state(1) = n
            state(2) = l
            state(3) = is
            end if
          end do
        end do
       end do
       sortenergies(k,iq) = energies(state(1),state(2),state(3),iq)
       do i = 1,3
         sortstates(k,i,iq) = state(i)
       end do
       energies(state(1),state(2),state(3),iq) = 0.0_wp
       if (state(2)==0) energies(state(1),state(2),:,iq) = 0.0_wp
         j = state(2) + spin(state(3))
       if (state(2)==0) j = 0.5
       nfill = nfill + 2*j+1
       k = k+1
     end do
  end do

  end subroutine energy_sort

  subroutine build_fields
  integer :: ir, iq, ir2
  real(wp), dimension(0:nbox,2) :: ucnew,umrnew,uddnew,usonew
  real(wp), dimension(0:nbox) :: ucoulnew
  real :: tot1=0.0d0,tot2=0.0d0
  real :: xmix, ymix

  xmix = 0.4
  ymix = 1.-xmix

  do iq = 1,2
   do ir = 0,nbox
  !!Central Field U(r)
           ucnew(ir,iq) = ucnew(ir,iq) + &
                & 2*(a0r0-a1r1)*rho(ir,3) + 4*a1r1 * rho(ir,iq)  &
               	& + (a0tau0-a1tau1) *tau(ir,3)+ 2 *a1tau1*tau(ir,iq) &
               	& + 2*( a0r0p-a1r1p )*ddrho(ir,3) + 4 *a1r1p * ddrho(ir,iq)
  !! Part of U(r) coming from d(M(r))
           umrnew(ir,iq) = umrnew(ir,iq)  &
                & + (cso0-cso1 ) *(djsc(ir,3) + 2 * jsc(ir,3)/meshpoints(ir) ) &
                & + 2 *cso1 * ( djsc(ir,iq) + 2 * jsc(ir,iq) / meshpoints(ir) )
  !! t3 part of U(r)
           uddnew(ir,iq) = uddnew(ir,iq) &
                & + ( 2 + sig ) * (cddr0-cddr1)*rho(ir,3)**(sig+1)  &
                & + 2*sig*cddr1*(rho(ir,1)**2+rho(ir,2)**2)*rho(ir,3)**(sig-1) &
                & + 4 * cddr1 * rho(ir,iq) * rho(ir,3)**sig
   !!spin-orbit part
           usonew(ir,iq) = usonew(ir,iq) &
                & - (cso0-cso1 )*drho(ir,3)/meshpoints(ir) &
                & - 2 *cso1 * drho(ir,iq) / meshpoints(ir)
           if (j2terms) then
           usonew(ir,iq) = usonew(ir,iq)&
                & -(a0t0-a1t1) *jsc(ir,3) / meshpoints(ir) &
               	& - 2 *a1t1 * jsc(ir,iq) / meshpoints(ir)
           end if
   !!coulomb
    if (iq==2) then
      do ir2=0,ir
       tot1=tot1+rho(ir2,2)*(meshpoints(ir)**2)
      end do
      do ir2=ir,nbox
       tot2=tot2+rho(ir2,2)*meshpoints(ir)
      end do
      ucoulnew(ir)=4.0d0*pi*e2*(tot1/meshpoints(ir) + tot2)
    end if
   end do
  end do

  uc = ucnew(:,:)*xmix + uc(:,:)*ymix
  umr = umrnew(:,:)*xmix + umr(:,:)*ymix
  udd = uddnew(:,:)*xmix + udd(:,:)*ymix
  uso = usonew(:,:)*xmix + uso(:,:)*ymix
  ucoul = ucoulnew(:)*xmix + ucoul(:)*ymix


  end subroutine build_fields

  subroutine build_densities
  integer :: npr,iq,ir,i,l,n,is
  real(wp) :: j

  do iq =1,2

      if (iq == 1) then
      npr = nn
      else
      npr = np
      end if
      rho(:,iq)=0.
      tau(:,iq)=0.
      jsc(:,iq)=0.
      do i = 1, npr
         if (sortenergies(i,iq) < - small) then
           n = sortstates(i,1,iq)
           l = sortstates(i,2,iq)
           is= sortstates(i,3,iq)
           j = sortstates(i,2,iq) + spin(sortstates(i,3,iq))
           if (sortstates(i,2,iq) == 0) j = 0.5
           do ir=1,nbox
             tau(ir,iq) = tau(ir,iq) + (2*j+1)*((dwavefunction(ir,n,l,is,iq)&
             -wfr(ir,n,l,is,iq)/meshpoints(ir))**2+l*(l+1)*(wfr(ir,n,l,is,iq)**2)&
             /meshpoints(ir)**2)/(4*pi*meshpoints(ir)**2)

             rho(ir,iq) = rho(ir,iq) + (2*j+1)*wfr(ir,sortstates(i,1,iq),sortstates(i,2,iq),sortstates(i,3,iq),iq)&
             *wfr(ir,sortstates(i,1,iq),sortstates(i,2,iq),sortstates(i,3,iq),iq) / (4*pi*meshpoints(ir)**2)

             jsc(ir,iq) = jsc(ir,iq) + (2*j+1)*(j*(j+1)-sortstates(i,2,iq)*(sortstates(i,2,iq)+1)-0.75)&
             *wfr(ir,sortstates(i,1,iq),sortstates(i,2,iq),sortstates(i,3,iq),iq)**2&
             /(4*pi*meshpoints(ir)**3)
            end do
            rho(0,iq) = rho(1,iq)
            tau(2,iq) = tau(3,iq)
            tau(1,iq) = tau(2,iq)
            tau(0,iq) = tau(1,iq)


          end if
      end do
   end do
   rho(:,3)=rho(:,1) + rho(:,2)
   rho(:,4)=rho(:,1) - rho(:,2)
   tau(:,3)=tau(:,1) + tau(:,2)
   tau(:,4)=tau(:,1) - tau(:,2)
   jsc(:,3)=jsc(:,1) + jsc(:,2)
   jsc(:,4)=jsc(:,1) - jsc(:,2)

   call ddensities
    do ir = 0,nbox
      write(14,*) ir*h, rho(ir,1),rho(ir,2),rho(ir,3),drho(ir,1),ddrho(ir,1),tau(ir,1),tau(ir,2),jsc(ir,1)!,rho(ir,4)
    end do
  end subroutine build_densities

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

    if ((ir <= nbox/2 - radius/2) .OR. (ir >= nbox/2 + radius/2)) then
      pot = Etrial/hbar22m
    else
      pot = (Etrial+v0)/hbar22m
    end if

  end function


  function woodsaxon_s(ir, Etrial) result(pot)
    real(wp) :: pot
    real(wp), intent(in) :: Etrial
    integer, intent(in) :: ir

  !pot = -v0 * 1 / (1 + exp((meshpoints(ir)-300*h)/0.67))
    if (ir < nbox/2) then
      pot = (Etrial + v0 * 1 / (1 + exp((meshpoints(nbox/2-ir)-radius/2*h)/0.67)))/hbar22m
    else
      pot = (Etrial + v0 * 1 / (1 + exp((meshpoints(ir-nbox/2)-radius/2*h)/0.67)))/hbar22m
  end if


  end function

  function woodsaxon(ir, Etrial) result(pot)
    real(wp) :: pot
    integer, intent(in) :: ir
    real(wp), intent(in) :: Etrial

    if (ir < (nbox/2 - radius)) then
      pot = Etrial/hbar22m
    else
      pot = ( v0 * 1 / (1 + exp((meshpoints(ir-(nbox/2 - radius))-radius*h)/0.67)))/ hbar22m
    end if
  end function


  function spinorbit(ir,l,is) result(pot)
    integer, intent(in) :: ir,l, is
    real(wp) ::  pot

    if (ir .EQ. 0) then
      pot = 0._wp
    else
      if (l .EQ. 0) then
        pot = 0._wp
      else
        pot = r0**2 * dfullwoodsaxon(ir)/meshpoints(ir)* 0.5 *((l+spin(is))*(l+spin(is)+1) - l*(l+1) - 0.75)
      end if
    end if

  end function

  !function coulomb(ir) result(pot)
  !  integer,intent(in) :: ir
  !  integer :: ir2
  !  real(wp) ::pot
  !  real :: tot1=0.0d0,tot2=0.0d0
  !       DO ir2=0,ir
  !	tot1=tot1+rho(ir2,2)*(meshpoints(ir2)**2)
  !	ENDDO
  !	DO ir2=ir,nbox
  !	tot2=tot2+rho(ir2,2)*meshpoints(ir2)
  !	ENDDO
  !	pot=4.0d0*pi*e2*(tot1/meshpoints(ir) + tot2)
  !end function
    function coulomb(ir) result(pot)
    integer,intent(in) :: ir
    real(wp) ::pot
        if(ir*h .lt. nrad ) pot= (np*e2/(2*nrad))*(3.0d0- (ir*h/nrad)**2)
        if(ir*h .ge. nrad ) pot= np*e2/(ir*h)

  end function
end module solver
