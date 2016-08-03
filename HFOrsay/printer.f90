!>Just a printing routine to the "hforsay.out" file.
subroutine printer(n,hfenergy,hfenergybcs,v2,esp,partnum,pr)
use constants
use bcs
use basis
use maths
implicit none
integer::n,k
logical::pr
double precision::hfenergy,hfenergybcs,v2(n),esp(n),partnum


open(22,file='hforsay.out')
write(22,*) "------------- System ----------------"
write(22,'(a,i5)') 'Number of neutrons......',Npart
write(22,*) "-------------- RESULTS --------------"
write(22,'(a,f16.9,a)') "Hartree-Fock Energy ....",hfenergy,' MeV'
write(22,'(a,f16.9,a)') "Particles Number........",partnum, '   '
if (flagbcs .eq. 1) then
write(22,'(a,f16.9,a)') "Pairing Gap.............",gap,'  '
write(22,'(a,f16.9,a)') "Pairing Energy..........",hfenergy - hfenergybcs,' Mev'
write(22,'(a,f16.9,a)') "Total Energy HF-BCS.....",hfenergybcs,' MeV'
endif
write(22,*)
if (pr) then
write(22,*) "--------- Single Particle Energies --------"
  write(22,'(a)') "   n  l  j          Energy(MeV)     Occ. Num"
  do k=1,red_size
  write(22,'(a,3i3,f18.6,a,f12.4)') ' ',n_red(k),l_red(k),j_red(k),esp(k),' MeV',v2(k)
  enddo
 write(22,*) "-------------------------------------------"
endif
close(22)
end subroutine

