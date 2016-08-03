module constants
!> This module is keeping all the needed physical and mathematical constants of the problem
implicit none
double precision,parameter::one = 1.d0,two=2.d0,half=0.5d0,zero=0.d0
double precision,parameter::pi=acos(-1.d0)
double precision,parameter::hc=197.327d0 ! (MeV*fm)
double precision,parameter::mneutron =939.57d0
!definition of the parameters of the Minesotta potential
double precision,parameter::kr=1.487d0 !in fm -2
double precision,parameter::kt=0.639d0 !in fm -2
double precision,parameter::ks=0.465d0 !in fm -2
double precision,parameter::v0r=200.00d0 !in MeV
double precision,parameter::v0t=178.00d0 !in MeV
double precision,parameter::v0s=91.85d0 !in MeV
double precision,parameter::ama = 20.736209412d0 !in (fm**2)
double precision,parameter::mc2 = 938.90590d0 !in MeV
integer::nbase,npart,maxit,ngauss,n_lines,ntx,flagbcs
integer::na_max,la_max,ja_max,nb_max,lb_max,jb_max,iplot,flagext
double precision::homega,bosc,g_pair
integer,allocatable::exttag(:,:,:)
contains
!> The subroutine used to read the input file
subroutine reader()
implicit none
integer,allocatable:: nr(:),nl(:),nj(:)
open(1,file='hforsay.dat',status='old')
read(1,'(10x,i5)') nbase
read(1,'(10x,i5)') ngauss
read(1,'(10x,i5)') npart
read(1,'(10x,i5)') maxit
read(1,'(10x,f10.4)') homega
read(1,'(10x,i5)') flagbcs
read(1,'(10x,f10.4)') g_pair
read(1,'(10x,i5)') iplot
read(1,'(10x,i5)') flagext
nbase = nbase
ntx = (nbase + 2)*(nbase +3)
ntx = ntx/2-1
bosc = dsqrt(2.d0*ama/homega)

write(*,*) "********** READER ***************"
write(*,'(a,i5)') " Base size         ",nbase
write(*,'(a,i5)') " Gauss Size         ",ngauss
write(*,'(a,i5)') " Particles Number  ",npart
write(*,'(a,i5)') " Max iteration     ",maxit
write(*,'(a,f10.4)') " Osc length       ",bosc
if (flagbcs .eq. 1) then
write(*,*) "*********** Pairing *************"
write(*,*) "Using BCS with Seniority Pairing"
write(*,'(a,f10.4)') " Pairing strength       ",g_pair
 if (g_pair .le. 0.01) then
   write(*,*) "Pairing strength too low, setting to 0.5"
   g_pair = 0.5d0
   endif
endif
write(*,*) "******** END READER *************"
close(1)
allocate(nr(ntx),nl(ntx),nj(ntx))
end subroutine

!> Just fancy printouts
subroutine fancy(lpr)
implicit none
logical::lpr
if (lpr) then
write(*,*) " _________________________________________________"
write(*,*) "|    _____ _   _    ___ _  _ _____                |"
write(*,*) "|   |_   _/_\ | |  | __| \| |_   _|               |"
write(*,*) "|     | |/ _ \| |__| _|| .` | | |                 |"
write(*,*) "|     |_/_/ \_\____|___|_|\_| |_|                 |" 
write(*,*) "|      __  ________   ____                        |"
write(*,*) "|     / / / / ____/  / __ \______________ ___  __ |"
write(*,*) "|    / /_/ / /_     / / / / ___/ ___/ __ `/ / / / |"
write(*,*) "|   / __  / __/    / /_/ / /  (__  ) /_/ / /_/ /  |"
write(*,*) "|  /_/ /_/_/       \____/_/  /____/\__,_/\__, /   |"
write(*,*) "|                                       /____/    |"
write(*,*) "|_________________________________________________|"
write(*,*) " __________________________________________________"
write(*,*) " Petar Marevic                                      "
write(*,*) "              Maxime Mougeot                       "
write(*,*) "                             Raphaël-David Lasseri  "
write(*,*) " __________________________________________________"
write(*,*) 
write(*,*)
endif
end subroutine






end module constants
