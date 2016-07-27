module constants
implicit none
double precision,parameter::one = 1.d0,two=2.d0,half=0.5d0
!double precision,parameter::pi=3.14159265359
double precision,parameter::pi=acos(-1.d0)
!double precision,parameter::bosc=2.03537781d0 !Oscilator length (fm)
double precision,parameter::hc=197.327d0 ! (MeV*fm)
double precision,parameter::mneutron =939.57d0
!definition of the parameters of the t
double precision,parameter::kr=1.487d0 !in fm -2
double precision,parameter::kt=0.639d0 !in fm -2
double precision,parameter::ks=0.465d0 !in fm -2
double precision,parameter::v0r=200.00d0 !in MeV
double precision,parameter::v0t=178.00d0 !in MeV
double precision,parameter::v0s=91.85d0 !in MeV
double precision,parameter::ama = 20.736209412d0 !in (fm**2)
double precision,parameter::mc2 = 938.90590d0 !in MeV
integer::nbase,npart,maxit,ngauss,n_lines,ntx
integer::na_max,la_max,ja_max,nb_max,lb_max,jb_max
double precision::homega,bosc
integer,allocatable::exttag(:,:,:)
contains
subroutine reader()
implicit none
integer,allocatable:: nr(:),nl(:),nj(:)
open(1,file='hforsay.dat',status='old')
read(1,'(10x,i5)') nbase
read(1,'(10x,i5)') ngauss
read(1,'(10x,i5)') npart
read(1,'(10x,i5)') maxit
read(1,'(10x,f10.4)') homega
nbase = nbase
ntx = (nbase + 2)*(nbase +3)
ntx = ntx/2-1
bosc = dsqrt(2.d0*ama/homega)

write(*,*) "******** READER *************"
write(*,'(a,i5)') " Base size         ",nbase
write(*,'(a,i5)') " Gauss Size         ",ngauss
write(*,'(a,i5)') " Paticules Number  ",npart
write(*,'(a,i5)') " Max iteration     ",maxit
write(*,'(a,f10.4)') " Osc length       ",bosc
write(*,*) "****** END READER ***********"
close(1)
allocate(nr(ntx),nl(ntx),nj(ntx))
end subroutine






end module constants
