module constants
implicit none
double precision,parameter::one = 1.d0,two=2.d0,half=0.5d0
double precision,parameter::pi=3.14159265359
double precision,parameter::bosc=10.d0 !Oscilator length

!definition of the parameters of the t
double precision,parameter::kr=1.487d0 !in fm -2
double precision,parameter::kt=0.639d0 !in fm -2
double precision,parameter::ks=0.465d0 !in fm -2
double precision,parameter::v0r=200.00d0 !in MeV
double precision,parameter::v0t=178.00d0 !in MeV
double precision,parameter::v0s=91.85d0 !in MeV

double precision,parameter::ama = 20.73d0 !in fm 2
double precision,parameter::mc2 = 938.90590d0 !in MeV
integer::nbase,npart,maxit
contains
subroutine reader()
implicit none
open(1,file='hforsay.dat',status='old')
read(1,'(10x,i5)') nbase
read(1,'(10x,i5)') npart
read(1,'(12x,i5)') maxit

write(*,*) "******** READER *************"
write(*,'(a,i5)') " Base size         ",nbase
write(*,'(a,i5)') " Paticules Number  ",npart
write(*,'(a,i5)') " Max iteration     ",maxit
write(*,*) "****** END READER ***********"



end subroutine



end module constants
