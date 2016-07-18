program matrix_diag


implicit none
double precision, dimension(:,:), allocatable :: M,VL,VR
integer::n,flag,i,j,k,l,lwork
integer,parameter:: seed = 7626
INTEGER::LDA, LDVL, LDVR,info
double precision::r
INTEGER,parameter::LWMAX=1000
double precision,dimension(:),allocatable::WR,WI,WORK

external dgeev
EXTERNAL         PRINT_EIGENVALUES, PRINT_EIGENVECTORS

!Input
write(*,*) "Random: 2, Unity: 1"
read(*,*) flag
write(*,*) "Dimension of the Matrix ?"
read(*,*) n
lda=n
ldvl=n
ldvr=n


allocate(WR(n))
allocate(WI(n))
allocate(WORK(lwmax))
allocate(VL(n,n))
allocate(VR(n,n))


!Filling the matrix
allocate(M(n,n))
M(:,:) = 0.d0
if (flag .eq. 1) then
   do i=1,n
    do j=1,n
       if (i == j) M(i,j) = 1.d0
    enddo
    enddo
endif
if (flag .eq. 2) then
    call srand(seed)
    do i=1,n
     do j = 1,n
     call random_number(r)
     M(i,j) = r
    enddo
   enddo
endif
    
!Diagonalisation
lwork=-1
call dgeev('Vectors', 'Vectors', N, M,N, WR, WI, VL, LDVL,VR, LDVR, WORK, LWORK, INFO )
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
call dgeev('Vectors', 'Vectors', N, M,N, WR, WI, VL, LDVL,VR, LDVR, WORK, LWORK, INFO )

if (info .gt. 0) then 
write(*,*) "Failing Miserably..."
stop
endif
! Printout of the Matrix
write(*,*) "Matrix values"
do i=1,n
  do j=1,n
      write (*,*) m(i,j)
  enddo
  write(*,*) "    "
enddo
! Printout Eigenvalues
CALL PRINT_EIGENVALUES('Eigenvalues', N, WR, WI )

end
