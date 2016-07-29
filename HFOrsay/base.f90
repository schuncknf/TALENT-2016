subroutine sphbasis(n,nr,nl,nj,nocc,lpr)
      use constants
      implicit none
      integer:: il,nnsph,nlsph,nrsph,mssph,njsph,n,nt
      logical:: lpr
      integer::nr(n),nl(n),nj(n),ms(n),nfull(n),nocc(n) 
      integer::mspin,noc
      integer::i,nf,ic
! nfull - maximal number of particles in specific state
! nocc  - true number of particles in specific state
      nr(1) = 0
      nl(1) = 0
      nj(1) = 1
      ms(1) = 1
      nfull(1) = 2
      il = 1
      do nnsph = 1, nbase
         if(mod(nnsph,2) .eq. 0) then
            do nlsph = 0,nnsph,2
            nrsph = (nnsph-nlsph)/2
            do mssph=0,1
               if(nlsph+mssph .gt. 0) then
                 il=il+1
                 nr(il)=nrsph
                 nl(il)=nlsph
                 nj(il)=nlsph+mssph
	         ms(il)=mssph	
                 nfull(il)= 2*(nlsph+mssph)
               endif
            enddo ! mssph
         enddo !nlsph

        elseif(mod(nnsph,2) .eq. 1) then
           do nlsph = 1,nnsph,2
              nrsph = (nnsph-nlsph)/2
              do mssph = 0,1
                 if(nlsph+mssph .gt. 0) then
                    il = il+1
                    nr(il) = nrsph
                    nl(il) = nlsph
                    nj(il) = nlsph+mssph
	            ms(il)=mssph
	            nfull(il) = 2*(nlsph+mssph) ! possible number of particles in a state
                endif
             enddo !nlsph
          enddo !mssph
        else
          stop 'error in base'
        endif

        if(nnsph .eq. nbase) nt = il

        enddo

! printout
        if(lpr) then

        write(*,*) ' ****** SPHERICAL BASIS ***************'
   
        do il=1,nt
        nrsph = nr(il)
        nlsph = nl(il)
        njsph = nj(il)
        nnsph = 2*nrsph + nlsph
	mspin = 2*(njsph-nlsph)-1
	noc = nfull(il)

        !write(*,110)'NN = ',nnsph,'nr = ',nrsph,'ml = ',nlsph,'(2*nj-1)/2 = ',2*njsph-1,'/2'
        write(*,110)'NN = ',nnsph,'nr = ',nrsph,'nl = ',nlsph,'ms = ', mspin, '/2', 'npart = ', noc

  110   format(5x,a,i2,3x,a,i2,3x,a,i2,3x,a,i2,a,3x,a,i2)

        enddo

        write(*,*) '****** END SPHERICAL BASIS ***********'

        endif


! determinaton of occupation number of each state
	
	nf = 0
	ic = 1
        do il = 1, nt
	   nf = nf + nfull(il)
	   if (nf .le. npart) then
	      nocc(il)=nfull(il)
	   elseif(nf .gt. npart .and. ic .eq. 1) then 
	      nocc(il)=nfull(il)+npart-nf
	      ic = 2
	   else 
	      nocc(il)=0
	   endif
        enddo

	if(lpr) then
	do il = 1, nt
	  write(*,*) il, nocc(il),nfull(il)
        enddo
	endif

!
end subroutine sphbasis

module basis
use constants
integer,allocatable::n_ext(:),m_ext(:),l_ext(:),j_ext(:),t_ext(:)
integer,allocatable::n_red(:),l_red(:),j_red(:),occ(:),nocc(:)
integer,allocatable::repick_base(:),tag_hf(:,:,:)
integer::red_size,size_full,occ_states
double precision,allocatable::tbme_ext(:,:,:,:)
contains 
subroutine external_basis()
implicit none
integer::tag,nn,nl,nm,nj,niso
integer::iho,stat,ifil
integer::n_lines,i,j
logical::fex
inquire(file='spM.dat', exist=fex)
if (fex) then
open(124,file='spM.dat')
read(124,*) n_lines
size_full = n_lines
allocate(n_ext(n_lines),m_ext(n_lines),l_ext(n_lines),j_ext(n_lines),t_ext(n_lines))
allocate(repick_base(n_lines))
n_ext=0;m_ext=0;l_ext=0;j_ext=0;t_ext=0
n_red=0;l_red=0;j_red=0;occ=0
repick_base=0
j=0
do iho=1,n_lines
!    ifil = iho-1
   read(124,'(I3,2X,5(1X,I3))', iostat=stat) tag,nn,nl,nj,nm,niso 
   n_ext(iho) = nn
   l_ext(iho) = nl
   j_ext(iho) = nj
   m_ext(iho) = nm
   t_ext(iho) = niso
   if (niso == 1 .and. nm .eq. nj) then    
     j = j + 1
   endif
enddo
close(124)
   red_size = j
   write(*,*) "Reduced Basis Size: ",red_size
   allocate(n_red(red_size),l_red(red_size),j_red(red_size),occ(red_size))
   j= 0
open(124,file='spM.dat')
read(124,*) n_lines
do iho=1,n_lines
   read(124,'(I3,2X,5(1X,I3))', iostat=stat) tag,nn,nl,nj,nm,niso 
   if (niso == 1 .and. nm .eq. nj) then    
     j = j + 1
     n_red(j) = nn
     l_red(j) = nl 
     j_red(j) = nj
     occ(j) = nj + 1
   endif
     repick_base(iho) = j
 !    write(*,*) "End Basis"
   if (stat /= 0) then
   write(*,*) "Problem while reading spM.dat"
   exit
   endif
enddo
allocate(tag_hf(minval(n_red):maxval(n_red),minval(l_red):maxval(l_red),minval(j_red):maxval(j_red)))
tag_hf = 0
do j=1,red_size
 tag_hf(n_red(j),l_red(j),j_red(j))=j
enddo
close(124)
else
write(*,*) "File spM.dat not found !"
stop
endif
end subroutine

subroutine external_tbme(lpr)
implicit none
integer::q1,q2,q3,q4
integer::n1,n2,n3,n4
double precision::tbme
logical::fex,lpr
integer::stat
inquire(file='VM-scheme.dat', exist=fex)
if (fex) then
write(*,*) "Reading external TBMEs"
allocate(tbme_ext(red_size,red_size,red_size,red_size))
open(125,file='VM-scheme.dat')
tbme_ext = 0.d0
!-------- Removing Legend
read(125,*)
read(125,*)
!-------
if (lpr) then
do
  read(125,*,iostat=stat) n1,n2,n3,n4,tbme
!  write(*,*) n1,n2,n3,n4,tbme
  if (t_ext(n1) .eq. 1 .and. t_ext(n2) .eq. 1 .and.t_ext(n3) .eq. 1 .and.t_ext(n4) .eq. 1) then !Keeping only neutrons elements
   if (m_ext(n1) .eq. j_ext(n1) .and. m_ext(n2) .eq. j_ext(n2) .and. m_ext(n3) .eq. j_ext(n3) .and. m_ext(n4) .eq. j_ext(n4)) then
   !Keeping only one projection of J
        q1 = repick_base(n1)
        q2 = repick_base(n2)
        q3 = repick_base(n3)
        q4 = repick_base(n4)
        write(12345,*) q1,q2,q3,q4,tbme
        tbme_ext(q1,q2,q3,q4) = 1.d0/(dble((1+j_ext(n1))*(1+j_ext(n2))))*tbme 
   !     tbme_ext(q1,q2,q3,q4) = tbme 
        if (tbme .ne. 0.d0) write(127,*) q1,q2,q3,q4,tbme_ext(q1,q2,q3,q4)
    endif ! Filtering over m
   endif !Filtering Isospin
   if (stat /= 0) then
   write(*,*) "End of VM-Scheme.dat"
   exit
   endif
enddo ! End of the tbme file 
endif

else
write(*,*) "File VM-Scheme.dat not found !"
stop
endif




end subroutine


subroutine filled_number()
implicit none
integer::il,nf,ic
allocate(nocc(red_size))
nf = 0
nocc=0
ic = 1
occ_states = 0
 do il = 1, red_size
   nf = nf + occ(il)
   if (nf .le. npart) then
      nocc(il)=occ(il)
   elseif(nf .gt. npart .and. ic .eq. 1) then
      nocc(il)=occ(il)+npart-nf
      ic = 2
   else
      nocc(il)=0
   endif
 enddo
 do il = 1, red_size
  if (nocc(il) .ne. 0) occ_states = occ_states +1
enddo

end subroutine


end

       

