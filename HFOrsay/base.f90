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

integer,allocatable::n_ext(:),m_ext(:),l_ext(:),j_ext(:),t_ext(:)

contains 
subroutine external_basis()
implicit none
integer::tag,nn,nl,nm,nj,niso
integer::i,stat
integer::n_lines
open(124,file='spM.dat')
read(124,*) n_lines
allocate(n_ext(n_lines),m_ext(n_lines),l_ext(n_lines),j_ext(n_lines),t_ext(n_lines))
n_ext=0;m_ext=0;l_ext=0;l_ext=0;j_ext=0;t_ext=0
do i=1,n_lines
   read(124,'(I3,2X,5(1X,I3))', iostat=stat) tag,nn,nl,nm,nj,niso 
   n_ext(i) = nn
   l_ext(i) = nl
   m_ext(i) = nm
   j_ext(i) = nj
   t_ext(i) = niso
   if (stat /= 0) then
   write(*,*) "Problem while reading spM.dat"
   exit
   endif
enddo
close(124)
end subroutine
end

       

