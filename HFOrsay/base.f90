subroutine sphbasis(n,nr,nl,nj,nocc,lpr)
      use constants
      implicit none
      integer:: il,nnsph,nlsph,nrsph,mssph,njsph,n,nt 
      logical:: lpr
      integer::nr(n),nl(n),nj(n),ms(n),nocc(n)
      integer::mspin,noc

      nr(1) = 0
      nl(1) = 0
      nj(1) = 1
      ms(1) = 1
      nocc(1) = 2
!
      il = 1
!
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
                 nocc(il)= 2*(nlsph+mssph)
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
	            nocc(il) = 2*(nlsph+mssph)
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
	noc = nocc(il)

        !write(*,110)'NN = ',nnsph,'nr = ',nrsph,'ml = ',nlsph,'(2*nj-1)/2 = ',2*njsph-1,'/2'
        write(*,110)'NN = ',nnsph,'nr = ',nrsph,'nl = ',nlsph,'ms = ', mspin, '/2', 'nocc = ', noc

  110   format(5x,a,i2,3x,a,i2,3x,a,i2,3x,a,i2,a,3x,a,i2)

        enddo

        write(*,*) '****** END SPHERICAL BASIS ***********'

        endif
end subroutine sphbasis

!subroutine external_basis()
!implicit none
!integer::na,nb,la,lb,ja,jb,a,b
!integer::i
!integer::n_lines
!allocate(resua(0:na_max-1,0:la_max-1,0:ja_max-1))
!allocate(resub(0:nb_max-1,0:lb_max-1,0:jb_max-1))
!open(124,file='base_config')
!resua=0
!resub=0
!do i=1,n_lines
!read(124,'(8i4)') na,la,ja,a,nb,lb,jb,b
!resua(na,la,ja) = a
!resub(nb,lb,jb) = b
!enddo
!end subroutine




