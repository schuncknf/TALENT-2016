module basis

use constants

implicit none 

integer:: nr(ntx),nl(ntx),nj(ntx),nt
logical:: lpr

contains 

subroutine sphbasis(lpr)
      integer:: il,nnsph,nlsph,nrsph,mssph 

      il = 1
      do nnsph = 1, n0base
 
         if(mod(nnsph,2) .eq. 0) then
            do nlsph = 0,nnsph,2
	    nrsph = (nnsph-nlsph)/2	
            do mssph=0,1
	       if(nlsph+mssph .gt. 0) then
	         il=il+1
	         nr(il)=nrsph
	         nl(il)=nlsph
	         nj(il)=nlsph+mssph
	       endif
            enddo
         enddo

        elseif(mod(nnsph,2) .eq. 1) then
           do nlsph = 1,nnsph,2
              nrsph = (nnsph-nlsph)/2
	      do mssph = 0,1
                 if(nlsph+mssph .gt. 0) then
	            il = il+1
		    nr(il) = nrsph
		    nl(il) = nlsph
	            nj(il) = nlsph+mssph
	        endif
	     enddo
	  enddo
c
	else
	  stop 'error in base'
	endif	

        if(nnsph .eq. n0base) nt = il

        nr(1) = 0
        nl(1) = 0
        nj(1) = 1

! printout
        if(lpr) then

        write(*,*) ' ****** SPHERICAL BASIS ***************'
   
        do il=1,nt

	nrsph = nr(il)
	nlsph = nl(il)
	njsph = nj(il)
	nnsph = 2*nrsph + nlsph

        write(*,110)'NN = ',nnsph,'nr = ',nrsph,'ml = ',
     &  nlsph,'(2*nj-1)/2 = ',2*njsph-1,'/2'

  110   format(5x,a,i2,3x,a,i2,3x,a,i2,3x,a,i2,a)

        enddo

        write(l6,*) '****** END SPHERICAL BASIS ***********'

        endif
	       
end subroutine sphbasis
	
end module basis



