module basis
use constants
integer,allocatable::n_ext(:),m_ext(:),l_ext(:),j_ext(:),t_ext(:)
integer,allocatable::n_red(:),l_red(:),j_red(:),occ(:),nocc(:)
integer,allocatable::repick_base(:),tag_hf(:,:,:),isospin(:)
integer::red_size,size_full,occ_states
integer::ired
integer::lmin,lmax,jmin,jmax
double precision,allocatable::tbme_ext(:,:,:,:)
contains
subroutine external_basis()
implicit none
integer::tag,nn,nl,nm,nj,niso
integer::iho,stat,ifil
integer::n_lines,i,j,count_base
logical::fex
inquire(file='spM.dat', exist=fex)
if (fex) then
open(124,file='spM.dat')
read(124,*) n_lines
size_full = n_lines
allocate(n_ext(n_lines),m_ext(n_lines),l_ext(n_lines),j_ext(n_lines),t_ext(n_lines))
allocate(repick_base(n_lines))
n_ext=0;m_ext=0;l_ext=0;j_ext=0;t_ext=0
j=0
do iho=1,n_lines
!    ifil = iho-1
   read(124,'(I3,2X,5(1X,I3))', iostat=stat) tag,nn,nl,nj,nm,niso
   n_ext(iho) = nn
   l_ext(iho) = nl
   j_ext(iho) = nj
   m_ext(iho) = nm
   t_ext(iho) = niso
   if (niso .eq. 1 .and. nm .eq. nj) then
     j = j + 1
   endif
enddo
close(124)
   red_size = j
   write(*,*) "Reduced Basis Size: ",red_size
   allocate(n_red(red_size),l_red(red_size),j_red(red_size),occ(red_size))
n_red=0;l_red=0;j_red=0;occ=0
repick_base=0
   j= 1
open(124,file='spM.dat')
read(124,*) n_lines
count_base = 0
do iho=1,n_lines 
   read(124,'(I3,2X,5(1X,I3))', iostat=stat) tag,nn,nl,nj,nm,niso
  if (niso ==1 .and. count_base == 0) then
   n_red(j) = nn
   l_red(j) = nl
   j_red(j) = nj
   occ(j) = nj + 1
   count_base = 1
   ired =iho
   endif
      
  if (count_base .eq. 1) then
   if (t_ext(iho) .eq. t_ext(ired) .and. t_ext(ired) .eq. 1) then
    if ((n_ext(ired) .ne. nn) .or. (nl .ne. l_ext(ired)) .or. (j_ext(ired) .ne. nj)) then
     j = j + 1
     n_red(j) = nn
     l_red(j) = nl
     j_red(j) = nj
     !m_red(j) = nm
     occ(j) = nj + 1
     ired=iho
     endif
   endif
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
lmin=minval(l_red);lmax=maxval(l_red);jmin=minval(j_red);jmax=maxval(j_red)
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
integer::m1,m2,m3,m4
integer::l1,l2,l3,l4
integer::j1,j2,j3,j4
double precision::tbme
logical::fex,lpr,tbme_exist
integer::stat
inquire(file='VM-scheme.dat', exist=fex)
inquire(file='tbme.bin', exist=tbme_exist)
if (fex) then
write(*,*) "Reading external TBMEs"
allocate(tbme_ext(red_size,red_size,red_size,red_size))
open(125,file='VM-scheme.dat')
open(unit = 10, file='tbme.bin',form='unformatted')
tbme_ext = 0.d0

!-------- Removing Legend
read(125,*)
read(125,*)
!-------
if (lpr) then
  if (tbme_exist) then
  read(10) tbme_ext
  else
      do
        read(125,*,iostat=stat) n1,n2,n3,n4,tbme
         m1=m_ext(n1);m2=m_ext(n2);m3=m_ext(n3);m4=m_ext(n4)
         j1=j_ext(n1);j2=j_ext(n2);j3=j_ext(n3);j4=j_ext(n4)
         l1=l_ext(n1);l2=l_ext(n2);l3=l_ext(n3);l4=l_ext(n4)
!         m1=m_ext(n1)**2;m2=m_ext(n2)**2;m3=m_ext(n3)**2;m4=m_ext(n4)**2
        !  write(*,*) n1,n2,n3,n4,tbme
        if (t_ext(n1) .eq. 1 .and. t_ext(n2) .eq. 1 .and.t_ext(n3) .eq. 1 .and.t_ext(n4) .eq. 1) then !Keeping only neutrons elements
           if (tbme .ne. 0.d0) then 
            if (m1 .eq. m3 .and. m2 .eq. m4) then
             if (j1 .eq. j3 .and. j2 .eq. j4) then
              if (l1 .eq. l3 .and. l2 .eq. l4) then
            !Keeping only one projection of J
            q1 = repick_base(n1)
            q2 = repick_base(n2)
            q3 = repick_base(n3)
            q4 = repick_base(n4)
         !   write(1234,*) n1,n2,n3,n4
   !         write(*,*) q1,q2,q3,q4,tbme
            tbme_ext(q1,q2,q3,q4) = tbme_ext(q1,q2,q3,q4) + 1.d0/(dble((1+j_ext(n1))*(1+j_ext(n2))))*tbme
            !     tbme_ext(q1,q2,q3,q4) = tbme
        !    if (tbme .ne. 0.d0) write(127,*) q1,q2,q3,q4,tbme_ext(q1,q2,q3,q4)
          endif ! Filtering over m
        endif !Filtering Isospin
        endif !Filtering Isospin
        endif !Filtering Isospin
        endif !Filtering Isospin
        if (stat /= 0) then
          write(*,*) "End of VM-Scheme.dat"
          exit
      endif
  enddo ! End of the tbme file
  write(10) tbme_ext
  close(10)
endif
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
