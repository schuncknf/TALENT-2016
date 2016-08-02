module maths
contains
function fac(n) result(fak)
implicit none
integer::i,n
double precision::lfak,fak
if (n==0) then
 fak=1.d0
else
 lfak = 0.d0
 do i=1,n
 lfak = lfak + log(dble(i))
 enddo
 fak = exp(lfak)
endif
end function fac
function ffac(n) result (ffak)
implicit none
integer::n,k
double precision::ffak
if (mod(n,2) == 0) then
  k=n/2
  ffak = 2.d0**(dble(k))*fac(k)
 else
  k=(n+1)/2
  ffak = fac(2*k)/(2**(dble(k))*fac(k))
endif
end function ffac
integer function delta(i,j)
  integer :: i, j
  delta=0; if(i==j) delta=1
end function delta
 subroutine sort(n,a)
  implicit none
  integer n,i,j
        double precision a(n),x
        do 30 i=2,n
        x=a(i)
        j=i
   10   j=j-1
        if(j.eq.0 .or. a(j).le.x) go to 20
        a(j+1)=a(j)
        go to 10
   20   a(j+1)=x
   30   continue
        end subroutine

subroutine sorteigv(n,indexarray,arraytosort)
implicit none
integer::n,i,j,k,l
double precision::val1,val2
double precision::indexarray(n),arraytosort(n),tempar(n),te(n)
integer::sorted(n)

tempar=indexarray
sorted =0
call sort(n,tempar)
do i=1,n
 val1=indexarray(i)
 do j=1,n
         val2=tempar(j)
 if (val1 .eq. val2) sorted(j)=i
 enddo
enddo
do i=1,n
  te(i) = arraytosort(sorted(i))
enddo
arraytosort= 0.d0
arraytosort=te
call sort(n,indexarray)
end subroutine

function trace(M,n) result(tr)
implicit none
integer::n,i
double precision::M(n,n),tr
tr = 0.d0
do i = 1,n
tr = tr + M(i,i)
enddo
end function

!function brentzeros(x1,x2,y1,y2,tol) result(rtbrent)
!use constants
!!     Using Brent's method find the root of the function FUNC(X)
!!     known to lie between X1 and X2.
!!     Y1 = FUNC(X1), Y2 = FUNC(X2)
!!     The root will be returned as RTBRENT with accuracy TOL
!      implicit none
!      double precision::a,b,c,d,fa,fb,fc,e,tol1,r
!      double precision::x1,x2,y1,y2,tol,xm,s,p,q
!      double precision::rtbrent
!      integer,parameter:: itmax = 100
!      integer::iter
!      double precision,parameter::eps = 1.d-12
!      a  = x1
!      b  = x2
!      c = a
!      fa = y1
!      fb = y2
!      if (fa*fb.gt.zero) stop ' in RTBRENT: root must be bracketed'
!      fc = fb
!      do 10 iter = 1,itmax
!         if (fb*fc.gt.zero) then
!            c  = a
!            fc = fa
!            d  = b - a
!            e  = d
!         endif
!         if (abs(fc).lt.abs(fb)) then
!            a  = b
!            b  = c
!            c  = a
!            fa = fb
!           fb = fc
!            fc = fa
!         endif
!         tol1 = tol*abs(b)
!         xm = half*(c-b)
!         if (abs(xm).le.tol1 .or. fb.eq.zero) then
!            rtbrent = b
!            return
!         endif
!        if (abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
!            s = fb/fa
!            if (a.eq.c) then
!               p = 2*xm*s
!               q = one - s
!            else
!               q = fa/fc
!               r = fb/fc
!               p = s*(2*xm*q*(q-r) - (b-a)*(r-one))
!               q = (q-one)*(r-one)*(s-one)
!            endif
!            if (p.gt.zero) q = -q
!            p = abs(p)
!            if (2*p.lt.min(3*xm*q-abs(tol1*q),abs(e*q))) then
!               e = d
!               d = p/q
!            else
!               d = xm
!               e = d
!            endif
!         else
!            d = xm
!            e = d
!                    endif
!         else
!            d = xm
!            e = d
!         endif
!         a  = b
!         fa = fb
!         if (abs(d).gt.tol1) then
!            b = b + d
!         else
!            b = b + sign(tol1,xm)
!         endif
!         fb = func(b)
!   10 continue
!      stop ' in RTBRENT: exceeding maximum number of iterations'
!      end
!
!

end module maths
