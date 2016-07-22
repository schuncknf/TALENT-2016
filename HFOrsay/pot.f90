module pot
  contains
  function potential(r1,r2,v0,mu) result(v)
  use constants
  implicit none
  double precision::r1,r2,v,mu,v0
  r1=dsqrt(r1)*bosc
  r2=dsqrt(r2)*bosc
  v=(v0*(exp(-mu*(r1**2+r2**2))*sinh(2.d0*mu*r1*r2)/(4.d0*r1*r2*mu)))
  end function
end module
