module pot  
  contains
  function potential(r1,r2,v0,mu) result(v)
  use constants
  implicit none
  double precision::r1,r2,v,mu,v0
  v=-(one/4.d0*(1-exp(-2*mu*r1*r2))*exp(-mu*(r1**2+r2**2))*v0)/(r1*r2*mu)
  end function
end module

