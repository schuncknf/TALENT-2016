module pot
  contains

  function potential(r1,r2,v0,mu) result(v)
  use constants
  implicit none
  double precision::r1,r2,v,mu,v0
  v=(one*v0*(exp(2*mu*r1*r2)-exp(-2*mu*r1*r2))*exp(-mu*(r1**2+r2**2)))/(4.d0*r1*r2*mu)
  end function

  function minnesota(r1,r2,vc,cc) result(res)
  use constants
  implicit none
  double precision::r1,r2,vc,cc,res
  res = vc/cc*half*exp(-cc*(r1**2+r2**2))*sinh(2*cc*r1*r2)
  end function

end module
