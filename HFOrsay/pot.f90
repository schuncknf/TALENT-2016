module pot
  contains

  function potential(r1,r2,v0,mu) result(v)
  use constants
  implicit none
  double precision::r1,r2,v,mu,v0
  r1=dsqrt(r1)*bosc
  r2=dsqrt(r2)*bosc
  v=(v0*(exp(-mu*(r1**2+r2**2))*sinh(2*mu*r1*r2)/(2.d0*r1*r2*mu)))
  end function

  function minnesota(r1,r2) result(res)
  use constants
  implicit none
  double precision::r1,r2,a,b,res
!  res = one/(r1*r2)*vc/cc*half*exp(-cc*(r1**2+r2**2))*sinh(2*cc*r1*r2)
  a = half*v0r/kr*(exp(-kr*(r1**2-2*r1*r2+r2**2))-exp(-kr*(r1**2+r2**2)))
  b = -half*v0s/ks*(exp(-ks*(r1**2-2*r1*r2+r2**2))-exp(-ks*(r1**2+r2**2)))
  res = a+b/(4.d0*r1*r2)
  end function

end module
