!> There's something Doxygen doesn't like about this particular subroutine.
  subroutine DME_fields(r,rho,tau,del_rho)
    implicit none
    real(dp), intent(in) :: r
    real(dp), intent(out) :: rho,tau,del_rho
    real(dp), dimension(0:2,0:5,0:2) :: Rnl
    real(dp) :: xi
    integer :: i,ji,j,nj,lj,jj,k,nk,lk,jk
    xi = r/b_ho
    call RadialHOALL(5,2,xi,b_ho,Rnl)
    rho = 0
    tau = 0
    del_rho = 0
    do i = 1,Noccupied
       ji = j_hf(i)
       do j = 1,Nsize
          nj = n_hf(j)
          lj = l_hf(j)
          jj = j_hf(j)
          do k = 1,Nsize
             nk = n_hf(k)
             lk = l_hf(k)
             jk = j_hf(k)
             if(lk.ne.lj.or.jj.ne.jk) cycle
             rho = rho + (ji+1)*D_mat( j,i)*D_mat( k,i)*&
                  Rnl(0,nj,lj)*Rnl(0,nk,lk)
             tau = tau + (ji+1)*D_mat( j,i)*D_mat( k,i)*&
                  (Rnl(1,nj,lj)*Rnl(1,nk,lk)+&
                  lj*(lj+1)/r**2*Rnl(0,nj,lj)*Rnl(0,nk,lk))
             del_rho = del_rho + (ji+1)*D_mat( j,i)*D_mat( k,i)*&
                  (2*(Rnl(0,nj,lj)*Rnl(1,nk,lk)/r + &
                      Rnl(0,nk,lk)*Rnl(1,nj,lj)/r + &
                      Rnl(1,nj,lj)*Rnl(1,nk,lk)) + & 
                   Rnl(0,nj,lj)*Rnl(2,nk,lk) + Rnl(0,nk,lk)*Rnl(2,nj,lj))

          enddo
       enddo
    enddo
    rho = rho/(4*pi)
    tau = tau/(4*pi)
    del_rho = del_rho/(4*pi)
 subroutine DME_fields
