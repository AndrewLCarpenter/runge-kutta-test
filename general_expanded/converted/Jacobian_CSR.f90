      subroutine Jacobian_CSR(uvec,dt,ep,akk,iprob,nvecLen,sigma,rho,beta)

      use precision_vars
      use CSR_Variables

      implicit none
 
      real(wp),                         intent(in   )  :: sigma
      real(wp),                         intent(in   )  :: rho
      real(wp),                         intent(in   )  :: beta

      real(wp),   dimension(nJac),      intent(in   )  :: uvec
      real(wp),                         intent(in   )  :: dt
      real(wp),                         intent(in   )  :: ep
      real(wp) ,                        intent(in   )  :: akk
      integer,                          intent(in   )  :: iprob
      integer,                          intent(in   )  :: nveclen

      real(wp),   dimension(nJac,nJac)                 :: xjac
      integer                                          :: i,j
      integer                                          :: icnt, jcnt

!     U_t = F(U)  ;  Jac = \frac{\partial F(U)}{\partial U}  ;  xjac = I - akk dt Jac

      if(iprob.eq.1)then
        xjac(1,1) = 1.0_wp-akk*dt*(0.0_wp)
        xjac(1,2) = 0.0_wp-akk*dt*(1.0_wp)
        xjac(2,1) = 0.0_wp-akk*dt*(-2*uvec(1)*uvec(2)-1)/ep
        xjac(2,2) = 1.0_wp-akk*dt*(+1-uvec(1)*uvec(1)  )/ep
      elseif(iprob.eq.2)then
        xjac(1,1) = 1.0_wp-akk*dt*(0.0_wp)
        xjac(1,2) = 0.0_wp-akk*dt*(0.0_wp)
        xjac(2,1) = 0.0_wp-akk*dt*(cos(uvec(1)) )/ep
        xjac(2,2) = 1.0_wp-akk*dt*(-1.0_wp)/ep
      elseif(iprob.eq.3)then
        xjac(1,1) = 1.0_wp-akk*dt*(-(1.0_wp/ep))
        xjac(1,2) = 0.0_wp-akk*dt*(+1.0_wp/ep)*2*uvec(2)
        xjac(2,1) = 0.0_wp-akk*dt*(0)
        xjac(2,2) = 1.0_wp-akk*dt*(0)
      elseif(iprob.eq.4)then
        xjac(1,1) = 1.0_wp-akk*dt
        xjac(1,2) = 0.0_wp-akk*dt/ep
        xjac(2,1) = 0.0_wp-akk*dt*(0.0_wp)
        xjac(2,2) = 1.0_wp+akk*dt/ep !possible error here
       elseif(iprob.eq.5)then
        xjac(1,1) = 1.0_wp-akk*dt*(-sigma)/ep
        xjac(1,2) = 0.0_wp-akk*dt*(sigma)/ep
        xjac(1,3) = 0.0_wp-akk*dt*(0.0_wp)/ep

        xjac(2,1) = 0.0_wp-akk*dt*(rho-uvec(3))
        xjac(2,2) = 1.0_wp-akk*dt*(-1.0_wp)
        xjac(2,3) = 0.0_wp-akk*dt*(-uvec(1))

        xjac(3,1) = 0.0_wp-akk*dt*(uvec(2))
        xjac(3,2) = 0.0_wp-akk*dt*(uvec(1))
        xjac(3,3) = 1.0_wp-akk*dt*(-beta)
      endif

      ! Initialize CSR
      iaJac(:) = 0  ;  iaJac(1) = 1
      jaJac(:) = 0  ;   aJac(:) = 0.0_wp

      ! Store dense matrix into CSR format
      icnt = 0   ;
      do i = 1,nJac
        jcnt = 0   ;
        do j = 1,nJac
          if(abs(xjac(i,j)) >= tolJac) then
            icnt = icnt + 1 ; jcnt = jcnt + 1 ; 
            jaJac(icnt) = j
             aJac(icnt) = xjac(i,j)
          endif
        enddo
        iaJac(i+1) = iaJac(i) + jcnt
      enddo

      return

      end subroutine Jacobian_CSR
