      subroutine Jacobian(uvec,xjac,dt,ep,akk,iprob,nvecLen,sigma,rho,beta)

      implicit none

      integer,   parameter                           :: wp=8

      integer,   parameter                           :: ivarlen=4
 
      real(wp),                                intent(in   )  :: sigma
      real(wp),                                intent(in   )  :: rho
      real(wp),                                intent(in   )  :: beta

      real(wp),   dimension(ivarlen),          intent(in   )  :: uvec
      real(wp),   dimension(ivarlen,ivarlen),  intent(  out)  :: xjac
      real(wp),                                intent(in   )  :: dt
      real(wp),                                intent(in   )  :: ep
      real(wp) ,                               intent(in   )  :: akk
      integer,                                 intent(in   )  :: iprob
      integer,                                 intent(in   )  :: nveclen

!     U_t = F(U)  ;  Jac = \frac{\partial F(U)}{\partial U}

!     Note:  
!     Jac is NOT complete LHS Implicit matrix.  
!     It needs the time term and proper scaling

!     if(iprob.eq.1)then
!        jac(1,1) = 0.0_wp
!        jac(1,2) = 1.0_wp
!        jac(2,1) = (-2.0_wp*uvec(1)*uvec(2)-1.0_wp)/ep
!        jac(2,2) = ( 1.0_wp-uvec(1)*uvec(1)       )/ep

!        ! Note: Storing the (1,1) zero 
!        ! final matrix matrix includes contribution from Identity matrix
!        ia(1) = 1 ; ia(2) = 3 ; ia(3) = 5 ;
!        ja(1) = 1        ; ja(2) = 2        ; 
!        ja(3) = 1        ; ja(4) = 2        ;
!         a(1) = jac(1,1) ;  a(2) = jac(1,2) ;  
!         a(3) = jac(2,1) ;  a(4) = jac(2,2) ;  
!     endif

      if(iprob.eq.1)then
        xjac(1,1) = 1.-akk*dt*(0.)
        xjac(1,2) = 0.-akk*dt*(1.)
        xjac(2,1) = 0.-akk*dt*(-2*uvec(1)*uvec(2)-1)/ep
        xjac(2,2) = 1.-akk*dt*(1-uvec(1)*uvec(1))/ep
      elseif(iprob.eq.2)then
        xjac(1,1) = 1.-akk*dt*(0.)
        xjac(1,2) = 0.-akk*dt*(0.)
        xjac(2,1) = 0.-akk*dt*(cos(uvec(1)) )/ep
        xjac(2,2) = 1.-akk*dt*(-1.)/ep
      elseif(iprob.eq.3)then
        xjac(1,1) = 1.-akk*dt*(-(1./ep))
        xjac(1,2) = 0.-akk*dt*(+1./ep)*2*uvec(2)
        xjac(2,1) = 0.-akk*dt*(0)
        xjac(2,2) = 1.-akk*dt*(0)
      elseif(iprob.eq.4)then
        xjac(1,1) = 1.-akk*dt
        xjac(1,2) = 0.-akk*dt/ep
        xjac(2,1) = 0.-akk*dt*(0)
        xjac(2,2) = 1.+akk*dt/ep !possible error here
       elseif(iprob.eq.5)then
        xjac(1,1) = 1.-akk*dt*(-sigma)/ep
        xjac(1,2) = 0.-akk*dt*(sigma)/ep
        xjac(1,3) = 0.-akk*dt*(0)/ep

        xjac(2,1) = 0.-akk*dt*(rho-uvec(3))
        xjac(2,2) = 1.-akk*dt*(-1.)
        xjac(2,3) = 0.-akk*dt*(-uvec(1))

        xjac(3,1) = 0.-akk*dt*(uvec(2))
        xjac(3,2) = 0.-akk*dt*(uvec(1))
        xjac(3,3) = 1.-akk*dt*(-beta)
      endif

      return

      end
