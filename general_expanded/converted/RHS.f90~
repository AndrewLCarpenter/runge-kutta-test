      subroutine RHS(uvec,res,dt,ep,iprob,sigma,rho,beta)

      implicit none

      integer,   parameter                           :: wp=8

      integer,   parameter                           :: is=9
      integer,   parameter                           :: ivarlen=4

      real(wp),                        intent(in   ) :: sigma
      real(wp),                        intent(in   ) :: rho
      real(wp),                        intent(in   ) :: beta

      real(wp),  dimension(ivarlen),    intent(in   ) :: uvec
      real(wp),  dimension(ivarlen),    intent(  out) :: res
      real(wp),                        intent(in   ) :: dt
      real(wp),                        intent(in   ) :: ep
      integer,                         intent(in   ) :: iprob

      
      if    (iprob.eq.1)then
        res(1) = dt*uvec(2)
        res(2) = dt*((1-uvec(1)*uvec(1))*uvec(2) - uvec(1))/ep
      elseif(iprob.eq.2)then
        res(1) = dt*(-uvec(2))
        res(2) = dt*( uvec(1) + (sin(uvec(1)) - uvec(2))/ep)
      elseif(iprob.eq.3)then
        res(1) = dt*(-(1./ep+2.)*uvec(1) + (1./ep)*uvec(2)*uvec(2))
        res(2) = dt*(uvec(1) - uvec(2) - uvec(2)*uvec(2) )
      elseif(iprob.eq.4)then
        res(1) = dt*(uvec(1) + uvec(2)/ep)
        res(2) = dt*(        - uvec(2))
      elseif(iprob.eq.5)then
        res(1)=dt*sigma*(uvec(2)-uvec(1))/ep
        res(2)=dt*(-uvec(1)*uvec(3)+rho*uvec(1)-uvec(2))
        res(3)=dt*(uvec(1)*uvec(2)-beta*uvec(3))
      endif

      return
      end
