      subroutine Kreiss(programStep,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE,resI,akk,xjac)
      implicit none

      integer,  parameter                      :: wp=8
      integer,  parameter                      :: vecl=4

      integer,                   intent(in   ) :: programStep

      !INIT vars
      real(wp), dimension(vecl), intent(inout) :: uvec
      real(wp),                  intent(in   ) :: ep
      real(wp), dimension(vecl), intent(  out) :: uexact
      real(wp),                  intent(inout) :: dt
      integer,                   intent(inout) :: nveclen
      real(wp),                  intent(  out) :: tfinal
      integer,                   intent(in   ) :: iDT

      real(wp)                                 :: tmp,tmpM,tmpP

      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE,resI
      
      !Jacob vars
      real(wp),                       intent(in   ) :: akk
      real(wp), dimension(vecl,vecl), intent(  out) :: xjac


      if (programStep==0) then
        dt = 0.25_wp/10**((iDT-1)/20.0_wp)
        nvecLen = 2

        tfinal = 1.0_wp
        tmpP = exp(+tfinal   )
        tmpM = exp(-tfinal/ep)
        tmp = exp(tfinal)
        uvec(1) = 1.0_wp
        uvec(2) = 1.0_wp
        uexact(1) = uvec(1)*tmpP + uvec(2)*(tmpP - tmpM)/(1.0_wp + ep)
        uexact(2) = uvec(2)*tmpM

      elseif (programStep==1 .or.programStep==2) then
        resE(1) = dt*uvec(1)
        resE(2) = dt*(-uvec(2))
        resI(1) = dt*(uvec(2)/ep)
        resI(2) = 0.0_wp
      elseif (programStep==3) then
        xjac(1,1) = 1.-akk*dt*(0.0_wp)
        xjac(1,2) = 0.-akk*dt/ep
        xjac(2,1) = 0.-akk*dt*( 0.0_wp)
        xjac(2,2) = 1.-akk*dt*(0.0_wp) !possible error here
      endif

      return
      end subroutine

      

