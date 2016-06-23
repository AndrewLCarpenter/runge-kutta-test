      subroutine Kaps(programStep,probname,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE,resI,akk,xjac)

      use precision_vars

      implicit none

      integer,  parameter                      :: vecl=4

      integer,                   intent(in   ) :: programStep

      !INIT vars
      character(len=9),          intent(  out) :: probname
      real(wp), dimension(vecl), intent(inout) :: uvec
      real(wp),                  intent(in   ) :: ep
      real(wp), dimension(vecl), intent(  out) :: uexact
      real(wp),                  intent(inout) :: dt
      integer,                   intent(inout) :: nveclen
      real(wp),                  intent(  out) :: tfinal
      integer,                   intent(in   ) :: iDT

      real(wp)                                 :: tmp

      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE,resI
      
      !Jacob vars
      real(wp),                       intent(in   ) :: akk
      real(wp), dimension(vecl,vecl), intent(  out) :: xjac


      if (programStep==0) then
        probname='Kaps     '
        dt = 0.5_wp/10**((iDT-1)/20.0_wp)
        nvecLen = 2

        tfinal = 1.0_wp
        tmp = exp(-tfinal)
        uexact(1) = tmp*tmp
        uexact(2) = tmp

!  IC: problem 1   :  equilibrium IC
        uvec(1) = 1.0_wp
        uvec(2) = 1.0_wp
      elseif (programStep==1 .or.programStep==2) then
        resE(1) = dt*(-2.0_wp*uvec(1))
        resE(2) = dt*(uvec(1) - uvec(2) - uvec(2)*uvec(2) )
        resI(1) = dt*(-1./ep*uvec(1) + (1./ep)*uvec(2)*uvec(2))
        resI(2) = 0.0_wp

      elseif (programStep==3) then
        xjac(1,1) = 1.-akk*dt*(-(1./ep))
        xjac(1,2) = 0.-akk*dt*(+1./ep)*2*uvec(2)
        xjac(2,1) = 0.-akk*dt*(0)
        xjac(2,2) = 1.-akk*dt*(0)
      endif

      return
      end subroutine

      

