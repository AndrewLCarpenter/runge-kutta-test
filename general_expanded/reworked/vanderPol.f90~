      subroutine vanderPol(stage,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE,resI,akk,xjac)
      implicit none

      integer,  parameter                      :: wp=8
      integer,  parameter                      :: vecl=4

      integer,                   intent(in   ) :: stage

      !INIT vars
      real(wp), dimension(vecl), intent(inout) :: uvec
      real(wp),                  intent(in   ) :: ep
      real(wp), dimension(vecl), intent(  out) :: uexact
      real(wp),                  intent(inout) :: dt
      integer,                   intent(inout) :: nveclen
      real(wp),                  intent(  out) :: tfinal
      integer,                   intent(in   ) :: iDT

      real(wp), dimension(81,vecl)             :: ExactTot
      real(wp)                                 :: diff
      integer                                  :: i

      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE,resI
      
      !Jacob vars
      real(wp),                       intent(in   ) :: akk
      real(wp), dimension(vecl,vecl), intent(  out) :: xjac


      if (stage==0) then
        open(unit=39,file='exact.vanderpol.data')
        rewind(39)
        do i=1,81
          read(39,*)ExactTot(i,1),ExactTot(i,2)
          ExactTot(i,3) = 1.0_wp/10**((i-1)/(10.0_wp))                  !  used for 81 values of ep
        enddo
        do i=1,81
          diff = abs(ExactTot(i,3) - ep)
          if(diff.le.1.0e-10_wp)then
            uexact(1) = ExactTot(i,1)
            uexact(2) = ExactTot(i,2)
            go to 100 
          endif
        enddo
 100    continue
        dt = 0.5_wp/10**((iDT-1)/20.0_wp)
        nvecLen = 2
        tfinal = 0.5_wp
        uvec(1) = 2.0_wp
        uvec(2) = -0.6666654321121172_wp

      elseif (stage==1.or.stage==2) then
        resE(1) = dt*uvec(2)
        resE(2) = 0.0_wp
        resI(1) = 0.0_wp
        resI(2) = dt*((1-uvec(1)*uvec(1))*uvec(2) - uvec(1))/ep

      elseif (stage==3) then
        xjac(1,1) = 1.-akk*dt*(0.)
        xjac(1,2) = 0.-akk*dt*(1.)
        xjac(2,1) = 0.-akk*dt*(-2*uvec(1)*uvec(2)-1)/ep
        xjac(2,2) = 1.-akk*dt*(1-uvec(1)*uvec(1))/ep
      endif

      return
      end subroutine

      

