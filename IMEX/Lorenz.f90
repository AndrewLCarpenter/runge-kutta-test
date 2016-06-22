      subroutine Lorenz(programStep,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE,resI,akk,xjac)

      use precision_vars

      implicit none

      integer, parameter                       :: vecl=4
      real(wp), parameter                      :: sigma=10,beta=8/3,rho=28
                 
      integer,                   intent(in   ) :: programStep
 
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
      integer                                  :: i,j

      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE,resI
      
      !Jacob vars
      real(wp),                       intent(in   ) :: akk
      real(wp), dimension(vecl,vecl), intent(  out) :: xjac


      if (programStep==0) then
        dt = 0.5_wp/10**((iDT-1)/20.0_wp)
        nveclen = 3 
        tfinal = 1.0_wp !!arbitrary

        open(unit=39,file='exact.lorenz.data')
        rewind(39)
        do i=1,81
          read(39,*)ExactTot(i,1),ExactTot(i,2),ExactTot(i,3)
          ExactTot(i,4) = 1.0_wp/10**((i-1)/(10.0_wp))                  !  used for 81 values of ep
        enddo
        do i=1,81
          diff = abs(ExactTot(i,4) - ep)
          if(diff.le.1.0e-10_wp)then
            uexact(1) = ExactTot(i,1)
            uexact(2) = ExactTot(i,2)
            uexact(3) = ExactTot(i,3)
            go to 100 
          endif
        enddo
 100    continue

        uvec(1) = 0.0_wp
        uvec(2) = 1.0_wp
        uvec(3) = 0.0_wp

      elseif (programStep==1 .or.programStep==2) then
        resE(1)=0.0_wp
        resE(2)=dt*(-uvec(1)*uvec(3)+rho*uvec(1)-uvec(2))
        resE(3)=dt*(uvec(1)*uvec(2)-beta*uvec(3))     
        resI(1)=dt*sigma*(uvec(2)-uvec(1))/ep
        resI(2)=0.0_wp
        resI(3)=0.0_wp

      elseif (programStep==3) then
        xjac(1,1) = 1.-akk*dt*(-sigma)/ep
        xjac(1,2) = 0.-akk*dt*(sigma)/ep
        xjac(1,3) = 0.-akk*dt*(0.0_wp)/ep

        xjac(2,1) = 0.-akk*dt*(0.0_wp)
        xjac(2,2) = 1.-akk*dt*(0.0_wp)
        xjac(2,3) = 0.-akk*dt*(0.0_wp)

        xjac(3,1) = 0.-akk*dt*(0.0_wp)
        xjac(3,2) = 0.-akk*dt*(0.0_wp)
        xjac(3,3) = 1.-akk*dt*(0.0_wp)
      endif

      return
      end subroutine

      

