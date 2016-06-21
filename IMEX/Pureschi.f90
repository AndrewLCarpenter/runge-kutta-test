      subroutine Pureschi(programStep,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE,resI,akk,xjac)

      use precision_vars

      implicit none

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

      real(wp), dimension(81,vecl)             :: ExactTot
      real(wp)                                 :: diff
      integer                                  :: i

      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE,resI
      
      !Jacob vars
      real(wp),                       intent(in   ) :: akk
      real(wp), dimension(vecl,vecl), intent(  out) :: xjac


      if (programStep==0) then
        open(unit=39,file='exact.pureschi.1.data')
        rewind(39)
        do i=1,81
          read(39,*)ExactTot(i,1),ExactTot(i,2)
          ExactTot(i,3) = 1.0_wp/10**((i-1)/(10.0_wp)) !  used for 81 values of ep
          enddo

        do i=1,81
          diff = abs(ExactTot(i,3) - ep)
          if(diff.le.1.0e-10_wp)then
            uexact(1) = ExactTot(i,1)
            uexact(2) = ExactTot(i,2)
            go to 2
          endif
        enddo
   2    continue

        dt = 0.5_wp/10**((iDT-1)/20._wp)
        nvecLen = 2
        tfinal = 5.0_wp

!       pi = acos(-1.0_wp)  brought in through ``precision_vars''

!  IC: problem 1   :  equilibrium IC
          uvec(1) = pi/2.0_wp
          uvec(2) = sin(pi/2.0_wp)

      elseif (programStep==1 .or.programStep==2) then
!may have created an error when splitting into E/I
        resE(1) = dt*(-uvec(2))
        resE(2) = dt* uvec(1)
        resI(1) = 0.0_wp
        resI(2) = dt*(sin(uvec(1)) - uvec(2))/ep

      elseif (programStep==3) then
        xjac(1,1) = 1.-akk*dt*(0.)
        xjac(1,2) = 0.-akk*dt*(0.)
        xjac(2,1) = 0.-akk*dt*(cos(uvec(1)) )/ep
        xjac(2,2) = 1.-akk*dt*(-1.)/ep
      endif

      return
      end subroutine

      

