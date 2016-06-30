      subroutine Kaps(programStep,nveclen,ep,dt, &
                     & tfinal,iDT,resE_vec,resi_vec,akk)

      use precision_vars
      use control_variables

      implicit none

      integer,  parameter                      :: vecl=2

      integer,                   intent(in   ) :: programStep

      !INIT vars
      !real(wp), dimension(vecl), intent(inout) :: uvec
      real(wp),                  intent(in   ) :: ep
   !   real(wp), dimension(vecl), intent(  out) :: uexact
      real(wp),                  intent(inout) :: dt
      integer,                   intent(inout) :: nveclen
      real(wp),                  intent(  out) :: tfinal
      integer,                   intent(in   ) :: iDT

      real(wp)                                 :: tmp, epI

      !RHS vars
      real(wp), dimension(vecl+1), intent(  out) :: resE_vec,resi_vec
      
      !Jacob vars
      real(wp),                       intent(in   ) :: akk
      !real(wp), dimension(vecl,vecl), intent(  out) :: xjac

      epI = 1.0_wp / ep 

      if (programStep==-1) then
        nvecLen = vecl
        probname='Kaps     '  
        tol=1.0e-12_wp
      elseif (programStep==0) then
        dt = 0.5_wp/10**((iDT-1)/20.0_wp)
        tfinal = 1.0_wp
        tmp = exp(-tfinal)
        uexact(1) = tmp*tmp
        uexact(2) = tmp

!  IC: problem 1   :  equilibrium IC
        uvec(1) = 1.0_wp
        uvec(2) = 1.0_wp

      elseif (programStep>=1 .and. programStep<=3) then
        select case (Temporal_Splitting)

          case('IMEX')
                if (programStep==1 .or.programStep==2) then
              resE_vec(1) = dt*(-2.0_wp*uvec(1))
              resE_vec(2) = dt*(uvec(1) - uvec(2) - uvec(2)*uvec(2) )
              resi_vec(1) = dt*(-epI*uvec(1) + epI*uvec(2)*uvec(2))
              resi_vec(2) = 0.0_wp
            elseif (programStep==3) then
              xjac(1,1) = 1.0_wp-akk*dt*(-epI)
              xjac(1,2) = 0.0_wp-akk*dt*(+epI)*2.0_wp*uvec(2)
              xjac(2,1) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(2,2) = 1.0_wp-akk*dt*(0.0_wp)
            endif
          case('IMPLICIT')
                if (programStep==1 .or.programStep==2) then
              resE_vec(1) = 0.0_wp
              resE_vec(2) = 0.0_wp
              resi_vec(1) = dt*(-(epI+2.0_wp)*uvec(1) + epI*uvec(2)*uvec(2))
              resi_vec(2) = dt*(uvec(1) - uvec(2) - uvec(2)*uvec(2) )
            elseif (programStep==3) then
              xjac(1,1) = 1.0_wp-akk*dt*(-(epI+2.0_wp))
              xjac(1,2) = 0.0_wp-akk*dt*(+epI*2.0_wp*uvec(2))
              xjac(2,1) = 0.0_wp-akk*dt*(1.0_wp)
              xjac(2,2) = 1.0_wp-akk*dt*(-(1.0_wp+2.0_wp*uvec(2)))
            endif
          end select
      endif
      
      return
      end subroutine

      

