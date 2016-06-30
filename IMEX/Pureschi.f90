      subroutine Pureschi(programStep,nveclen,ep,dt, &
                         &  tfinal,iDT,resE_vec,resI_vec,akk)

      use precision_vars
      use control_variables

      implicit none

      integer,  parameter                      :: vecl=2

      integer,                   intent(in   ) :: programStep

      !INIT vars
      !real(wp), dimension(vecl), intent(inout) :: uvec
      real(wp),                  intent(in   ) :: ep
    !  real(wp), dimension(vecl), intent(  out) :: uexact
      real(wp),                  intent(inout) :: dt
      integer,                   intent(inout) :: nveclen
      real(wp),                  intent(  out) :: tfinal
      integer,                   intent(in   ) :: iDT

      real(wp), dimension(81,vecl+1)             :: ExactTot
      real(wp)                                 :: diff
      integer                                  :: i

      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE_vec,resI_vec
      
      !Jacob vars
      real(wp),                       intent(in   ) :: akk
      !real(wp), dimension(vecl,vecl), intent(  out) :: xjac


      if (programStep==-1) then
        nvecLen = vecl
        probname='Pureschi '   
        tol=1.0e-12_wp
      elseif (programStep==0) then
        probname='Pureschi '
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

      elseif (programStep>=1 .and. programStep<=3) then

        select case (Temporal_Splitting)

          case('IMEX')
                if (programStep==1 .or.programStep==2) then
              resE_vec(1) = dt*(-uvec(2))
              resE_vec(2) = dt*(+uvec(1))
              resI_vec(1) = 0.0_wp
              resI_vec(2) = dt*(sin(uvec(1)) - uvec(2))/ep
            elseif (programStep==3) then
              xjac(1,1) = 1.0_wp-akk*dt*(0.0_wp)
              xjac(1,2) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(2,1) = 0.0_wp-akk*dt*(cos(uvec(1)) )/ep
              xjac(2,2) = 1.0_wp-akk*dt*(-1.0_wp)/ep
            endif
          case('IMPLICIT')
                if (programStep==1 .or.programStep==2) then
              resE_vec(1) = 0.0_wp
              resE_vec(2) = 0.0_wp
              resI_vec(1) = dt*(-uvec(2))
              resI_vec(2) = dt*(+uvec(1) + (sin(uvec(1)) - uvec(2))/ep)
            elseif (programStep==3) then
              xjac(1,1) = 1.0_wp-akk*dt*(+0.0_wp)
              xjac(1,2) = 0.0_wp-akk*dt*(-1.0_wp)
              xjac(2,1) = 0.0_wp-akk*dt*(+1.0_wp + cos(uvec(1))/ep)
              xjac(2,2) = 1.0_wp-akk*dt*(-1.0_wp)/ep
            endif
        end select
      endif
      
      return
      end subroutine

