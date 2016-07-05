!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Pureschi and Russo problem 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
      subroutine Pureschi(programStep,nveclen,ep,dt, &
     &                    tfinal,iDT,resE_vec,resI_vec,akk)

      use precision_vars
      use control_variables

      implicit none

      integer,  parameter    :: vecl=2
      integer, intent(in   ) :: programStep

      !INIT vars
      real(wp),                  intent(in   ) :: ep
      real(wp),                  intent(inout) :: dt
      integer,                   intent(inout) :: nveclen
      real(wp),                  intent(  out) :: tfinal
      integer,                   intent(in   ) :: iDT

      real(wp), dimension(81,vecl+1)           :: ExactTot
      real(wp)                                 :: diff
      integer                                  :: i

      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE_vec,resI_vec
      
      !Jacob vars
      real(wp), intent(in   ) :: akk
!------------------------------------------------------------------------------

      !**Pre-initialization. Get problem name and vector length**
      if (programStep==-1) then
        nvecLen = vecl
        probname='Pureschi '   
        tol=1.0e-12_wp
        dt_error_tol=1.0e-13_wp
        
      !**Initialization of problem information**        
      elseif (programStep==0) then
      
        dt = 0.5_wp/10**((iDT-1)/20._wp) ! timestep
        tfinal = 5.0_wp                  ! final time

        !**Exact Solution**
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
            exit
          endif
        enddo

        !**Equilibrium IC**
        uvec(1) = pi/2.0_wp
        uvec(2) = sin(pi/2.0_wp)

      !**RHS and Jacobian**
      elseif (programStep>=1) then

        select case (Temporal_Splitting)

          case('IMEX') ! For IMEX schemes
            !**RHS**          
            if (programStep==1 .or.programStep==2) then
              resE_vec(1) = dt*(-uvec(2))
              resE_vec(2) = dt*(+uvec(1))
              resI_vec(1) = 0.0_wp
              resI_vec(2) = dt*(sin(uvec(1)) - uvec(2))/ep
            !**Jacobian**
            elseif (programStep==3) then
              xjac(1,1) = 1.0_wp-akk*dt*(0.0_wp)
              xjac(1,2) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(2,1) = 0.0_wp-akk*dt*(cos(uvec(1)) )/ep
              xjac(2,2) = 1.0_wp-akk*dt*(-1.0_wp)/ep
            endif
            
          case('IMPLICIT') ! For fully implicit schemes
            !**RHS**
            if (programStep==1 .or.programStep==2) then
              resE_vec(:) = 0.0_wp
              resI_vec(1) = dt*(-uvec(2))
              resI_vec(2) = dt*(+uvec(1) + (sin(uvec(1)) - uvec(2))/ep)
            !**Jacobian**
            elseif (programStep==3) then
              xjac(1,1) = 1.0_wp-akk*dt*(+0.0_wp)
              xjac(1,2) = 0.0_wp-akk*dt*(-1.0_wp)
              xjac(2,1) = 0.0_wp-akk*dt*(+1.0_wp + cos(uvec(1))/ep)
              xjac(2,2) = 1.0_wp-akk*dt*(-1.0_wp)/ep
            endif
        
          case default ! To catch invald inputs
            write(*,*)'Invaild case entered. Enter "IMEX" or "IMPLICIT"'
            write(*,*)'Exiting'
            stop
            
        end select
          
      endif
      
      end subroutine Pureschi

