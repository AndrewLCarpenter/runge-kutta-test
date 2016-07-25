!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate Brusselator
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
      module Brusselator_mod
      private
      public :: Brusselator
      contains
      subroutine Brusselator(nveclen,ep,dt,tfinal,iDT,resE_vec,resI_vec,akk)

      use precision_vars,    only: wp
      use control_variables, only: temporal_splitting,probname,xjac, &
     &                             tol,dt_error_tol,uvec,uexact,programstep

      implicit none; save
!-----------------------VARIABLES----------------------------------------------
      integer,  parameter    :: vecl=2
      real(wp), parameter    :: aa=1.0_wp, bb=1.7_wp

      !INIT vars
      real(wp), intent(in   ) :: ep
      real(wp), intent(inout) :: dt
      integer,  intent(  out) :: nveclen
      real(wp), intent(  out) :: tfinal
      integer,  intent(in   ) :: iDT

      real(wp), dimension(81,vecl+1) :: ExactTot
      real(wp)                       :: diff
      integer                        :: i

      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE_vec,resI_vec
      
      !Jacob vars
      real(wp), intent(in   ) :: akk
!------------------------------------------------------------------------------
      
      Program_Step_Select: select case(programstep)
        !**Pre-initialization. Get problem name and vector length**
        case('INITIALIZE_PROBLEM_INFORMATION')
          nvecLen = vecl
          probname='Brusselat'     
          tol=9.9e-11_wp  
          dt_error_tol=1.0e-13_wp
          
        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')
        
          !Time information
          dt = 0.25_wp/10**((iDT-1)/20.0_wp) ! timestep   
          tfinal = 1.0_wp                   ! final time
          
          !**Exact Solution** 
          open(unit=39,file='exact.brusselator4.data')
          rewind(39)
          do i=1,81
            read(39,*)ExactTot(i,1),ExactTot(i,2)
            ExactTot(i,3) = 1.0_wp/10**((i-1)/(10.0_wp))  !  used for 81 values of ep
          enddo
          do i=1,81
            diff = abs(ExactTot(i,3) - ep)
            if(diff.le.1.0e-10_wp)then
              uexact(:) = ExactTot(i,:vecl)
              exit
            endif
          enddo

          !**IC**
          uvec(1) = 1.0_wp
          uvec(2) = 1.0_wp
      
        case('BUILD_RHS')
          choose_RHS_type: select case (Temporal_Splitting)
            case('IMEX') ! For IMEX schemes
              resE_vec(1) = dt*(aa-uvec(1)*bb-uvec(1))
              resE_vec(2) = dt*(uvec(1)*bb)
              resI_vec(1) = dt*( uvec(1)*uvec(1)*uvec(2)*ep)
              resI_vec(2) = dt*(-uvec(1)*uvec(1)*uvec(2)*ep)
            case('IMPLICIT') ! For fully implicit schemes
              resE_vec(:) = 0.0_wp
              resI_vec(1) = dt*(aa+uvec(1)*uvec(1)*uvec(2)/ep-uvec(1)*bb-uvec(1))
              resI_vec(2) = dt*(  -uvec(1)*uvec(1)*uvec(2)/ep+uvec(1)*bb)
          end select choose_RHS_type
       
        case('BUILD_JACOBIAN')
          choose_Jac_type: select case(temporal_splitting)
            case('IMEX') ! For IMEX schemes
              xjac(1,1) = 1.0_wp-akk*dt*(2.0_wp*uvec(1)*uvec(2)*ep)
              xjac(1,2) = 0.0_wp-akk*dt*(       uvec(1)*uvec(1)*ep)
              xjac(2,1) = 0.0_wp-akk*dt*(2.0_wp*uvec(1)*uvec(2)*ep)
              xjac(2,2) = 1.0_wp-akk*dt*(      -uvec(1)*uvec(1)*ep)
            case('IMPLICIT') ! For fully implicit schemes
              xjac(1,1) = 1.0_wp-akk*dt*(   2.0_wp*uvec(1)*uvec(2)/ep-bb-1.0_wp)
              xjac(1,2) = 0.0_wp-akk*dt*(          uvec(1)*uvec(1)/ep)
              xjac(2,1) = 0.0_wp-akk*dt*(bb-2.0_wp*uvec(1)*uvec(2)/ep)
              xjac(2,2) = 1.0_wp-akk*dt*(         -uvec(1)*uvec(1)/ep)
          end select choose_Jac_type
        
      end select Program_Step_select
      end subroutine Brusselator
      end module Brusselator_mod
