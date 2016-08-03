!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Kaps problem (Hairer II, pp 403)
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
      MODULE vanderPol_mod
      private
      public :: vanderPol
      contains

      subroutine vanderPol(nveclen,neq,ep,dt,tfinal,iDT,resE_vec,resI_vec,akk)

      use precision_vars,    only: wp
      use control_variables, only: temporal_splitting,probname,xjac,var_names,&
     &                             tol,dt_error_tol,uvec,uexact,programstep

      implicit none; save
!-----------------------VARIABLES----------------------------------------------
      integer,  parameter    :: vecl=2

      !INIT vars
      real(wp), intent(in   ) :: ep
      real(wp), intent(inout) :: dt
      integer,  intent(  out) :: nveclen,neq
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
          neq = vecl
          probname='vanderPol'     
          tol=1.0e-12_wp  
          dt_error_tol=5.0e-14_wp
          
          allocate(var_names(neq))
          var_names(:)=(/'Differential', 'Algebraic   '/)
        
        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')
      
          ! Time information
          dt = 0.5_wp/10**((iDT-1)/20.0_wp) ! timestep   
          tfinal = 0.5_wp                   ! final time
          
          !**Exact Solution** 
          open(unit=39,file='exact.vanderpol.data')
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
          uvec(1) = 2.0_wp
          uvec(2) = -0.6666654321121172_wp
!          uvec(2) = -2.0_wp/3.0_wp+10.0_wp/81.0_wp*ep-292.0_wp/2187.0_wp*ep**2 &
!     &            -1814.0_wp/19683.0_wp*ep**3

        case('BUILD_RHS')
          choose_RHS_type: select case (Temporal_Splitting)
            case('IMEX') ! For IMEX schemes
              resE_vec(1) = dt*uvec(2)
              resE_vec(2) = 0.0_wp
              resI_vec(1) = 0.0_wp
              resI_vec(2) = dt*((1-uvec(1)*uvec(1))*uvec(2) - uvec(1))/ep
            case('IMPLICIT') ! For fully implicit schemes
              resE_vec(:) = 0.0_wp
              resI_vec(1) = dt*uvec(2)
              resI_vec(2) = dt*((1-uvec(1)*uvec(1))*uvec(2) - uvec(1))/ep
          end select choose_RHS_type
        
        case('BUILD_JACOBIAN')
          choose_Jac_type: select case(temporal_splitting)
            case('IMEX') ! For IMEX schemes
              xjac(1,1) = 1.0_wp-akk*dt*(0.0_wp)
              xjac(1,2) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(2,1) = 0.0_wp-akk*dt*(-2*uvec(1)*uvec(2)-1)/ep
              xjac(2,2) = 1.0_wp-akk*dt*(1-uvec(1)*uvec(1))/ep
            case('IMPLICIT') ! For fully implicit schemes
              xjac(1,1) = 1.0_wp-akk*dt*(0.0_wp)
              xjac(1,2) = 0.0_wp-akk*dt*(1.0_wp)
              xjac(2,1) = 0.0_wp-akk*dt*(-2*uvec(1)*uvec(2)-1)/ep
              xjac(2,2) = 1.0_wp-akk*dt*(+1-uvec(1)*uvec(1))/ep
          end select choose_Jac_type
          
      end select Program_Step_Select
      end subroutine vanderPol
      end module vanderPol_mod
