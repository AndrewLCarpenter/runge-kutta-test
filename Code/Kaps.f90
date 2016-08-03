!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Kaps problem (Dekker 7.5.2 pp. 215 (Kaps problem   : Index 1)) 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
      module Kaps_mod
      private
      public :: Kaps
      contains
      subroutine Kaps(nveclen,neq,ep,dt,tfinal,iDT,resE_vec,resi_vec,akk)

      use precision_vars,    only: wp
      use control_variables, only: temporal_splitting,probname,xjac,var_names,&
     &                             tol,dt_error_tol,uvec,uexact,programstep

      implicit none; save
!-----------------------VARIABLES----------------------------------------------
      integer, parameter     :: vecl=2

      !INIT vars
      real(wp), intent(in   ) :: ep
      real(wp), intent(inout) :: dt
      integer,  intent(  out) :: nveclen,neq
      real(wp), intent(  out) :: tfinal
      integer,  intent(in   ) :: iDT

      real(wp)                :: tmp, epI

      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE_vec,resi_vec
      
      !Jacob vars
      real(wp), intent(in   ) :: akk
!------------------------------------------------------------------------------
    
      Program_Step_Select: select case(programstep)
        !**Pre-initialization. Get problem name and vector length**
        case('INITIALIZE_PROBLEM_INFORMATION')
          nvecLen = vecl
          neq = vecl
          probname='Kaps     '  
          tol=1.0e-12_wp
          dt_error_tol=1.0e-11_wp
              
          allocate(var_names(neq))
          var_names(:)=(/'Algebraic   ', 'Differential'/)    
          
        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')
       
        !Time information
        dt = 0.5_wp/10**((iDT-1)/20.0_wp) ! timestep
        tfinal = 1.0_wp                   ! final time

        !*Exact Solution**        
        tmp = exp(-tfinal)
        uexact(1) = tmp*tmp
        uexact(2) = tmp

        !**Equilibrium IC**
        uvec(1) = 1.0_wp
        uvec(2) = 1.0_wp

        case('BUILD_RHS')
          epI = 1.0_wp / ep  !**Initialize 1/epsilon
          choose_RHS_type: select case (Temporal_Splitting)
            case('IMEX') ! For IMEX schemes
              resE_vec(1) = dt*(-2.0_wp*uvec(1))
              resE_vec(2) = dt*(uvec(1) - uvec(2) - uvec(2)*uvec(2) )
              resi_vec(1) = dt*(-epI*uvec(1) + epI*uvec(2)*uvec(2))
              resi_vec(2) = 0.0_wp
            case('IMPLICIT') ! For fully implicit schemes
              resE_vec(:) = 0.0_wp
              resi_vec(1) = dt*(-(epI+2.0_wp)*uvec(1) + epI*uvec(2)*uvec(2))
              resi_vec(2) = dt*(uvec(1) - uvec(2) - uvec(2)*uvec(2) )
          end select choose_RHS_type
        
        case('BUILD_JACOBIAN')
          epI = 1.0_wp / ep  !**Initialize 1/epsilon
          choose_Jac_type: select case (Temporal_Splitting)
            case('IMEX') ! For IMEX schemes
              xjac(1,1) = 1.0_wp-akk*dt*(-epI)
              xjac(1,2) = 0.0_wp-akk*dt*( epI)*2.0_wp*uvec(2)
              xjac(2,1) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(2,2) = 1.0_wp-akk*dt*(0.0_wp)            
            case('IMPLICIT') ! For fully implicit schemes
              xjac(1,1) = 1.0_wp-akk*dt*(-(epI+2.0_wp))
              xjac(1,2) = 0.0_wp-akk*dt*(+epI*2.0_wp*uvec(2))
              xjac(2,1) = 0.0_wp-akk*dt*(1.0_wp)
              xjac(2,2) = 1.0_wp-akk*dt*(-(1.0_wp+2.0_wp*uvec(2)))
          end select choose_Jac_type
          
      end select Program_Step_Select
      end subroutine Kaps
      end module Kaps_mod
