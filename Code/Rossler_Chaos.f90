!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Rossler_Chaos problem 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
      module Rossler_mod
      private
      public :: Rossler_Chaos
      contains
      subroutine Rossler_Chaos(nveclen,neq,ep,dt,tfinal,iDT,resE_vec,resI_vec,akk)
      use precision_vars,    only: wp
      use control_variables, only: temporal_splitting,probname,xjac,var_names,&
     &                             tol,dt_error_tol,uvec,uexact,programstep

!     Rossler system
!     x_t = - y - z
!     y_t = + x + aa y
!     z_t = + bb * x + (-cc + x) z
!
      implicit none; save
!-----------------------VARIABLES----------------------------------------------
      integer, parameter     :: vecl=3
      real(wp), parameter    :: aa= 0.15_wp
      real(wp), parameter    :: bb= 0.20_wp
      real(wp), parameter    :: cc=10.00_wp
 
      !INIT vars
      real(wp),        intent(in   ) :: ep
      real(wp),        intent(inout) :: dt
      integer,         intent(  out) :: nveclen,neq
      real(wp),        intent(  out) :: tfinal
      integer,         intent(in   ) :: iDT

      real(wp), dimension(81,vecl+1) :: ExactTot
      real(wp)                       :: diff
      integer                        :: i

      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE_vec,resI_vec
      
      !Jacob vars
      real(wp), intent(in   ) :: akk
      real(wp)                :: cct
!------------------------------------------------------------------------------

      Program_Step_Select: select case(programstep)
        !**Pre-initialization. Get problem name and vector length**
        case('INITIALIZE_PROBLEM_INFORMATION')
          nvecLen = vecl
          neq = vecl
          probname='Rossler_3'         
          tol=1.0e-14_wp
          dt_error_tol=2.0e-12_wp
          
          allocate(var_names(neq))
          var_names(:)=(/'Differential', 'Differential', 'Algebraic   '/)
          
        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')

         !Time information
         !dt = 0.25_wp*0.00001_wp/10**((iDT-1)/20.0_wp) !used for exact solution
          dt = 1.0_wp/10**((iDT-1)/20.0_wp) ! timestep 
          tfinal = 1.0_wp                    ! final time

          !**Exact Solution**
          open(unit=39,file='./Exact_Data/exact.Rossler_Chaos.data')
          rewind(39)
          do i=1,81
            read(39,*)ExactTot(i,1),ExactTot(i,2),ExactTot(i,3)
            ExactTot(i,4) = 1.0_wp/10**((i-1)/(10.0_wp)) !  used for 81 values of ep
          enddo
          do i=1,81
            diff = abs(ExactTot(i,4) - ep)
            if(diff.le.1.0e-10_wp)then
              uexact(:) = ExactTot(i,:vecl)
              exit
            endif
          enddo

          !**IC**
          uvec(1) = 0.0_wp
          uvec(2) = 1.0_wp
          uvec(3) = 0.0_wp
        
        case('BUILD_RHS')
          !  Stiff component is cc
          cct = cc / ep
          choose_RHS_type: select case (Temporal_Splitting)
            case('IMEX') ! For IMEX schemes
              resE_vec(1)=dt * (-uvec(2)-uvec(3)) ;
              resE_vec(2)=dt * (+uvec(1)+aa*uvec(2)) ;
              resE_vec(3)=dt * (bb*uvec(1) + uvec(3)*(    +uvec(1))) ;
              resI_vec(1:2)=0.0_wp
              resI_vec(3)=dt * (           + uvec(3)*(-cct        )) ;
            case('IMPLICIT','FIRK') ! For fully implicit schemes
              resE_vec(:)=0.0_wp
              resI_vec(1)=dt * (-uvec(2)-uvec(3)) ;
              resI_vec(2)=dt * (+uvec(1)+aa*uvec(2)) ;
              resI_vec(3)=dt * (bb*uvec(1) + uvec(3)*(-cct+uvec(1))) ;
          end select choose_RHS_type

         case('BUILD_JACOBIAN')
          !  Stiff component is cc
          cct = cc / ep
          choose_Jac_type: select case (Temporal_Splitting)
            case('IMEX') ! For IMEX schemes
              xjac(1,1) = 1.0_wp-akk*dt*(+0.0_wp)
              xjac(1,2) = 0.0_wp-akk*dt*(+0.0_wp)
              xjac(1,3) = 0.0_wp-akk*dt*(+0.0_wp)

              xjac(2,1) = 0.0_wp-akk*dt*(+0.0_wp)
              xjac(2,2) = 1.0_wp-akk*dt*(+0.0_wp)
              xjac(2,3) = 0.0_wp-akk*dt*(+0.0_wp)

              xjac(3,1) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(3,2) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(3,3) = 1.0_wp-akk*dt*( -cct )
            case('IMPLICIT') ! For fully implicit schemes
              xjac(1,1) = 1.0_wp-akk*dt*(+0.0_wp)
              xjac(1,2) = 0.0_wp-akk*dt*(-1.0_wp)
              xjac(1,3) = 0.0_wp-akk*dt*(-1.0_wp)

              xjac(2,1) = 0.0_wp-akk*dt*(+1.0_wp)
              xjac(2,2) = 1.0_wp-akk*dt*(aa )
              xjac(2,3) = 0.0_wp-akk*dt*(0.0_wp)

              xjac(3,1) = 0.0_wp-akk*dt*(+bb + uvec(3))
              xjac(3,2) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(3,3) = 1.0_wp-akk*dt*(-cct + uvec(1))
            case('FIRK') ! For fully implicit schemes
              xjac(1,1) = (+0.0_wp)
              xjac(1,2) = (-1.0_wp)
              xjac(1,3) = (-1.0_wp)

              xjac(2,1) = (+1.0_wp)
              xjac(2,2) = (aa )
              xjac(2,3) = (0.0_wp)

              xjac(3,1) = (+bb + uvec(3))
              xjac(3,2) = (0.0_wp)
              xjac(3,3) = (-cct + uvec(1))
          end select choose_Jac_type
          
      end select Program_step_select
      end subroutine Rossler_Chaos
      end module Rossler_mod
