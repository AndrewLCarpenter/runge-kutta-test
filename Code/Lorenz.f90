!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Lorenz problem 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
      module Lorenz_mod
      private
      public :: Lorenz
      contains
      subroutine Lorenz(nveclen,ep,dt,tfinal,iDT,rese_vec,resi_vec,akk)
      use precision_vars,    only: wp
      use control_variables, only: temporal_splitting,probname,xjac, &
     &                             tol,dt_error_tol,uvec,uexact,programstep

      implicit none; save
!-----------------------VARIABLES----------------------------------------------
      integer, parameter     :: vecl=3
      real(wp), parameter    :: sigma=5.0_wp
      real(wp), parameter    :: beta=1.0_wp/3.0_wp
      real(wp), parameter    :: rho=2.0_wp
                 
      !INIT vars
      real(wp),        intent(in   ) :: ep
      real(wp),        intent(inout) :: dt
      integer,         intent(  out) :: nveclen
      real(wp),        intent(  out) :: tfinal
      integer,         intent(in   ) :: iDT

      real(wp), dimension(81,vecl+1) :: ExactTot
      real(wp)                       :: diff
      integer                        :: i,j

      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: rese_vec,resi_vec
      
      !Jacob vars
      real(wp), intent(in   ) :: akk
!------------------------------------------------------------------------------
      
      Program_Step_Select: select case(programstep)
        !**Pre-initialization. Get problem name and vector length**
        case('INITIALIZE_PROBLEM_INFORMATION')
          nvecLen = vecl
          probname='Lorenz   '         
          tol=1.0e-14_wp
          dt_error_tol=1.0e-11_wp
          
        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')
      
          !Time information
          !dt = 0.25_wp*0.00001_wp/10**((iDT-1)/20.0_wp) !used for exact solution
          dt = 0.25_wp/10**((iDT-1)/20.0_wp) ! timestep 
          tfinal = 1.0_wp                    ! final time

          !**Exact Solution**
          open(unit=39,file='exact.lorenz.data')
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
          choose_RHS_type: select case (Temporal_Splitting)
            case('EXPLICIT') ! For fully explicit schemes
              resI_vec(:)=0.0_wp
              resE_vec(1)=dt*sigma*(uvec(2)-uvec(1))/ep
              resE_vec(2)=dt*(-uvec(1)*uvec(3)+rho*uvec(1)-uvec(2))
              resE_vec(3)=dt*(uvec(1)*uvec(2)-beta*uvec(3))
            case('IMEX') ! For IMEX schemes
              rese_vec(1)=0.0_wp
              rese_vec(2)=dt*(-uvec(1)*uvec(3)+rho*uvec(1)-uvec(2))
              rese_vec(3)=dt*(uvec(1)*uvec(2)-beta*uvec(3))     
              resi_vec(1)=dt*sigma*(uvec(2)-uvec(1))/ep
              resi_vec(2)=0.0_wp
              resi_vec(3)=0.0_wp            
            case('IMPLICIT') ! For fully implicit schemes
              rese_vec(:)=0.0_wp
              resi_vec(1)=dt*sigma*(uvec(2)-uvec(1))/ep
              resi_vec(2)=dt*(-uvec(1)*uvec(3)+rho*uvec(1)-uvec(2))
              resi_vec(3)=dt*(uvec(1)*uvec(2)-beta*uvec(3))
          end select choose_RHS_type
          
        case('BUILD_JACOBIAN')          
          choose_Jac_type: select case (Temporal_Splitting)
            case('EXPLICIT') ! For fully implicit schemes
              xjac(1,1) = 1.0_wp
              xjac(1,2) = 0.0_wp
              xjac(1,3) = 0.0_wp

              xjac(2,1) = 0.0_wp
              xjac(2,2) = 1.0_wp
              xjac(2,3) = 0.0_wp

              xjac(3,1) = 0.0_wp
              xjac(3,2) = 0.0_wp
              xjac(3,3) = 1.0_wp
            case('IMEX') ! For IMEX schemes
              xjac(1,1) = 1.0_wp-akk*dt*(-sigma)/ep
              xjac(1,2) = 0.0_wp-akk*dt*(+sigma)/ep
              xjac(1,3) = 0.0_wp-akk*dt*(0.0_wp)/ep

              xjac(2,1) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(2,2) = 1.0_wp-akk*dt*(0.0_wp)
              xjac(2,3) = 0.0_wp-akk*dt*(0.0_wp)

              xjac(3,1) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(3,2) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(3,3) = 1.0_wp-akk*dt*(0.0_wp)
            case('IMPLICIT') ! For fully implicit schemes
              xjac(1,1) = 1.0_wp-akk*dt*(-sigma)/ep
              xjac(1,2) = 0.0_wp-akk*dt*(+sigma)/ep
              xjac(1,3) = 0.0_wp-akk*dt*(0)/ep

              xjac(2,1) = 0.0_wp-akk*dt*(rho-uvec(3))
              xjac(2,2) = 1.0_wp-akk*dt*(-1.0_wp)
              xjac(2,3) = 0.0_wp-akk*dt*(-uvec(1))

              xjac(3,1) = 0.0_wp-akk*dt*(uvec(2))
              xjac(3,2) = 0.0_wp-akk*dt*(uvec(1))
              xjac(3,3) = 1.0_wp-akk*dt*(-beta)
          end select choose_Jac_type
          
      end select Program_Step_Select
      end subroutine Lorenz
      end module Lorenz_mod
