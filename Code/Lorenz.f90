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

      subroutine Lorenz(nveclen,neq,ep,dt,tfinal,iDT,rese_vec,resi_vec,akk)

      use precision_vars,    only: wp
      use control_variables, only: temporal_splitting,probname,xjac,var_names,&
     &                             tol,dt_error_tol,uvec,uexact,programstep

      implicit none; save
!-----------------------VARIABLES----------------------------------------------
      integer,  parameter    :: vecl  = 3
      real(wp), parameter    :: sigma0= 10.0_wp
      real(wp), parameter    :: beta0 = 8.0_wp/3.0_wp
      real(wp), parameter    :: rho0  = 28.0_wp
                 
      !INIT vars
      real(wp),        intent(in   ) :: ep
      real(wp),        intent(inout) :: dt
      integer,         intent(  out) :: nveclen,neq
      real(wp),        intent(  out) :: tfinal
      integer,         intent(in   ) :: iDT

      real(wp), dimension(81,vecl+1) :: ExactTot
      real(wp)                       :: diff
      real(wp)                       :: x0,y0,z0,x1,y1,z1
      real(wp)                       :: sigma,beta,rho
      integer                        :: i

      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: rese_vec,resi_vec
      
      !Jacob vars
      real(wp), intent(in   ) :: akk

      sigma =  sigma0/ ep
      beta  =  beta0
      rho   =  rho0
      
!------------------------------------------------------------------------------
      Program_Step_Select: select case(programstep)
        !**Pre-initialization. Get problem name and vector length**
        case('INITIALIZE_PROBLEM_INFORMATION')
          nvecLen = vecl
          neq = vecl
          probname='Lorenz   '         
          tol=1.0e-12_wp
          dt_error_tol=1.0e-11_wp
          
          allocate(var_names(neq))
          var_names(:)=(/'Algebraic   ', 'Differential', 'Differential'/)
          
        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')
      
          tfinal = 1.0_wp                    ! final time

          !Timestep information
          choose_dt: select case(temporal_splitting)
            case('EXPLICIT'); dt = 0.005_wp/10**((iDT-1)/20.0_wp) ! explicit timestep
            case default    ; dt = 0.250_wp/10**((iDT-1)/20.0_wp) ! implicit timestep      
          end select choose_dt

!         dt = 0.25_wp*0.00001_wp/10**((iDT-1)/20.0_wp) !used for exact solution
!         dt = 0.25_wp/10**((iDT-1)/20.0_wp) ! timestep 

          !**Exact Solution**
          open(unit=39,file='exact.Lorenz.data')
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
!         uvec(1) = 0.0_wp
!         uvec(2) = 1.0_wp
!         uvec(3) = 0.0_wp

          z0 = rho - 1.0_wp
          x0 = sqrt(beta*z0)
          y0 = x0
          z1 =  0.0_wp
          x1 = +1.0_wp
          y1 = x1

          uvec(1) = x0 + ep * x1
          uvec(2) = y0 + ep * y1
          uvec(3) = z0 + ep * z1
        
        case('BUILD_RHS')
          choose_RHS_type: select case (Temporal_Splitting)
            case('EXPLICIT') ! For fully explicit schemes
              resE_vec(1)=dt*sigma*(uvec(2)-uvec(1))
              resE_vec(2)=dt*(-uvec(1)*uvec(3)+rho*uvec(1)-uvec(2))
              resE_vec(3)=dt*(+uvec(1)*uvec(2)-beta*uvec(3))
              resI_vec(:)=0.0_wp
            case('IMEX') ! For IMEX schemes
              resE_vec(1)=0.0_wp
              resE_vec(2)=dt*(-uvec(1)*uvec(3)+rho*uvec(1)-uvec(2))
              resE_vec(3)=dt*(+uvec(1)*uvec(2)-beta*uvec(3))     
              resI_vec(1)=dt*sigma*(uvec(2)-uvec(1))
              resI_vec(2)=0.0_wp
              resI_vec(3)=0.0_wp            
            case('IMPLICIT') ! For fully implicit schemes
              resE_vec(:)=0.0_wp
              resI_vec(1)=dt*(sigma*(uvec(2)-uvec(1)))
              resI_vec(2)=dt*( (rho - uvec(3))*uvec(1) - uvec(2) )
              resI_vec(3)=dt*(+uvec(1)*uvec(2)-beta*uvec(3))
          end select choose_RHS_type
          
        case('BUILD_JACOBIAN')          
          choose_Jac_type: select case (Temporal_Splitting)
            case('EXPLICIT') ! For fully explicit schemes
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
              xjac(1,1) = 1.0_wp-akk*dt*(-sigma)
              xjac(1,2) = 0.0_wp-akk*dt*(+sigma)
              xjac(1,3) = 0.0_wp-akk*dt*(0.0_wp)

              xjac(2,1) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(2,2) = 1.0_wp-akk*dt*(0.0_wp)
              xjac(2,3) = 0.0_wp-akk*dt*(0.0_wp)

              xjac(3,1) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(3,2) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(3,3) = 1.0_wp-akk*dt*(0.0_wp)
            case('IMPLICIT') ! For fully implicit schemes
              xjac(1,1) = 1.0_wp-akk*dt*(-sigma)
              xjac(1,2) = 0.0_wp-akk*dt*(+sigma)
              xjac(1,3) = 0.0_wp-akk*dt*(0)

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
