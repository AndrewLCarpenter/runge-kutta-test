!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Kreiss problem (Dekker 7.5.1 pp. 214 (Kreiss' problem: Index 2)) 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
      module Kreiss_mod
      
      use precision_vars,    only: wp

      implicit none; save

      private
      public :: Kreiss

      real(wp), parameter :: tfinal_parameter = 1.0_wp     

      contains

      subroutine Kreiss(nveclen,neq,eps,dt,tfinal,iDT,time,rese_vec,resi_vec,akk)

      use precision_vars,    only: wp
      use control_variables, only: temporal_splitting,probname,xjac,var_names,&
     &                             tol,dt_error_tol,uvec,uexact,programstep

!-----------------------VARIABLES----------------------------------------------
      integer, parameter     :: vecl=2

      !INIT vars
      real(wp),                  intent(in   ) :: eps
      real(wp),                  intent(inout) :: dt
      integer,                   intent(  out) :: nveclen,neq
      real(wp),                  intent(  out) :: tfinal
      integer,                   intent(in   ) :: iDT
      real(wp),                  intent(in   ) :: time

      real(wp) :: epi,lambda_p,lambda_m,sin_t,cos_t, ustore,ep
      real(wp), dimension(vecl,vecl) :: E_mat,E_inv,wrk_matrix
      real(wp), dimension(vecl)         :: wrk_vec   
      
      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: rese_vec,resi_vec
      
      !Jacob vars
      real(wp), intent(in   ) :: akk
!------------------------------------------------------------------------------
      
      if(eps>(1.0_wp/3.0_wp)) then
        ep=1.0_wp/3.0_wp
      else
        ep=eps
      endif

      Program_Step_Select: select case(programstep)
        !**Pre-initialization. Get problem name and vector length**
        case('INITIALIZE_PROBLEM_INFORMATION')
          nvecLen = vecl
          neq = vecl
          probname='Kreiss   '   
          tol=5.0e-9_wp
          dt_error_tol=1.0e-10_wp
          
          allocate(var_names(neq))
          var_names(:)=(/'Differential', 'Algebraic   '/)
        
        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')
          !Initialize temporary constants
          epi =1.0_wp/ep
          cos_t=cos(time)
          sin_t=sin(time)
          
          !Time information
          choose_dt: select case(temporal_splitting)
          case('IMEX')    ; dt = 0.02_wp/10**((iDT-1)/30.0_wp) ! timestep
          case('IMPLICIT'); dt =  1.0_wp/10**((iDT-1)/20.0_wp) ! timestep
          case('FIRK')    ; dt =  1.0_wp/10**((iDT-1)/20.0_wp) ! timestep
          end select choose_dt
          tfinal=tfinal_parameter                    ! final time

          !**Exact Solution**
          lambda_p=0.5_wp*(-1.0_wp-epi) + 0.5_wp*sqrt((1.0_wp-epi)**2-4.0_wp)
          lambda_m=0.5_wp*(-1.0_wp-epi) - 0.5_wp*sqrt((1.0_wp-epi)**2-4.0_wp)
         
          E_mat(1,1)=cos(tfinal)
          E_mat(1,2)=-sin(tfinal)
          E_mat(2,1)=sin(tfinal)
          E_mat(2,2)=cos(tfinal)
          
          wrk_matrix(1,1)=ep
          wrk_matrix(1,2)=ep
          wrk_matrix(2,1)=1+lambda_p*ep
          wrk_matrix(2,2)=1+lambda_m*ep
          
          wrk_vec(1)=exp(lambda_p*tfinal)
          wrk_vec(2)=exp(lambda_m*tfinal)
          
          uexact=matmul(matmul(E_mat,wrk_matrix),wrk_vec)
          
          !**IC**
          uvec(1) = 2*ep
          uvec(2) = -1.0_wp*ep+1.0_wp
        
        case('BUILD_RHS') 
          !Initialize temporary constants
          cos_t=cos(time)
          sin_t=sin(time)
          choose_RHS_type: select case (Temporal_Splitting)
            case('IMEX') ! For IMEX schemes
              rese_vec(1) = dt*(-sin_t*sin_t*uvec(1)+sin_t*cos_t*uvec(2))
              rese_vec(2) = dt*( sin_t*cos_t*uvec(1)-cos_t*cos_t*uvec(2))
              resi_vec(1) = dt*(-cos_t*cos_t*uvec(1)-sin_t*cos_t*uvec(2))/ep
              resi_vec(2) = dt*(-sin_t*cos_t*uvec(1)-sin_t*sin_t*uvec(2))/ep
            case('IMPLICIT','FIRK') ! For fully implicit schemes
              rese_vec(:) = 0.0_wp
              resi_vec(1) =  dt*(-sin_t*sin_t*uvec(1)+sin_t*cos_t*uvec(2))+ &
     &                       dt*(-cos_t*cos_t*uvec(1)-sin_t*cos_t*uvec(2))/ep
              resi_vec(2) =  dt*( sin_t*cos_t*uvec(1)-cos_t*cos_t*uvec(2))+ &
     &                       dt*(-sin_t*cos_t*uvec(1)-sin_t*sin_t*uvec(2))/ep
          end select choose_RHS_type
          
        case('BUILD_JACOBIAN')
          cos_t=cos(time)
          sin_t=sin(time)
          choose_Jac_type: select case (Temporal_Splitting)
            case('IMEX') ! For IMEX schemes
              xjac(1,1) = 1.0_wp-akk*dt*(-cos_t*cos_t)/ep
              xjac(1,2) = 0.0_wp-akk*dt*(-sin_t*cos_t)/ep
              xjac(2,1) = 0.0_wp-akk*dt*(-sin_t*cos_t)/ep
              xjac(2,2) = 1.0_wp-akk*dt*(-sin_t*sin_t)/ep
            case('IMPLICIT') ! For fully implicit schemes
              xjac(1,1) = 1.0_wp-akk*dt*(-cos_t*cos_t-ep*sin_t*sin_t)/ep
              xjac(1,2) = 0.0_wp-akk*dt*(-cos_t*sin_t+ep*sin_t*cos_t)/ep
              xjac(2,1) = 0.0_wp-akk*dt*(-cos_t*sin_t+ep*cos_t*sin_t)/ep
              xjac(2,2) = 1.0_wp-akk*dt*(-sin_t*sin_t-ep*cos_t*cos_t)/ep
            case('FIRK') ! For fully implicit schemes
              xjac(1,1) = (-cos_t*cos_t-ep*sin_t*sin_t)/ep
              xjac(1,2) = (-cos_t*sin_t+ep*sin_t*cos_t)/ep
              xjac(2,1) = (-cos_t*sin_t+ep*cos_t*sin_t)/ep
              xjac(2,2) = (-sin_t*sin_t-ep*cos_t*cos_t)/ep
          end select choose_Jac_type
        
        case('ROTATE_VARS')
          tfinal=tfinal_parameter
          E_mat(1,1)= +cos(tfinal)
          E_mat(1,2)= -sin(tfinal)
          E_mat(2,1)= +sin(tfinal)
          E_mat(2,2)= +cos(tfinal)
          E_inv(:,:) = Transpose(E_mat) 
          uexact(:)  = MatMul(E_inv(:,:),uexact(:))
          uvec(:)    = MatMul(E_inv(:,:),uvec(:))

          ustore =uvec(1)
          uvec(1)=uvec(2)
          uvec(2)=ustore

          ustore   =uexact(1)
          uexact(1)=uexact(2)
          uexact(2)=ustore          

      end select Program_Step_select      
      end subroutine Kreiss
      end module Kreiss_mod
