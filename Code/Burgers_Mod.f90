!******************************************************************************
! Module containing routines to describe Burger's equation 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
! SBP_COEF_MODULE.F90           *DEFINES CSR OPERATORS 
! JACOBIAN_CSR_MOD.F90          *ALLOCATE AND STORE CSR JACOBIAN VARIABLES
! UNARY_MOD.F90                 *PERFORMS SPARSE MATRIX OPERATIONS
! MATVEC_MODULE.F90             *PERFORMS SPARSE MATRIX*VECTOR OPERATIONS
!******************************************************************************

      module Burgers_Module 

      use precision_vars, only: third, half, twothird, wp

      implicit none; save
  
      private
      public    ::  Burgers

!--------------------------------VARIABLES-------------------------------------         
      real(wp)  :: dx
      character(120),         parameter :: exact_solution = 'tanh'
!      character(120),         parameter :: exact_solution = 'Olver'
      real(wp),               parameter :: sig0 = -1.0_wp
      real(wp),               parameter :: sig1 = +1.0_wp
      
      integer,  parameter    :: vecl=128   
      integer,  parameter    :: neq=1
      real(wp), dimension(vecl) :: x      
      real(wp), parameter :: xL=0.0_wp,xR=1.0_wp

      real(wp), dimension(4), parameter :: d1vec0= (/-24.0_wp/17.0_wp,  &
                                                  &  +59.0_wp/34.0_wp,  &
                                                  &   -4.0_wp/17.0_wp,  &
                                                  &   -3.0_wp/34.0_wp/)
                                                                      
      real(wp), dimension(4), parameter :: d1vec1= (/+ 3.0_wp/34.0_wp,  &
                                                  &  + 4.0_wp/17.0_wp,  &
                                                  &  -59.0_wp/34.0_wp,  &
                                                  &  +24.0_wp/17.0_wp/)
!------------------------------------------------------------------------------      
      contains    
!==============================================================================
!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Burgers problem 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
! JACOBIAN_CSR_MOD.F90      *CONTAINS CSR JACOBIAN VARIABLES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From control_variables:
!   temporal_splitting -> string used in choose_RHS_type and choose_Jac_type,         character(len=80),                       not modified
!   probname           -> string used to set output file names,                       character(len=9),                            set
!   Jac_case           -> string used in to determine CSR Jacobian or dense Jacobian, character(len=6),                            set
!   tol                -> Newton iteration exit tolerance,                            real(wp),                                    set
!   dt_error_tol       -> L2 error tolerance for convergence plots,                   real(wp),                                    set
!   uvec               -> Array containing variables,                                 real(wp), dimension(u-vector length),        set & modified
!   uexact             -> Array containing exact solution to variables,               real(wp), dimension(u-vector length),        set
!   programstep        -> string used in program_step_select,                         character(len=80),                       not modified 
!   var_names          -> string array used to label output graphs,                   character(len=12), dimension(num. eq's),     set
! From SBP_Coef_Module:
!   Define_CSR_Operators -> Subroutine to define CSR derivative operators needed
!   D1                   -> First derivative operator, used to check if allocated, not modified
!   nnz_D2               -> Number of non zeros in 2nd derivative operator,        not modified
! From Jacobian_CSR_Mod:
!   Allocate_Jac_CSR_Storage -> Subroutine to create Jacobian and LU decomposition arrays for CSR problems
!
! MODULE VARIABLES/ROUTINES:
! vecl           -> total length of u-vector used for problem,      integer,                                                   not modified
! neq            -> number of equations in problem,                 integer,                                                   not modified
! x              -> stores x-grid points,                           real(wp),  dimension(u-vector length / num. eq's),         not modified
! dx             -> separation between grid points,                 real(wp),                                                  not modified
! grid           -> Subroutine to create grid
! exact_Burg     -> Subroutine to build exact solution
! Burgers_dudt   -> Subroutine to build dudt (LHS)
! Build_Jac      -> Subroutine to build Jacobian and store it
!******************************************************************************
! INPUTS:
! ep   -> Stiffness epsilon value,                                   real(wp)
! time -> Current solution time,                                     real(wp)
! iDT  -> Timestep counter from timestep loop to define dt,          integer
! akk  -> Diagonal term from RK scheme,                              real(wp)
! INOUTS:
! dt   -> timestep: created and then used later in the problem case, real(wp)
! OUTPUTS:
! nveclen  -> total length of u-vector used for problem              integer
! eq       -> number of equations in problem                         integer
! tfinal   -> final time for iterative solver                        real(wp)
! resE_vec -> explicit residual for a particular stage               real(wp), dimension(u-vector length)
! resI_vec -> implicit residual for a particular stage               real(wp), dimension(u-vector length)
!******************************************************************************
      subroutine Burgers(nveclen,eq,ep,dt,tfinal,iDT,time,resE_vec,resI_vec,akk)

      use control_variables, only: temporal_splitting,probname,Jac_case,     &
     &                             tol,dt_error_tol,uvec,uexact,programstep, &
     &                             var_names
      use SBP_Coef_Module,   only: Define_CSR_Operators,nnz_D2,D1
      use Jacobian_CSR_Mod,  only: Allocate_Jac_CSR_Storage

!-----------------------VARIABLES----------------------------------------------
      !INIT vars     
      integer,  intent(  out) :: nveclen,eq
      real(wp), intent(in   ) :: ep
      real(wp), intent(inout) :: dt
      real(wp), intent(  out) :: tfinal
      real(wp), intent(in   ) :: time
      integer,  intent(in   ) :: iDT

      real(wp)                :: tinitial
      
      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE_vec,resI_vec
      
      !Jacob vars
      real(wp), intent(in   ) :: akk
!------------------------------------------------------------------------------
      
      Program_Step_Select: select case(programstep)
        !**Pre-initialization. Get problem name, vector length and grid**
        case('INITIALIZE_PROBLEM_INFORMATION')
          nveclen = vecl
          eq = neq
          probname='Burgers  '     
          tol=1.0e-12_wp  
          dt_error_tol=5.0e-14_wp
          jac_case='SPARSE'
          
          allocate(var_names(neq))
          var_names(:)='Differential'
          
          call grid()           
 
        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')

          !Allocate derivative operators
          if(.not. allocated(D1)) call Define_CSR_Operators(vecl,dx)   

          !Time information
          tinitial=0.0_wp  ! initial time
          tfinal = 0.5_wp  ! final time    
          choose_dt: select case(temporal_splitting)       
            case('EXPLICIT')
              dt = 0.00005_wp*0.1_wp/10**((iDT-1)/20.0_wp) ! timestep  explicit 
            case('IMEX')
              print*,'Wrong case: IMEX not implemented.  Stopping'
              stop
            case('IMPLICIT')
              dt = 0.05_wp*0.1_wp/10**((iDT-1)/20.0_wp)    ! timestep  implcit
            case default
          end select choose_dt
          
          call exact_Burg(uvec,ep,tinitial) !set initial conditions   
          call exact_Burg(uexact,ep,tfinal) !set exact solution at tfinal
              
        case('BUILD_RHS')
          choose_RHS_type: select case (Temporal_Splitting)   
            case('EXPLICIT')
              call Burgers_dUdt(uvec,resE_vec,time,ep,dt)
              resI_vec(:)=0.0_wp
            case('IMPLICIT','FIRK') ! For fully implicit schemes
              call Burgers_dUdt(uvec,resI_vec,time,ep,dt)
              resE_vec(:)=0.0_wp
          end select choose_RHS_type
      !  PRINT*,'BUILDRHS'
        case('BUILD_JACOBIAN')
          choose_Jac_type: select case(temporal_splitting)
            case('IMPLICIT') ! For fully implicit schemes
              call Allocate_Jac_CSR_Storage(vecl,nnz_D2)
              call Build_Jac(uvec,ep,dt,akk,time)          
          end select choose_Jac_type
      end select Program_Step_select
      end subroutine Burgers
!==============================================================================
!==============================================================================
!==============================================================================
!******************************************************************************
! Subroutine to initialize x-grid
!******************************************************************************
! MODULE VARIABLES/ROUTINES:
! vecl -> total length of u-vector used for problem, integer,                                           not modified
! neq  -> number of equations in problem,            integer,                                           not modified
! x    -> stores x-grid points,                      real(wp),  dimension(u-vector length / num. eq's),     set
! dx   -> separation between grid points,            real(wp),                                              set
! xL   -> Left x-bound                               real(wp),                                          not modified
! xR   -> Right x-bound                              real(wp),                                          not modified
!******************************************************************************
      subroutine grid()
      integer :: i
      
      do i=1,vecl
        x(i)= xL + (xR-xL)*(i-1.0_wp)/(vecl-1.0_wp)
      enddo
      
      dx=x(2)-x(1)

      end subroutine grid
!==============================================================================
!******************************************************************************
! Subroutine to return 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
!
! MODULE VARIABLES/ROUTINES:
! vecl           -> total length of u-vector used for problem,      integer,                                                   not modified
! neq            -> number of equations in problem,                 integer,                                                   not modified
! x              -> stores x-grid points,                           real(wp),  dimension(u-vector length / num. eq's),         not modified
! NL_Burg_exactsolution -> Function to return solution at point     real(wp)
!******************************************************************************
! INPUTS:
! eps  -> Stiffness epsilon value, real(wp)
! time -> Current solution time,   real(wp)
! OUTPUTS:
! u    -> exact solution vector,   real(wp), dimension(u-vector length)
!******************************************************************************

      subroutine exact_Burg(u,eps,time)

      real(wp),                  intent(in   ) :: eps,time
      real(wp), dimension(vecl), intent(  out) :: u
      integer                                  :: i

      do i = 1,vecl
        u(i) =  NL_Burg_exactsolution(x(i),time,eps)
      enddo

      end subroutine exact_Burg
!==============================================================================
!******************************************************************************
! Subroutine to set RHS
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! SBP_COEF_MODULE.F90       *DEFINES CSR OPERATORS 
! MATVEC_MODULE.F90         *PERFORMS SPARSE MATRIX*VECTOR OPERATIONS
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp       -> working precision
!   third    -> exact 1/3, real(wp)
!   half     -> exact 1/2, real(wp)
!   twothird -> exact 2/3, real(wp)
! From SBP_Coef_Module:
!   Pinv -> Used for creating boundary conditions, real(wp), dimension(u-vector length),     not modified
!    D1  -> First derivative  a matrix operator,   real(wp), dimension(4*u-vector length),   not modified
!   jD1  -> First derivative ja matrix operator,   integer,  dimension(4*u-vector length),   not modified
!   iD1  -> First derivative ia matrix operator,   integer,  dimension(u-vector length + 1), not modified
!    D2  -> Second derivative  a matrix operator,  real(wp), dimension(4*u-vector length),   not modified
!   jD2  -> Second derivative ja matrix operator,  integer,  dimension(4*u-vector length),   not modified
!   iD2  -> Second derivative ia matrix operator,  integer,  dimension(u-vector length + 1), not modified
! From matvec_module:
!   amux -> multiplies a vector into a CSR matrix
!
! MODULE VARIABLES/ROUTINES:
! vecl   -> total length of u-vector used for problem, integer,                              not modified
! sig0   -> BC parameter,                              real(wp),                             not modified
! sig1   -> BC parameter,                              real(wp),                             not modified
! d1vec0 -> BC derivate parameter,                     real(wp), dimension(4),               not modified
! d1vec1 -> BC derivate parameter,                     real(wp), dimension(4),               not modified
! NL_Burg_exactsolution -> Function to return solution at point  real(wp)
!******************************************************************************
! INPUTS:
! u    -> variable vector,         real(wp), dimension(u-vector length)
! time -> Current solution time,   real(wp)
! eps  -> Stiffness epsilon value, real(wp)
! dt   -> timestep,                real(wp)
! OUTPUTS:
! dudt  -> RHS vector,             real(wp), dimension(u-vector length)
!******************************************************************************

      subroutine Burgers_dUdt(u,dudt,time,eps,dt)
      
      use SBP_Coef_Module, only: Pinv,D1,D2,jD1,jD2,iD1,iD2 
      use matvec_module,   only: amux

      real(wp), dimension(vecl), intent(in   ) ::  u
      real(wp), dimension(vecl), intent(  out) :: dudt
      real(wp),                  intent(in   ) :: time, eps, dt

      real(wp), dimension(vecl)                :: f, df, dfv, gsat, wrk    
      real(wp)                                 :: uL,uR,du0,du1,a0,a1,g0,g1
      
      gsat = 0.0_wp ; df  = 0.0_wp ; dfv = 0.0_wp ;

!--------------------------FLUX------------------------------------------------
      !     Canonical splitting for Inviscid flux in Burgers eqn
      !     f(u)_x = 2/3 (u u/2)_x + 1/3 u u_x

      f(:) = half*u(:)*u(:)
      call amux(vecl,f,wrk,D1,jD1,iD1)

      df(:) = df(:) + twothird*wrk(:)

      call amux(vecl,u,wrk,D1,jD1,iD1)
      df(:) = df(:) + third*u(:)*wrk(:)

      !  Viscous Flux in Burgers eqn

      call amux(vecl,u,dfv,D2,jD2,iD2)

!--------------Inflow Boundary Condition---------------------------------------

      uR = u(1)

      a0 = third*(uR + sqrt(uR**2+1.0e-28_wp))    !a0 = third*(uR + abs(uR))
      g0 = a0 * NL_Burg_exactsolution(x(1),time,eps)                &
         - eps* NL_Burg_exact_derivative(x(1),time,eps)

      du0 = dot_product(d1vec0(1:4),u(1:4)) / dx

      gsat(1) = gsat(1) + sig0 * Pinv(1) * (a0*uR - eps*du0 - g0) 

!---------------Outflow Boundary Condition-------------------------------------

      uL = u(vecl)

      a1 = third*(uL - sqrt(uL**2+1.0e-28_wp)) !a1 = third*(uL - abs(uL))
      g1 = a1 * NL_Burg_exactsolution(x(vecl),time,eps)                &
         - eps* NL_Burg_exact_derivative(x(vecl),time,eps)
      du1 = dot_product(d1vec1(1:4),u(vecl-3:vecl)) / dx

      gsat(vecl) = gsat(vecl) + sig1 * Pinv(vecl) * (a1*uL-eps*du1-g1)

!--------------------Sum all terms---------------------------------------------

      dudt(:) = dt*(eps*dfv(:) - df(:) + gsat(:))

      end subroutine Burgers_dUdt
      
!==============================================================================
!******************************************************************************
! Subroutine to finish building Jacobian and set it to global variables
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! SBP_COEF_MODULE.F90       *DEFINES CSR OPERATORS 
! MATVEC_MODULE.F90         *PERFORMS SPARSE MATRIX*VECTOR OPERATIONS
! UNARY_MOD.F90             *PERFORMS SPARSE MATRIX OPERATIONS
! JACOBIAN_CSR_MOD.F90      *ALLOCATE AND STORE CSR JACOBIAN VARIABLES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
!   third    -> exact 1/3, real(wp)
! From SBP_Coef_Module:
!   Pinv -> Used for creating boundary conditions, real(wp), dimension(u-vector length),     not modified
!    D1  -> First derivative  a matrix operator,   real(wp), dimension(4*u-vector length),   not modified
!   jD1  -> First derivative ja matrix operator,   integer,  dimension(4*u-vector length),   not modified
!   iD1  -> First derivative ia matrix operator,   integer,  dimension(u-vector length + 1), not modified
!    D2  -> Second derivative  a matrix operator,  real(wp), dimension(4*u-vector length),   not modified
!   jD2  -> Second derivative ja matrix operator,  integer,  dimension(4*u-vector length),   not modified
!   iD2  -> Second derivative ia matrix operator,  integer,  dimension(u-vector length + 1), not modified
! From matvec_module:
!   amux -> multiplies a vector into a CSR matrix
! From unary_mod:
!   aplb   -> Subroutine to add two CSR matricies
!   aplsca -> Subroutine to add diagonal constant to a CSR matrix
!   amudia -> Subroutine to add 
! From Jacobian_CSR_Mod:
!   iaJac -> ia matrix for global storage of Jacobian, integer,  dimension(u-vector length + 1),                               set
!   jaJac -> ja matrix for global storage of Jacobian, integer,  dimension(dependant on temporal_splitting, see main routine), set
!    aJac ->  a matrix for global storage of Jacobian, real(wp), dimension(dependant on temporal_splitting, see main routine), set
!
! MODULE VARIABLES/ROUTINES:
! vecl   -> total length of u-vector used for problem, integer, not modified
! sig0   -> BC parameter,                              real(wp),                             not modified
! sig1   -> BC parameter,                              real(wp),                             not modified
! d1vec0 -> BC derivate parameter,                     real(wp), dimension(4),               not modified
! d1vec1 -> BC derivate parameter,                     real(wp), dimension(4),               not modified
! NL_Burg_exactsolution -> Function to return solution at point  real(wp)
!******************************************************************************
! INPUTS:
! u    -> variable vector,         real(wp), dimension(u-vector length)
! eps  -> Stiffness epsilon value, real(wp)
! dt  -> timestep,                     real(wp)
! akk -> Diagonal term from RK scheme, real(wp)
! time -> Current solution time,   real(wp)
!******************************************************************************

      subroutine Build_Jac(u,eps,dt,akk,time)

      use SBP_Coef_Module,  only: Pinv,D1,D2,jD1,jD2,iD1,iD2,nnz_D2       
      use matvec_module,    only: amux
      use unary_mod,        only: aplb,aplsca,amudia,diamua,apldia
      use Jacobian_CSR_Mod, only: iaJac,jaJac,aJac
      
      real(wp), dimension(vecl), intent(in) :: u
      real(wp),                  intent(in) :: eps,dt,akk,time    
        
      integer,  dimension(vecl+1) :: iwrk1,iwrk2,iwrk3,iJac
      integer,  dimension(nnz_D2) :: jwrk1,jwrk2,jwrk3,jJac
      real(wp), dimension(nnz_D2) :: wrk1,wrk2,wrk3,Jac,wrk4,eps_d2      
      
      integer,  dimension(vecl)  :: iw
      real(wp), dimension(vecl)  :: diag
      integer,  dimension(2)     :: ierr
      real(wp)                   :: uL,a1_d,uR,a0_d,a0,a1
     
!---------------------dgsat/dt-------------------------------------------------
      uL=u(vecl)
      a1 = third*(uL - sqrt(uL**2+1.0e-28_wp))
      a1_d=third*(1-uL/sqrt(uL**2+1.0e-28_wp))
   
      uR=u(1)
      a0 = third*(uR + sqrt(uR**2+1.0e-28_wp))
      a0_d=third*(1+uR/sqrt(uR**2+1.0e-28_wp))

!******Steps to make xjac******
!       jac=eps*d2-2/3*d1*u-1/3*u*d1-1/3*diag(d1*u)+dgsat/du
!1      wrk1=-2/3*d1*u
!2      wrk2=-1/3*u*d1
!3      wrk3=wrk1+wrk2
!4      eps_d2=eps*d2
!5      jac=eps_d2+wrk3-1/3*diag(d1*u)+dgsat/du
!6      wrk4=-akk*dt*jac
!
!       xjac=I-akk*dt*jac=-akk*dt*jac+(1)*I
!7      xjac=wrk4+(1)*I
! 
!8      xjac-> wrk4, iJac, jJac
!******************************
! --------------Make xjac -----------------------------------------------------
      call amudia(vecl,1,D1,jD1,iD1,u,wrk1,jwrk1,iwrk1)                   !1
      wrk1(:)=-twothird*wrk1(:)                                              !1
      
      call diamua(vecl,1,D1,jD1,iD1,u,wrk2,jwrk2,iwrk2)                   !2
      wrk2(:)=-third*wrk2(:)                                                 !2
      
      call aplb(vecl,vecl,1,wrk1,jwrk1,iwrk1,wrk2,jwrk2,iwrk2,wrk3,&   !3
     &          jwrk3,iwrk3,nnz_D2,iw,ierr(1))                               !3

      eps_d2=eps*D2(:)                                                       !4
      call aplb(vecl,vecl,1,eps_d2,jD2,iD2,wrk3,jwrk3,iwrk3,Jac,&     !5a
     &          jJac,iJac,nnz_D2,iw,ierr(2))                                !5a

      call amux(vecl,-third*u,diag,D1,jD1,iD1) !First-derivative of u vec!5b
      call apldia(vecl,0,Jac,jJac,iJac,diag,Jac,jJac,iJac,iw)            !5b

     ! R - 5c
      Jac(1)=Jac(1)+sig0*Pinv(1)*(a0+a0_d*(uR- &
     &       NL_Burg_exactsolution(x(1),time,eps))) 
      Jac(1:4)=Jac(1:4)+sig0*Pinv(1)*(-eps)*d1vec0(:)/dx
      
      ! L - 5c
      Jac(nnz_D2)=Jac(nnz_D2)+sig1*Pinv(vecl)*(a1+a1_d*(uL- &
     &            NL_Burg_exactsolution(x(vecl),time,eps))) 
      Jac(nnz_D2-3:nnz_D2)=Jac(nnz_D2-3:nnz_D2)+ &
     &                     sig1*Pinv(vecl)*(-eps)*d1vec1(:)/dx
     
      wrk4=-akk*dt*Jac                                                       !6
      call aplsca(vecl,wrk4,jJac,iJac,1.0_wp,iw)                          !7

      !8
      iaJac=iJac
      jaJac=jJac
      aJac=wrk4
!-----------------------------------------------------------------------------  
      if (sum(ierr)/=0) then ! Catch errors
        print*,'Error building Jacobian'
        stop
      endif
   
      end subroutine Build_Jac
      
!==============================================================================
! 
!!      subroutine error(u,time,eps,dt)
!!
!!      use SBP_Coef_Module, only: Pmat
!!
!!      real(wp),                     intent(in   ) :: time, eps, dt
!!      real(wp), dimension(vecl), intent(inout) :: u
!!
!!      integer                                    :: i
!!      real(wp)                                   :: errL2,errLinf,wrk,psum
!!
!!!     calculate the rms residual over the domain                        
!!
!!      psum    = 0.0_wp
!!      errL2   = 0.0_wp
!!      errLinf = 0.0_wp
!!      do i=1,vecl
!!              wrk = abs(u(i) - NL_Burg_exactsolution(x(i),time,eps))
!!            errL2 =  errL2 +  pmat(i)*wrk * wrk
!!          errLinf =  max(errLinf,wrk)
!!             psum = psum + pmat(i)
!!      enddo
!!      errL2  = ( sqrt(errL2/vecl/psum) )
!!      errLinf= ( errLinf)
!!!!      write( *,89)ixd-1,errL2,errLinf
!!   89 format('P=',I4,' ,,  L2',e18.10,' ,,  Linf',e18.10,' ,,')
!!
!!      end subroutine error
!==============================================================================
!
!!      subroutine plot(punit,u,time,eps)
!!
!!      integer,                      intent(in   ) :: punit
!!      real(wp),                     intent(in   ) :: time, eps
!!      real(wp), dimension(vecl), intent(inout) :: u
!!
!!      integer                                    :: i
!!      real(wp)                                   :: wrk,err
!!
!!!     write to plotter file                                             
!!      do i=1,vecl
!!        wrk = NL_Burg_exactsolution(x(i),time,eps)
!!        err = ( abs( wrk - u(i) ) + 1.0e-15 )
!!        write(punit,2)x(i),wrk,u(i),err
!!      enddo
!!    2 format(4(1x,e15.7))
!!
!!      end subroutine plot
!==============================================================================
!******************************************************************************
! Function to return exact solution
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp -> working precision
!
! MODULE VARIABLES/ROUTINES:
! exact_soluton -> String to determine which exact solution 
!******************************************************************************
! INPUTS:
! xin -> x value of point,        real(wp)
! tin -> current time,            real(wp)
! eps -> Stiffness epsilon value, real(wp)
! OUTPUTS:
! NL_Burg_exactsolution -> Function to return solution at point  real(wp)
!******************************************************************************
      function NL_Burg_exactsolution(xin,tin,eps)

      real(wp), intent(in) :: xin, tin, eps
      real(wp), parameter  ::  a = -0.40_wp, b = 1.0_wp , d = +0.50_wp
      real(wp)             ::  c = (a+b)/2

      real(wp) :: NL_Burg_exactsolution
      real(wp) :: t1,t2, x0

      select case (exact_solution)
        case('Olver')         !          pp. 1190   Peter J. Olver

          t1  = (xin - c*tin - d )
          t2  = exp((b-a)*t1/(2.0_wp*eps))

          NL_Burg_exactsolution = (a * t2  + b ) / (1.0_wp * t2 + 1.0_wp)

        case('tanh')

          x0 =  0.1_wp
          NL_Burg_exactsolution = 1.1_wp - tanh(0.5_wp*(xin-1.1_wp*tin-x0)/eps)

      end select

      end function NL_Burg_exactsolution
!==============================================================================
!******************************************************************************
! Function to return exact derivative
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp -> working precision
!
! MODULE VARIABLES/ROUTINES:
! exact_soluton -> String to determine which exact solution 
!******************************************************************************
! INPUTS:
! xin -> x value of point,        real(wp)
! tin -> current time,            real(wp)
! eps -> Stiffness epsilon value, real(wp)
! OUTPUTS:
! NL_Burg_exact_derivative -> Function to return derivative at point, real(wp)
!******************************************************************************
      function NL_Burg_exact_derivative(xin,tin,eps)

       real(wp), intent(in) :: xin, tin, eps
       real(wp), parameter  ::  a = -0.40_wp, b = 1.0_wp , d = +0.50_wp
       real(wp)             ::  c = (a+b)/2

       real(wp) :: NL_Burg_exact_derivative
       real(wp) :: t1,t2,t3, x0

       select case (exact_solution)
         case('Olver')         !          pp. 1190   Peter J. Olver

           t1  = (xin - c*tin - d )
           t2  = ((a - b)*t1)/(4.0_wp*eps)
           t3  = 2.0_wp / (exp(t2)+exp(-t2))

           NL_Burg_exact_derivative =  -((a - b)**2* t3**2) /(8.0_wp*eps)

         case('tanh')

           x0 =  0.1_wp
           NL_Burg_exact_derivative =  -1.0_wp/cosh((-1.1_wp*tin+xin-x0)/ &
     &                                 (2.0_wp*eps))**2 / (2.0_wp * eps)

       end select

      end function NL_Burg_exact_derivative
!==============================================================================
      end module Burgers_Module
