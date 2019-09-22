!******************************************************************************
! Module containing routines to describe Boldrighini-Franceschini linear system in 
! Commun. Math. PHy. 64, 159-170 (1979) ;
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
! SBP_COEF_MODULE.F90           *DEFINES CSR OPERATORS 
! JACOBIAN_CSR_MOD.F90          *ALLOCATE AND STORE CSR JACOBIAN VARIABLES
! UNARY_MOD.F90                 *PERFORMS SPARSE MATRIX OPERATIONS
! MATVEC_MODULE.F90             *PERFORMS SPARSE MATRIX*VECTOR OPERATIONS
!******************************************************************************
      module Boldrighini_Franceschini_mod

      use precision_vars,    only: wp, eyeN
      use control_variables, only: temporal_splitting,probname,xjac,var_names,&
     &                             tol,dt_error_tol,uvec,uexact,programstep

      implicit none; save 
  
      private
      public  :: Boldrighini_Franceschini

!--------------------------------VARIABLES-------------------------------------         
      integer,  parameter :: neq  =   5 
      integer,  parameter :: vecl =   5
!  R01 =  5   sqrt(3/2) ~  6.123724357
!  R02 = 80/9 sqrt(3/2) ~ 10.886621079
!  R03 =  22.8537016  (soft transition)
!  R04 =  28.41
!  R05 =  28.64  
!  R06 =  28.6660
!  R07 =  28.6662
!  R08 = ~28.700, 
!  R09 = ~28.716, 
!  R10 = ~28.719, 
!  R11 = ~28.720

      real(wp), parameter ::   Re = + 28.680_wp
      real(wp), parameter ::   a1 = + 4.0_wp
      real(wp), parameter ::   a2 = + 4.0_wp
      real(wp), parameter ::   b1 = + 3.0_wp
      real(wp), parameter ::   c1 = + 7.0_wp
      real(wp), parameter ::   d1 = + 1.0_wp
      real(wp), parameter ::   e1 = + 3.0_wp

      logical :: update_Jac

      contains    

!==============================================================================
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Boldrighini_Franceschini problem 
!==============================================================================
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
! SBP_COEF_MODULE.F90       *DEFINES CSR OPERATORS 
! UNARY_MOD.F90             *PERFORMS SPARSE MATRIX OPERATIONS
! JACOBIAN_CSR_MOD.F90      *ALLOCATE AND STORE CSR JACOBIAN VARIABLES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
!   pi  -> mathematical Pi,  real(wp)
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
!  Define_CSR_Operators -> Subroutine to define CSR derivative operators needed
!  D1_per               -> First derivative operator, used to check if allocated, not modified
! From unary_mod:
!   aplb -> Subroutine to add one CSR matrix to another
! From Jacobian_CSR_Mod:
!   Allocate_Jac_CSR_Storage -> Subroutine to create Jacobian and LU decomposition arrays for CSR problems
!
! MODULE VARIABLES/ROUTINES:
! neq           -> number of equations in problem,                 integer,                                                             not modified
! vecl          -> total length of u-vector used for problem,      integer,                                                             not modified
! update_Jac    -> logical flag set to decide when to update Jac,  logical,                                                                 set & modified
! exact_Boldrighini_Franceschini             -> Subroutine to build exact solution
! Boldrighini_Franceschini_dudt              -> Subroutine to build dudt (LHS)
! Boldrighini_Franceschini_Build_Jac -> Subroutine to build spatial Jacobian
!******************************************************************************
! INPUTS:
! ep   -> Stiffness epsilon value,                                   real(wp)
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

      subroutine Boldrighini_Franceschini(nveclen,eq,Cep,dt,tfinal,iDT,resE_vec,resI_vec,akk)

      use control_variables, only: temporal_splitting,probname,Jac_case,     &
     &                             tol,dt_error_tol,uvec,uexact,programstep, &
     &                             var_names

!-----------------------VARIABLES----------------------------------------------
      !INIT vars     
      real(wp), intent(in   ) :: Cep
      integer,  intent(in   ) :: iDT
      integer,  intent(  out) :: eq

      real(wp), intent(inout) :: dt

      integer,  intent(  out) :: nveclen
      real(wp), intent(  out) :: tfinal

      real(wp)                :: dt_max
      
      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE_vec,resI_vec
      real(wp), dimension(vecl)                :: dudt_E,dudt_I
      real(wp), dimension(vecl,vecl)           :: Jac_E, Jac_I
      real(wp), dimension(vecl)                :: ran
      
      !Jacob vars
      real(wp), intent(in   )  :: akk

!------------------------------------------------------------------------------   
      
      Program_Step_Select: select case(programstep)

        !**Pre-initialization. Get problem name, vector length and grid**

        case('INITIALIZE_PROBLEM_INFORMATION')
          nveclen     = vecl
          eq          = vecl
          probname    = 'Boldrighini_Franceschini'     
          Jac_case    = 'DENSE'
          tol         = 1.0e-13_wp  
          dt_error_tol= 5.0e-14_wp
          
          allocate(var_names(neq))
          var_names(:)=(/ 'Differential', 'Differential',  &
                          'Differential', 'Differential',  &
                          'Differential'/)
          
        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')
          
          !Update=true for Jac/RHS updates
          update_Jac=.true. !reset every new epsilon/dt


          dt_max =   0000.010_wp                                  ! maximum dt
          tfinal =   0100.00_wp                                   ! final time   
          choose_dt: select case(temporal_splitting)
            case('EXPLICIT'); dt = dt_max/10**((iDT-1)/20.0_wp)   ! explicit timestep
            case default    ; dt = dt_max/10**((iDT-1)/20.0_wp)   ! implicit timestep      
          end select choose_dt
          
          call random_number(ran(:)) ;
!         uvec(:) = 2*(-0.5+ran(:))
          uvec(1) = +1.0_wp
          uvec(2) = -2.0_wp
          uvec(3) = +3.0_wp
          uvec(4) = -4.0_wp
          uvec(5) = +5.0_wp

          call exact_Boldrighini_Franceschini(uexact,Cep)         !set exact solution at tfinal

        case('BUILD_RHS')

          call Boldrighini_Franceschini_dUdt(uvec,Cep,dudt_E,dudt_I)

          choose_RHS_type: select case (Temporal_Splitting)
            case('EXPLICIT') ! For fully explicit schemes
              resE_vec(:)= dt * (+ dudt_E(:) + dudt_I(:))
              resI_vec(:)= dt * 0.0_wp
            case('IMEX') ! For IMEX schemes
              resE_vec(:)= dt * (+ dudt_E(:)            ) 
              resI_vec(:)= dt * (            + dudt_I(:))
            case('IMPLICIT','FIRK') ! For fully implicit schemes
              resE_vec(:)= dt * 0.0_wp
              resI_vec(:)= dt * (+ dudt_E(:) + dudt_I(:))
          end select choose_RHS_type
          
        case('BUILD_JACOBIAN')

          call Boldrighini_Franceschini_Build_Jac(uvec,Cep,Jac_E,Jac_I)

          choose_Jac_type: select case(temporal_splitting)

            case('EXPLICIT') ! For fully implicit schemes

              xjac(:,:) = eyeN(vecl)

            case('IMEX') ! For IMEX schemes

              xjac(:,:) = eyeN(vecl) - akk*dt* (             + Jac_I(:,:))

            case('IMPLICIT') ! For fully implicit schemes

              xjac(:,:) = eyeN(vecl) - akk*dt* (+ Jac_E(:,:) + Jac_I(:,:))

            case('FIRK') ! For fully implicit schemes

              xjac(:,:) = + Jac_E(:,:) + Jac_I(:,:)

          end select choose_Jac_type

      end select Program_Step_Select     

      end subroutine Boldrighini_Franceschini

!==============================================================================
! Subroutine to return exact solution vector
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp -> working precision
!
! MODULE VARIABLES/ROUTINES:
! vecl -> total length of u-vector used for problem, integer, not modified
!******************************************************************************
! INPUTS:
! eps -> Stiffness epsilon value, real(wp)
! OUTPUTS:
! u   -> exact solution vector,   real(wp), dimension(u-vector length)
!******************************************************************************

      subroutine exact_Boldrighini_Franceschini(u,ep)

      real(wp),                  intent(in   ) :: ep
      real(wp), dimension(vecl), intent(  out) :: u

      integer                                  :: i
      real(wp), dimension(81,vecl+1)           :: ExactTot
      real(wp) :: diff

      u(:)=0.0_wp*ep
      
      !**Exact Solution** 
      open(unit=39,file=&
                &  './Exact_Data/exact.Boldrighini_Franceschini.data')
      rewind(39)
      do i=1,81
        read(39,*)ExactTot(i,1:vecl)
        ExactTot(i,vecl+1) = 1.0_wp/10**((i-1)/(10.0_wp))  !  used for 81 values of ep
      enddo
      do i=1,81
        diff = abs(ExactTot(i,vecl+1) - ep)
        if(diff <= 1.0e-10_wp)then
          u(:) = ExactTot(i,:vecl)
          exit
        endif
      enddo

      end subroutine exact_Boldrighini_Franceschini

!==============================================================================
! Subroutine to set RHS
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! SBP_COEF_MODULE.F90       *DEFINES CSR OPERATORS 
! MATVEC_MODULE.F90         *PERFORMS SPARSE MATRIX*VECTOR OPERATIONS
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
!
! MODULE VARIABLES/ROUTINES:
! vecl -> total length of u-vector used for problem, integer, not modified
!******************************************************************************
! INPUTS:
! u   -> variable vector,         real(wp), dimension(u-vector length)
! Cep -> Stiffness epsilon value, real(wp)
! OUTPUTS:
! dudt        -> RHS vector involving just spacial derivatives, real(wp), dimension(u-vector length)
!******************************************************************************

      subroutine Boldrighini_Franceschini_dUdt(u,Cep,dudt_E,dudt_I)

      real(wp), dimension(:), intent(in   ) :: u
      real(wp),               intent(in   ) :: Cep
      real(wp), dimension(:), intent(  out) :: dudt_E, dudt_I
           
      real(wp)                              :: x1,x2,x3,x4,x5

        x1 = u(1); x2 = u(2); x3 = u(3); x4 = u(4); x5 = u(5);
  
        dudt_E(:) = 0.0_wp
        dudt_E(3) = Re

        dudt_I(1) = -2.0_wp*x1 + a1*x2*x3 + a2*x4*x5
        dudt_I(2) = -9.0_wp*x2 + b1*x1*x3
        dudt_I(3) = -5.0_wp*x3 - c1*x1*x2            ! + Re
        dudt_I(4) = -5.0_wp*x4 - d1*x1*x5
        dudt_I(5) = -1.0_wp*x5 - e1*x1*x4


      end subroutine Boldrighini_Franceschini_dUdt

!==============================================================================
! Subroutine to set spatial Jacobian
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! MODULE VARIABLES/ROUTINES
! vecl -> total length of u-vector used for problem, integer, not modified
!******************************************************************************
! OUTPUTS:
! a  ->  a combined derivative matrix for output to main routine, real(wp), dimension(neq,neq)
!******************************************************************************

      subroutine Boldrighini_Franceschini_Build_Jac(u,Cep,Jac_E,Jac_I)

      real(wp), dimension(:  ), intent(in   ) :: u
      real(wp),                 intent(in   ) :: Cep
      real(wp), dimension(:,:), intent(  out) :: Jac_E, Jac_I
      
      real(wp)                                :: x1,x2,x3,x4,x5

        x1 = u(1); x2 = u(2); x3 = u(3); x4 = u(4); x5 = u(5);
  
        Jac_E(:,:) = 0.0_wp

        Jac_I(:,:) = 0.0_wp

        Jac_I(1,1) = -2.0_wp
        Jac_I(1,2) = +a1*x3
        Jac_I(1,3) = +a1*x2
        Jac_I(1,4) = +a2*x5
        Jac_I(1,5) = +a2*x4
        Jac_I(2,1) = +b1*x3
        Jac_I(2,2) = -9.0_wp
        Jac_I(2,3) = +b1*x1
        Jac_I(2,4) = +0.0_wp
        Jac_I(2,5) = +0.0_wp
        Jac_I(3,1) = -c1*x2
        Jac_I(3,2) = -c1*x1
        Jac_I(3,3) = -5.0_wp
        Jac_I(3,4) = +0.0_wp
        Jac_I(3,5) = +0.0_wp
        Jac_I(4,1) = -d1*x5
        Jac_I(4,2) = +0.0_wp
        Jac_I(4,3) = +0.0_wp
        Jac_I(4,4) = -5.0_wp
        Jac_I(4,5) = -d1*x1
        Jac_I(5,1) = -e1*x4
        Jac_I(5,2) = +0.0_wp
        Jac_I(5,3) = +0.0_wp
        Jac_I(5,4) = -e1*x1
        Jac_I(5,5) = -1.0_wp

      end subroutine Boldrighini_Franceschini_Build_Jac      

!==============================================================================

      end module Boldrighini_Franceschini_Mod
