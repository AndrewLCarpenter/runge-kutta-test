!******************************************************************************
! Module containing routines to describe Charney-Devore linear system in 
! H.E. De Swart, Acta Aplicandae Mathematicae 11 (1988) 49-96 ;
! Eqns are copine from D.T. Crommelin, A.J. Majda, J.Atmos.Sci,(2004):
! Strategies for Model reduction: comparing different optimal bases
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
! SBP_COEF_MODULE.F90           *DEFINES CSR OPERATORS 
! JACOBIAN_CSR_MOD.F90          *ALLOCATE AND STORE CSR JACOBIAN VARIABLES
! UNARY_MOD.F90                 *PERFORMS SPARSE MATRIX OPERATIONS
! MATVEC_MODULE.F90             *PERFORMS SPARSE MATRIX*VECTOR OPERATIONS
!******************************************************************************
      module Charney_DeVore6_mod

      use precision_vars,    only: wp, pi
      use control_variables, only: temporal_splitting,probname,xjac,var_names,&
     &                             tol,dt_error_tol,uvec,uexact,programstep

      implicit none; save 
  
      private
      public  :: Charney_DeVore6

!--------------------------------CONSTANTS-------------------------------------         
      real(wp), parameter :: sqrt2 = sqrt(2.0_wp)

!--------------------------------VARIABLES-------------------------------------         
      integer,  parameter :: neq  =   6 
      integer,  parameter :: vecl =   6

      real(wp), parameter :: x1s = +0.95000_wp
      real(wp), parameter :: x4s = -0.76095_wp
      real(wp), parameter :: CC  = +0.10_wp
      real(wp), parameter :: be  = +1.25_wp
      real(wp), parameter :: ga  = +0.20_wp
      real(wp), parameter :: b   = +0.50_wp

!     ============  Original definitions of coefficients  ===================
!     real(wp), parameter :: am = ( (8*sqrt2/pi) * (m*m/(4*m*m-1)) * &
!                                   ((b*b+m*m - 1)/(b*b+m*m)) )
!     real(wp), parameter :: bm = ( (be*b*b)/ (b*b+m*m) )
!     real(wp), parameter :: dm = ( (64.0_wp*sqrt2)/(15*pi) *   &
!                               &   (b*b-m*m + 1)/ (b*b+m*m) )
!     real(wp), parameter :: gm = ( ga * (4*m*m*m)/ (4*m*m-1) *      &
!                               &   (sqrt2*b/(pi * (b*b+m*m))) )
!     real(wp), parameter :: gms= ( ga * (4*m / (4*m*m-1)) * (sqrt2*b/pi) )
!     ============  Original definitions of coefficients  ===================

      real(wp), parameter :: ep = ( (16*sqrt2)/(5*pi) )

      real(wp), parameter :: a1 = ( (8*sqrt2/pi) * (1*1/(4*1*1-1)) * &
                                    ((b*b+1*1 - 1)/(b*b+1*1)) )
      real(wp), parameter :: a2 = ( (8*sqrt2/pi) * (2*2/(4*2*2-1)) * &
                                    ((b*b+2*2 - 1)/(b*b+2*2)) )

      real(wp), parameter :: b1 = ( (be*b*b)/ (b*b+1*1) )
      real(wp), parameter :: b2 = ( (be*b*b)/ (b*b+2*2) )

      real(wp), parameter :: d1 = ( (64.0_wp*sqrt2)/(15*pi) *   &
                                &   (b*b-1*1 + 1)/ (b*b+1*1) )
      real(wp), parameter :: d2 = ( (64.0_wp*sqrt2)/(15*pi) *   &
                                &   (b*b-2*2 + 1)/ (b*b+2*2) )

      real(wp), parameter :: g1 = ( ga * (4*1*1*1)/ (4*1*1-1) *      &
                                &   (sqrt2*b/(pi * (b*b+1*1))) )
      real(wp), parameter :: g2 = ( ga * (4*2*2*2)/ (4*2*2-1) *      &
                                &   (sqrt2*b/(pi * (b*b+2*2))) )

      real(wp), parameter :: g1s= ( ga * (4*1 / (4*1*1-1)) * (sqrt2*b/pi) )
      real(wp), parameter :: g2s= ( ga * (4*2 / (4*2*2-1)) * (sqrt2*b/pi) )
                             
      real(wp), dimension(neq,neq)  :: eye= reshape(              &
                                          & (/+1,+0,+0,+0,+0,+0,  &
                                          &   +0,+1,+0,+0,+0,+0,  &
                                          &   +0,+0,+1,+0,+0,+0,  &
                                          &   +0,+0,+0,+1,+0,+0,  &
                                          &   +0,+0,+0,+0,+1,+0,  &
                                          &   +0,+0,+0,+0,+0,+1/),&
                                          &      (/6,6/) ) 

      logical :: update_Jac

      contains    

!==============================================================================
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Charney_DeVore problem 
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
! exact_Charney_DeVore             -> Subroutine to build exact solution
! Charney_DeVore6_dudt              -> Subroutine to build dudt (LHS)
! Charney_DeVore6_Build_Jac -> Subroutine to build spatial Jacobian
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

      subroutine Charney_DeVore6(nveclen,eq,Cep,dt,tfinal,iDT,resE_vec,resI_vec,akk)

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
      
      !Jacob vars
      real(wp), intent(in   )  :: akk

!------------------------------------------------------------------------------   
      
      Program_Step_Select: select case(programstep)

        !**Pre-initialization. Get problem name, vector length and grid**

        case('INITIALIZE_PROBLEM_INFORMATION')
          nveclen = vecl
          eq      = vecl
          probname='Charney_DeVore6'     
          Jac_case='DENSE'
          tol=1.0e-14_wp  
          dt_error_tol=5.0e-14_wp
          
          allocate(var_names(neq))
          var_names(:)=(/ 'Differential', 'Differential',  &
                          'Differential', 'Differential',  &
                          'Differential', 'Differential'/)
          
        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')
          
          !Update=true for Jac/RHS updates
          update_Jac=.true. !reset every new epsilon/dt


          dt_max =   0001.00_wp                                   ! maximum dt
          tfinal =   4000.00_wp                                   ! final time   
          choose_dt: select case(temporal_splitting)
            case('EXPLICIT'); dt = dt_max/10**((iDT-1)/20.0_wp)   ! explicit timestep
            case default    ; dt = dt_max/10**((iDT-1)/20.0_wp)   ! implicit timestep      
          end select choose_dt
          
!         uvec(1) = +1.14_wp ;
!         uvec(2) = +0.00_wp ;
!         uvec(3) = +0.00_wp ;
!         uvec(4) = -0.91_wp ;
!         uvec(5) = +0.00_wp ;
!         uvec(6) = +0.00_wp ;
          uvec(1) = +0.14_wp ;
          uvec(2) = +2.00_wp ;
          uvec(3) = +0.00_wp ;
          uvec(4) = -1.91_wp ;
          uvec(5) = +1.00_wp ;
          uvec(6) = +0.00_wp ;

          call exact_Charney_DeVore6(uexact,Cep)                    !set exact solution at tfinal

        case('BUILD_RHS')

          call Charney_DeVore6_dUdt(uvec,Cep,dudt_E,dudt_I)        

          choose_RHS_type: select case (Temporal_Splitting)
            case('EXPLICIT') ! For fully explicit schemes
              resE_vec(:)= dt * (+ dudt_E(:) + dudt_I(:))
              resI_vec(:)= dt * 0.0_wp
            case('IMEX') ! For IMEX schemes
              resE_vec(:)= dt * (+ dudt_E(:)            ) 
              resI_vec(:)= dt * (            + dudt_I(:))
            case('IMPLICIT') ! For fully implicit schemes
              resE_vec(:)= dt * 0.0_wp
              resI_vec(:)= dt * (+ dudt_E(:) + dudt_I(:))
          end select choose_RHS_type
          
        case('BUILD_JACOBIAN')

          call Charney_DeVore6_Build_Jac(uvec,Cep,Jac_E,Jac_I)

          choose_Jac_type: select case(temporal_splitting)

            case('EXPLICIT') ! For fully implicit schemes

              xjac(:,:) = eye(:,:)

            case('IMEX') ! For IMEX schemes

              xjac(:,:) = eye(:,:) - akk*dt* (             + Jac_I(:,:))

            case('IMPLICIT') ! For fully implicit schemes

              xjac(:,:) = eye(:,:) - akk*dt* (+ Jac_E(:,:) + Jac_I(:,:))

          end select choose_Jac_type

      end select Program_Step_Select     

      end subroutine Charney_DeVore6

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

      subroutine exact_Charney_DeVore6(u,ep)

      real(wp),                  intent(in   ) :: ep
      real(wp), dimension(vecl), intent(  out) :: u

      integer                                  :: i
      real(wp), dimension(81,vecl+1)           :: ExactTot
      real(wp) :: diff

      u(:)=0.0_wp*ep
      
      !**Exact Solution** 
      open(unit=39,file='exact.Charney_DeVore6.data')
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

      end subroutine exact_Charney_DeVore6

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

      subroutine Charney_DeVore6_dUdt(u,Cep,dudt_E,dudt_I)

      real(wp), dimension(:), intent(in   ) :: u
      real(wp),               intent(in   ) :: Cep
      real(wp), dimension(:), intent(  out) :: dudt_E, dudt_I
           
      real(wp)                              :: x1,x2,x3,x4,x5,x6

        x1 = u(1); x2 = u(2); x3 = u(3); x4 = u(4); x5 = u(5); x6 = u(6);
  
        dudt_E(1) =                 + g1s*x3              
        dudt_E(2) = - (a1*x1-b1)*x3                        - d1*x4*x6
        dudt_E(3) = + (a1*x1-b1)*x2 - g1 *x1               + d1*x4*x5
        dudt_E(4) =                 + g2s*x6               + ep*(x2*x6 - x3*x5)
        dudt_E(5) = - (a2*x1-b2)*x6                        - d2*x3*x4
        dudt_E(6) = + (a2*x1-b2)*x5 - g2 *x4               + d2*x2*x4

        dudt_I(1) =                           - Cep*(x1-x1s)
        dudt_I(2) =                           - Cep* x2    
        dudt_I(3) =                           - Cep* x3    
        dudt_I(4) =                           - Cep*(x4-x4s)
        dudt_I(5) =                           - Cep* x5    
        dudt_I(6) =                           - Cep* x6   


      end subroutine Charney_DeVore6_dUdt

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

      subroutine Charney_DeVore6_Build_Jac(u,Cep,Jac_E,Jac_I)

      real(wp), dimension(:  ), intent(in   ) :: u
      real(wp),                 intent(in   ) :: Cep
      real(wp), dimension(:,:), intent(  out) :: Jac_E, Jac_I
      
      real(wp)                                :: x1,x2,x3,x4,x5,x6

        x1 = u(1); x2 = u(2); x3 = u(3); x4 = u(4); x5 = u(5); x6 = u(6);
  
        Jac_E(:,:) = 0.0_wp

        Jac_E(1,3) = g1s
        Jac_E(2,1) = -(a1*x3)
        Jac_E(2,3) = b1 - a1*x1
        Jac_E(2,4) = -(d1*x6)
        Jac_E(2,6) = -(d1*x4)
        Jac_E(3,1) = -g1 + a1*x2
        Jac_E(3,2) = -b1 + a1*x1
        Jac_E(3,4) = d1*x5
        Jac_E(3,5) = d1*x4
        Jac_E(4,2) = ep*x6
        Jac_E(4,3) = -(ep*x5)
        Jac_E(4,5) = -(ep*x3)
        Jac_E(4,6) = g2s + ep*x2
        Jac_E(5,1) = -(a2*x6)
        Jac_E(5,3) = -(d2*x4)
        Jac_E(5,4) = -(d2*x3)
        Jac_E(5,6) = b2 - a2*x1
        Jac_E(6,1) = a2*x5
        Jac_E(6,2) = d2*x4
        Jac_E(6,4) = -g2 + d2*x2
        Jac_E(6,5) = -b2 + a2*x1

        Jac_I(:,:) = 0.0_wp

        Jac_I(1,1) = -Cep
        Jac_I(2,2) = -Cep
        Jac_I(3,3) = -Cep
        Jac_I(4,4) = -Cep
        Jac_I(5,5) = -Cep
        Jac_I(6,6) = -Cep

      end subroutine Charney_DeVore6_Build_Jac      

!==============================================================================

      end module Charney_DeVore6_Mod
