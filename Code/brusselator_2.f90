!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate Brusselator
!
! x'=a+x*x*y-b*x-x
! y'=b*x-x*x*y
!
! (du/dt)=(df/dx df/dy)*(u)
! (dv/dt) (dg/dx dg/dy) (v)
!
! J=(b-1 a)
!   (-b -a)
!
! u'=1-4*u+u*u*v+ep*d^2u/dx^2
! v'=  3*u-u*u*v+ep*d^2v/dx^2
! where:
! u(x,t), v(x,t)
! u(0,t)=u(1,t)=1       v(0,t)=v(1,t)=3
! u(x,0)=1+sin(2*pi*x)  v(x,0)=3
! 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
      module Brusselator_mod
      
      use precision_vars,    only: wp,two,pi
      
      implicit none; save            
      
      private
      public :: Brusselator
      
      integer,  parameter :: vecl=32     ! total uvec length == 2*vecl
      integer,  parameter :: neq=2
      real(wp), parameter :: xL=0.0_wp, xR=1.0_wp
      
      real(wp), dimension(vecl) :: x
      real(wp)                  :: dx
      
      integer,  dimension(10*vecl)  :: jDeriv2_p
      real(wp), dimension(10*vecl)  :: Deriv2_p
      integer,  dimension(vecl*2+1) :: iDeriv2_p
      
      contains
      
!==============================================================================      
!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Brusselator problem 
!******************************************************************************
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
!   two -> exact 2.0_wp,     real(wp)
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
!   Define_CSR_Operators -> Subroutine to define CSR derivative operators needed
!   D2_per               -> Second derivative operator, used to check if allocated,          not modified
!   nnz_D2_per           -> Second derivative operator number of non zeros,         integer, not modified
! From unary_mod:
!   aplb -> Subroutine to add one CSR matrix to another
! From Jacobian_CSR_Mod:
!   Allocate_Jac_CSR_Storage -> Subroutine to create Jacobian and LU decomposition arrays for CSR problems
!
! MODULE VARIABLES/ROUTINES:
! vecl      -> total length of u-vector used for problem,      integer,                                                   not modified
! neq       -> number of equations in problem,                 integer,                                                   not modified
! x         -> stores x-grid points,                           real(wp),  dimension(u-vector length / num. eq's),         not modified
! dx        -> separation between grid points,                 real(wp),                                                  not modified
! jDeriv2_p -> permuted ja matrix for D1 operators             integer,   dimension(2 * 5 * u-vector length / num. eq's), not modified
! Deriv2_p  -> permuted  a matrix for D1 operators             real(wp),  dimension(2 * 5 * u-vector length / num. eq's), not modified
! iDeriv2_p -> permuted ia matrix for D1 operators             integer,   dimension(u-vector length + 1),                 not modified
! grid              -> Subroutine to create grid
! exact_Bruss       -> Function to build exact solution
! Bruss_dudt        -> Subroutine to build dudt (LHS)
! Build_Spatial_Jac -> Subroutine to build spatial Jacobian
! Build_Source_Jac  -> Subroutine to build source  Jacobian
! Build_Jac         -> Subroutine to build Jacobian and store it
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
      subroutine Brusselator(nveclen,eq,ep,dt,tfinal,iDT,resE_vec,resI_vec,akk)

      use control_variables, only: temporal_splitting,probname,Jac_case,     &
     &                             tol,dt_error_tol,uvec,uexact,programstep, &
     &                             var_names
      use SBP_Coef_Module,   only: Define_CSR_Operators,D2_per,nnz_D2_per
      use unary_mod,         only: aplb
      use Jacobian_CSR_Mod,  only: Allocate_Jac_CSR_Storage
!-----------------------VARIABLES----------------------------------------------

      !INIT vars
      real(wp), intent(in   ) :: ep
      real(wp), intent(inout) :: dt
      integer,  intent(  out) :: nveclen,eq
      real(wp), intent(  out) :: tfinal
      integer,  intent(in   ) :: iDT

      !RHS vars
      real(wp), dimension(vecl*2), intent(  out) :: resE_vec,resI_vec
      real(wp), dimension(vecl*2)                :: dudt_Dx,dudt_Source
            
      !Jacob vars
      real(wp), intent(in   )       :: akk
      integer                       :: nnz_Jac
      real(wp), dimension(4*vecl)   :: Source      
      integer,  dimension(4*vecl)   :: jSource
      integer,  dimension(2*vecl+1) :: iSource,iwrk_Jac
      real(wp), dimension(12*vecl)  :: wrk_Jac
      integer,  dimension(12*vecl)  :: jwrk_Jac
      integer,  dimension(vecl)     :: iw
      integer                       :: ierr=0
!------------------------------------------------------------------------------
      
      Program_Step_Select: select case(programstep)
        !**Pre-initialization. Get problem name and vector length**
        case('INITIALIZE_PROBLEM_INFORMATION')
          nvecLen = neq*vecl
          eq = neq
          probname='Brusselat'     
          tol=9.9e-11_wp  
          dt_error_tol=1.0e-13_wp
          Jac_case='SPARSE'
          
          allocate(var_names(neq))
          var_names(:)=(/'Differential', 'Algebraic   '/)
          
          call grid()

        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')
          
          !Allocate derivative operators
          if(.not. allocated(D2_per))call Define_CSR_Operators(vecl,dx) 
           
          !Time information
          choose_dt: select case(temporal_splitting)
            case('IMPLICIT');dt = 0.25_wp/10**((iDT-1)/28.0_wp) ! timestep   
            case('IMEX')    ;dt = 0.0025_wp/10**((iDT-1)/40.0_wp)
          end select choose_dt
          tfinal = 10.0_wp                   ! final time
 
          !**IC**
          uvec(1:2*vecl:2) = 1.0_wp+sin(two*pi*x(:))
          uvec(2:2*vecl:2) = 3.0_wp
      
          !**Exact Solution**
          uexact=exact_Bruss(ep)
          
        case('BUILD_RHS')
          call Bruss_dudt(uvec,dudt_Dx,dudt_Source,ep)
          choose_RHS_type: select case (Temporal_Splitting)
            case('IMPLICIT')
              resE_vec(:)=0.0_wp
              resI_vec(:)=dt*(dudt_Dx(:)+dudt_Source(:))
            case('IMEX')
              resE_vec(:)=dt*dudt_Dx(:)
              resI_vec(:)=dt*dudt_Source(:)
          end select choose_RHS_type
       
        case('BUILD_JACOBIAN')
          choose_Jac_type: select case(temporal_splitting)
            case('IMPLICIT')
              nnz_Jac=2*nnz_D2_per+vecl*2
              call Allocate_Jac_CSR_Storage(vecl*2,nnz_Jac)
              
              call Build_Source_Jac(uvec,Source,jSource,iSource)
              
              call aplb(vecl*2,vecl*2,1,Source,jSource,iSource,               &
     &                  Deriv2_p,jDeriv2_p,iDeriv2_p,                         &
     &                  wrk_Jac,jwrk_Jac,iwrk_Jac,nnz_Jac,iw,ierr)
              if (ierr/=0) then; print*,'Build Jac ierr=',ierr;stop;endif
              
              call Build_Jac(dt,akk,wrk_Jac,jwrk_Jac,iwrk_Jac)   
                  
            case('IMEX')
              nnz_Jac=vecl*4
              call Allocate_Jac_CSR_Storage(vecl*2,nnz_Jac)
              call Build_Source_Jac(uvec,Source,jSource,iSource)
              call Build_Jac(dt,akk,Source,jSource,iSource)
          end select choose_Jac_type
          
      end select Program_Step_select
      end subroutine Brusselator
!==============================================================================
!==============================================================================
!==============================================================================
!******************************************************************************
! Subroutine to initialize x-grid
!******************************************************************************
! MODULE VARIABLES/ROUTINES:
! vecl -> total length of u-vector used for problem, integer,                                           not modified
! x    -> stores x-grid points,                      real(wp),  dimension(u-vector length / num. eq's),     set
! dx   -> separation between grid points,            real(wp),                                              set
! xL   -> Left x-bound                               real(wp),                                          not modified
! xR   -> Right x-bound                              real(wp),                                          not modified
!******************************************************************************
      subroutine grid()
      integer :: i
      
      do i=1,vecl
        x(i)= xL + (xR-xL)*(i-1.0_wp)/(vecl)
      enddo
      
      dx=x(2)-x(1)

      end subroutine grid
!==============================================================================
!******************************************************************************
! Function to return exact solution vector
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
! eps         -> Stiffness epsilon value, real(wp)
! OUTPUTS:
! exact_Bruss -> exact solution vector,   real(wp), dimension(u-vector length)
!******************************************************************************
      function exact_Bruss(eps)

      real(wp),                  intent(in) :: eps

      integer                               :: i
      real(wp), dimension(81,vecl*2+1)      :: ExactTot
      real(wp)                              :: diff
      real(wp), dimension(vecl*2)           :: exact_Bruss
      character(len=2)                      :: vstr

      exact_Bruss(:)=0.0_wp
      
      !**Exact Solution** 
      write(vstr,"(I2.1)")vecl
      open(unit=39,file='exact.Brusselator_'//vstr//'.data')
      rewind(39)
      do i=1,81
        read(39,*)ExactTot(i,1:vecl*2)
        ExactTot(i,vecl*2+1) = 1.0_wp/10**((i-1)/(10.0_wp))  !  used for 81 values of ep
      enddo
      do i=1,81
        diff = abs(ExactTot(i,vecl*2+1) - eps)
        if(diff <= 1.0e-10_wp)then
          exact_Bruss(:) = ExactTot(i,:vecl*2)
          exit
        endif
      enddo

      return
      end function exact_Bruss

!==============================================================================
!******************************************************************************
! Subroutine to set RHS and store spacial part of Jacobian
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! SBP_COEF_MODULE.F90       *DEFINES CSR OPERATORS 
! UNARY_MOD.F90             *PERFORMS SPARSE MATRIX OPERATIONS
! MATVEC_MODULE.F90         *PERFORMS SPARSE MATRIX*VECTOR OPERATIONS
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From SBP_Coef_Module:
!    D2_per    -> Second derivative  a matrix periodic operator,  real(wp), dimension(5*u-vector length/num. eq's), not modified
!   jD2_per    -> Second derivative ja matrix periodic operator,  integer,  dimension(5*u-vector length/num. eq's), not modified
!   iD2_per    -> Second derivative ia matrix periodic operator,  integer,  dimension(u-vector length + 1),         not modified
!   nnz_D2_per -> Second derivative operator number of non zeros, integer,                                          not modified
! From unary_mod:
!   aplb  -> Subroutine to add one CSR matrix to another
!   dperm -> Subroutine to permute a matrix
!   csort -> subroutine to sort a CSR matrix
! From matvec_module:
!   amux -> multiplies a vector into a CSR matrix
!
! MODULE VARIABLES/ROUTINES:
! vecl      -> total length of u-vector used for problem,      integer,                                                   not modified
! jDeriv2_p -> permuted ja matrix for D1 operators             integer,   dimension(2 * 5 * u-vector length / num. eq's), not modified
! Deriv2_p  -> permuted  a matrix for D1 operators             real(wp),  dimension(2 * 5 * u-vector length / num. eq's), not modified
! iDeriv2_p -> permuted ia matrix for D1 operators             integer,   dimension(u-vector length + 1),                 not modified
!******************************************************************************
! INPUTS:
! vec -> variable vector,         real(wp), dimension(u-vector length)
! eps -> Stiffness epsilon value, real(wp)
! OUTPUTS:
! dudt_deriv2 -> RHS vector involving just spacial derivatives, real(wp), dimension(u-vector length)
! dudt_source -> RHS vector involving just source terms         real(wp), dimension(u-vector length)
!******************************************************************************
      subroutine Bruss_dUdt(vec,dudt_deriv2,dudt_source,eps)
      
      use SBP_Coef_Module, only: D2_per,jD2_per,iD2_per,nnz_D2_per
      use unary_mod,       only: dperm
      use matvec_module,   only: amux

      real(wp), dimension(:),      intent(in   ) :: vec
      real(wp), dimension(vecl*2), intent(  out) :: dudt_deriv2,dudt_source
      real(wp),                    intent(in   ) :: eps
     
      real(wp), dimension(vecl)                  :: u,v
      integer,  dimension(vecl*2)                :: j_perm
      integer                                    :: i      
      integer,  dimension(vecl*2+1)              :: iDeriv2      
      real(wp), dimension(2*nnz_D2_per)          :: Deriv2
      integer,  dimension(2*nnz_D2_per)          :: jDeriv2
            
      u=vec(1:vecl*2:2); v=vec(2:vecl*2:2)

! ---------------------- Dx part of dudt --------------------------------------
      Deriv2(:nnz_D2_per)    = eps*D2_per(:)
      Deriv2(nnz_D2_per+1:)  = eps*D2_per(:)
      iDeriv2(:vecl)         = iD2_per(:)
      iDeriv2(vecl+1:)       = iD2_per(:)+iD2_per(vecl+1)-1
      jDeriv2(:nnz_D2_per)   = jD2_per(:)
      jDeriv2(nnz_D2_per+1:) = jD2_per(:)+vecl

! permute matrix into:
! |x| |
! | |x|
      j_perm(1)=1
      j_perm(vecl+1)=2
      do i=2,vecl*2
        if (i==vecl+1) cycle
        j_perm(i)=j_perm(i-1) + 2
      enddo
 
      ! permute D2's matrix   
      call dperm(vecl*2,Deriv2,jDeriv2,iDeriv2,Deriv2_p, &
     &           jDeriv2_p,iDeriv2_p,j_perm,j_perm,1)

      ! get dudt     
      call amux(vecl*2,vec,dudt_Deriv2,Deriv2_p,jDeriv2_p,iDeriv2_p)   

!----------------------  Source part of dudt  ---------------------------------      
      dudt_Source(1:2*vecl:2)=1.0_wp-4.0_wp*u(:)+u(:)*u(:)*v(:)
      dudt_Source(2:2*vecl:2)=       3.0_wp*u(:)-u(:)*u(:)*v(:)

      end subroutine Bruss_dUdt
!==============================================================================
!******************************************************************************
! Subroutine to set source Jacobian
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! UNARY_MOD.F90             *PERFORMS SPARSE MATRIX OPERATIONS
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From unary_mod:
!   dperm -> Subroutine to permute a matrix
!
! MODULE VARIABLES/ROUTINES
! vecl -> total length of u-vector used for problem, integer, not modified
!******************************************************************************
! INPUTS:
! vec       -> variable vector,                                      real(wp), dimension(u-vector length)
! ep        -> stiffness constant (epsilon),                         real(wp)
! OUTPUTS:
! iSource_p -> ia combined source matrix for output to main routine, integer,  dimension(u-vector length + 1)
! jSource_p -> ja combined source matrix for output to main routine, integer,  dimension(u-vector length)
!  Source_p ->  a combined source matrix for output to main routine, real(wp), dimension(u-vector length)
!******************************************************************************
      subroutine Build_Source_Jac(vec,Source_p,jSource_p,iSource_p)
      
      use unary_mod, only: dperm

      real(wp), dimension(:),        intent(in   ) :: vec
      real(wp), dimension(4*vecl),   intent(  out) :: Source_p
      integer,  dimension(4*vecl),   intent(  out) :: jSource_p
      integer,  dimension(2*vecl+1), intent(  out) :: iSource_p
      
      real(wp), dimension(4*vecl)   :: Source      
      integer,  dimension(4*vecl)   :: jSource
      integer,  dimension(2*vecl+1) :: iSource
      
      integer,  dimension(vecl*2)   :: j_perm
      real(wp), dimension(vecl)     :: u,v
      integer                       :: i
!------------------------------------------------------------------------------
      
      u=vec(1:vecl*2:2); v=vec(2:vecl*2:2)
     
! Set Source Jacobian
      Source(1:vecl*2:2)=-4.0_wp+2.0_wp*u(:)*v(:)
      Source(2:vecl*2:2)=2.0_wp*u(:)
      Source(vecl*2+1:vecl*4:2)=3.0_wp-2.0_wp*u(:)*v(:)
      Source(vecl*2+2:vecl*4:2)=-u(:)*u(:)
      iSource(:) = (/ (i, i=1, vecl*4+1,2) /)
      do i=1,vecl
        jSource(2*i-1)        = i
        jSource(2*i)          = vecl+i
        jSource(2*i-1+vecl*2) = i
        jSource(2*i+vecl*2)   = vecl+i
      enddo
    
      j_perm(1)=1
      j_perm(vecl+1)=2
      do i=2,vecl*2
        if (i==vecl+1) cycle
        j_perm(i)=j_perm(i-1) + 2
      enddo

      ! permute source's matrix    
      call dperm(vecl*2,Source,jSource,iSource,Source_p,jSource_p,iSource_p,  &
     &             j_perm,j_perm,1)
     
      end subroutine Build_Source_Jac
!==============================================================================
!******************************************************************************
! Subroutine to finish building Jacobian and set it to global variables
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! UNARY_MOD.F90             *PERFORMS SPARSE MATRIX OPERATIONS
! JACOBIAN_CSR_MOD.F90      *ALLOCATE AND STORE CSR JACOBIAN VARIABLES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From unary_mod:
!   aplsca -> Subroutine to add diagonal constant to a CSR matrix
! From Jacobian_CSR_Mod:
!   iaJac -> ia matrix for global storage of Jacobian, integer,  dimension(u-vector length + 1),                               set
!   jaJac -> ja matrix for global storage of Jacobian, integer,  dimension(dependant on temporal_splitting, see main routine), set
!    aJac ->  a matrix for global storage of Jacobian, real(wp), dimension(dependant on temporal_splitting, see main routine), set
!
! MODULE VARIABLES/ROUTINES:
! vecl  -> total length of u-vector used for problem, integer, not modified
!******************************************************************************
! INPUTS:
! dt  -> timestep,                     real(wp)
! akk -> Diagonal term from RK scheme, real(wp)
! a   ->  a matrix for input Jacobian, real(wp), dimension(dependant on temporal_splitting, see main routine)
! ja  -> ja matrix for input Jacobian, integer,  dimension(dependant on temporal_splitting, see main routine)
! ia  -> ia matrix for input Jacobian, intenger, dimension(u-vector length + 1)
!******************************************************************************    
      subroutine Build_Jac(dt,akk,a,ja,ia)
      
      use unary_mod,        only: aplsca
      use Jacobian_CSR_Mod, only: iaJac,jaJac,aJac
      
      real(wp),               intent(in) :: dt,akk
      real(wp), dimension(:), intent(in) :: a
      integer,  dimension(:), intent(in) :: ja,ia
      
      integer,  dimension(vecl*2)   :: iw    
      real(wp), dimension(size(a))  :: wrk
      integer,  dimension(size(a))  :: jwrk
      integer,  dimension(vecl*2+1) :: iwrk
      
      !Store in temporary variables to allow a,ja,ia to be intent(in)
      wrk(:)=-akk*dt*a(:)
      jwrk(:)=ja(:)
      iwrk(:)=ia(:)  
      
      call aplsca(vecl*2,wrk,jwrk,iwrk,1.0_wp,iw)      

      iaJac=iwrk
      jaJac=jwrk
      aJac = wrk

      end subroutine Build_Jac
!==============================================================================    
      end module Brusselator_mod
