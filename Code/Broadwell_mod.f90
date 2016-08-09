!******************************************************************************
! Module containing routines to describe Broadwell's linear system in "ON A 
! CLASS OF UNIFORMLY ACCURATE IMEX RUNGE–KUTTA SCHEMES AND APPLICATIONS TO 
! HYPERBOLIC SYSTEMS WITH RELAXATION" pg 12, equations 34
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
! SBP_COEF_MODULE.F90           *DEFINES CSR OPERATORS 
! JACOBIAN_CSR_MOD.F90          *ALLOCATE AND STORE CSR JACOBIAN VARIABLES
! UNARY_MOD.F90                 *PERFORMS SPARSE MATRIX OPERATIONS
! MATVEC_MODULE.F90             *PERFORMS SPARSE MATRIX*VECTOR OPERATIONS
!******************************************************************************
      module Broadwell_Mod 

      use precision_vars, only: wp,pi

      implicit none; save 
  
      private
      public  :: Broadwell

!--------------------------------VARIABLES-------------------------------------         
      integer,  parameter :: nx   =  32 
      integer,  parameter :: neq  =   3 

      integer,  parameter :: vecl = nx * neq

      real(wp), parameter :: Length =  20.0_wp
      real(wp), parameter :: xL     = -10.0_wp
      real(wp), parameter :: xR     = +10.0_wp      
      real(wp), parameter :: a_u1   = + 0.3_wp
      real(wp), parameter :: a_u2   = + 0.1_wp
   
      real(wp), dimension(nx)     :: x      
      real(wp)                    :: dx

      integer,  dimension(vecl+1) :: iSource_p,iDeriv_comb_p,iwrk_Jac
      integer,  dimension(4*vecl) :: jDeriv_comb_p
      integer,  dimension(5*vecl) :: jwrk_Jac
      integer,  dimension(vecl)   :: jSource_p

      real(wp), dimension(4*vecl) ::  Deriv_comb_p
      real(wp), dimension(5*vecl) :: wrk_Jac
      
      real(wp), dimension(vecl)   ::  Source_p
      logical :: update_Jac
!------------------------------------------------------------------------------      
      contains    
!==============================================================================
!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Broadwell problem 
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
! nx            -> number of grid points                           integer,                                                             not modified
! neq           -> number of equations in problem,                 integer,                                                             not modified
! vecl          -> total length of u-vector used for problem,      integer,                                                             not modified
! Length        -> Length of x-grid                                real(wp),                                                            not modified
! a_u1          -> problem parameter                               real(wp),                                                            not modified
! a_u2          -> problem parameter                               real(wp),                                                            not modified
! x             -> stores x-grid points,                           real(wp),  dimension(num. points),                                   not modified
! dx            -> separation between grid points,                 real(wp),                                                            not modified
! jDeriv_comb_p -> permuted ja matrix for D1 operators             integer,   dimension(num. eq. * 4 * num. points),                    not modified
! Deriv_comb_p  -> permuted  a matrix for D1 operators             real(wp),  dimension(num. eq. * 4 * num. points),                    not modified
! iDeriv_comb_p -> permuted ia matrix for D1 operators             integer,   dimension(u-vector length) + 1),                          not modified
! jwrk_Jac      -> combined ja matrix for Source + D1              integer,   dimension(num. eq's * 4 * num. points + u-vector length),     set
! wrk_Jac       -> combined  a matrix for Source + D1              real(wp),  dimension(num. eq's * 4 * num. points + u-vector length),     set
! iwrk_Jac      -> combined ia matrix for Source + D1              integer,   dimension(u-vector length + 1),                               set
! jSource_p     -> permuted ja matrix for Source terms             integer,   dimension(u-vector length),                               not modified
! Source_p      -> permuted  a matrix for Source terms             real(wp),  dimension(u-vector length),                               not modified
! iSource_p     -> permuted ia matrix for Source terms             integer,   dimension(u-vector length + 1),                           not modified
! update_Jac    -> logical flag set to decide when to update Jac,  logical,                                                                 set & modified
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
      subroutine Broadwell(nveclen,eq,ep,dt,tfinal,iDT,resE_vec,resI_vec,akk)

      use control_variables, only: temporal_splitting,probname,Jac_case,     &
     &                             tol,dt_error_tol,uvec,uexact,programstep, &
     &                             var_names
      use SBP_Coef_Module,   only: Define_CSR_Operators,D1_per
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
      real(wp), dimension(vecl), intent(  out) :: resE_vec,resI_vec
      real(wp), dimension(vecl)                :: dudt_D1,dudt_Source
      
      !Jacob vars
      real(wp), intent(in   )  :: akk
      integer                  :: nnz_Jac,ierr=0
      integer, dimension(vecl) :: iw

      real(wp), dimension(nx)             :: r0,u0,m0,z0
      real(wp), dimension(nx)             :: zE,z1,H,z1_test
      real(wp), dimension(nx)             :: dzEdr, dzEdm, dr0dx, dm0dx
      real(wp)                            :: f

!------------------------------------------------------------------------------   
      
      Program_Step_Select: select case(programstep)
        !**Pre-initialization. Get problem name, vector length and grid**
        case('INITIALIZE_PROBLEM_INFORMATION')
          nveclen = vecl
          eq = neq
          probname='Broadwell'     
          Jac_case='SPARSE'
          tol=1.0e-14_wp  
          dt_error_tol=5.0e-14_wp
          
          allocate(var_names(neq))
          var_names(:)=(/'Differential', 'Differential', 'Algebraic   '/)
          
          call grid()           

        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')
          
          !Update=true for Jac/RHS updates
          update_Jac=.true. !reset every new epsilon/dt

          !Allocate derivative operators
          if(.not. allocated(D1_per))call Define_CSR_Operators(nx,dx)   
          
          !Time information
          tfinal = 0.2_wp  ! final time   
          choose_dt: select case(temporal_splitting)
            case('EXPLICIT'); dt = 0.000025_wp*0.1_wp/10**((iDT-1)/20.0_wp) ! explicit timestep
            case default    ; dt = 0.2_wp/10**((iDT-1)/20.0_wp)             ! implicit timestep      
          end select choose_dt
          
          f = 2.0_wp*pi/Length
          ! Set IC's
          r0(:) = (1.0_wp + a_u1 * sin(f*x(:)) )
          u0(:) = (0.5_wp + a_u2 * sin(f*x(:)) )
          m0(:) = r0(:)*u0(:)
          zE(:) = 0.5_wp * (r0(:)*r0(:) + m0(:)*m0(:)) / r0(:)

          dzEdr = 0.5_wp * ( 1.0_wp - u0(:)*u0(:) )
          dzEdm = u0(:)

          dr0dx(:) = f * a_u1 * cos(f * x(:)) 
          dm0dx(:) = r0(:) * f * a_u2 * cos(f * x(:)) & 
                   + u0(:) * f * a_u1 * cos(f * x(:))
              H(:) = (1.0_wp - dzEdr(:) + dzEdm(:)*dzEdm(:))*dm0dx(:) + dzEdr(:)*dzEdm(:)*dr0dx(:)
             z1(:) = -H(:) / r0(:) /2.0_wp !Boscarino's paper is missing 1/2 multiplier on z1
             z0(:) = zE(:) + ep*z1(:)
           
             !z1(:)=-1.0_wp/2.0_wp*(dm0dx(:)/2.0_wp*(1+u0(:)**2)+u0(:)*(dzEdr(:)*dr0dx(:)+dzEdm(:)*dm0dx(:)))/r0(:)
          uvec(1:vecl-2:neq) = r0(:)
          uvec(2:vecl-1:neq) = m0(:)
          uvec(3:vecl-0:neq) = z0(:)
!------------------------------------------------
        
          
        !  print*,'test',z1_test(:)-z1(:)/2.0_wp



          !set exact solution at tfinal
          call exact_Broadwell(uexact,ep) 
              
        case('BUILD_RHS')
          call Broadwell_dUdt(uvec,dudt_D1,dudt_Source,ep)        
          choose_RHS_type: select case (Temporal_Splitting)       
            case('EXPLICIT')
              resE_vec(:)=dt*(dudt_D1(:)+dudt_Source(:))
              resI_vec(:)=0.0_wp                   
            case('IMPLICIT')
              resE_vec(:)=0.0_wp
              resI_vec(:)=dt*(dudt_D1(:)+dudt_Source(:))
            case('IMEX')
              resE_vec(:)=dt*dudt_D1(:)
              resI_vec(:)=dt*dudt_Source(:)
          end select choose_RHS_type
          
        case('BUILD_JACOBIAN')
          choose_Jac_type: select case(temporal_splitting)
            case('IMPLICIT')
              if (update_Jac) then
                nnz_Jac=5*vecl + 2*nx
                call Allocate_Jac_CSR_Storage(vecl,nnz_Jac)
             
                call Broadwell_Build_Spatial_Jac(Deriv_comb_p,jDeriv_comb_p,iDeriv_comb_p)
                update_Jac=.false.  !no need to update matrix that forms spacial Jacobian
              endif

              call Broadwell_Build_Source_Jac(ep,uvec,Source_p,jSource_p,iSource_p)

              call aplb(vecl,vecl,1,Source_p,jSource_p,iSource_p,    &
     &                  Deriv_comb_p,jDeriv_comb_p,iDeriv_comb_p,    &
     &                  wrk_Jac,jwrk_Jac,iwrk_Jac,5*vecl,iw,ierr) 
              if (ierr/=0) then; print*,'Build Jac ierr=',ierr; stop; endif 
              
              call Broadwell_Add_Diag_Jac(dt,akk,wrk_Jac,jwrk_Jac,iwrk_Jac)

            case('IMEX')
              nnz_Jac=vecl + 2*nx
              call Allocate_Jac_CSR_Storage(vecl,nnz_Jac)
              call Broadwell_Build_Source_Jac(ep,uvec,Source_p,jSource_p,iSource_p)
              call Broadwell_Add_Diag_Jac(dt,akk,Source_p,jSource_p,iSource_p)
              
          end select choose_Jac_type
      end select Program_Step_Select     

      end subroutine Broadwell
!==============================================================================
!==============================================================================
!==============================================================================
!******************************************************************************
! Subroutine to initialize x-grid
!******************************************************************************
! MODULE VARIABLES/ROUTINES:
! neq  -> number of equations in problem,            integer,                          not modified
! x    -> stores x-grid points,                      real(wp), dimension(num. points),     set
! dx   -> separation between grid points,            real(wp),                             set
! xL   -> Left x-bound                               real(wp),                         not modified
! xR   -> Right x-bound                              real(wp),                         not modified
!******************************************************************************
      subroutine grid()
      integer :: i
      
      do i=1,nx
        x(i)= xL + (xR-xL)*(i-1.0_wp)/(nx)
      enddo
      
      dx=x(2)-x(1)

      end subroutine grid
!==============================================================================
!******************************************************************************
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
      subroutine exact_Broadwell(u,ep)

      real(wp),                  intent(in   ) :: ep
      real(wp), dimension(vecl), intent(  out) :: u

      integer                                  :: i
      real(wp), dimension(81,vecl+1)           :: ExactTot
      real(wp) :: diff

      u(:)=0.0_wp*ep
      
      !**Exact Solution** 
      open(unit=39,file='exact.Broadwell_32.data')
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

      end subroutine exact_Broadwell

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
!   wp  -> working precision
! From SBP_Coef_Module:
!    D1_per -> First derivative  a matrix periodic operator, real(wp), dimension(4*u-vector length/num. eq's), not modified
!   jD1_per -> First derivative ja matrix periodic operator, integer,  dimension(4*u-vector length/num. eq's), not modified
!   iD1_per -> First derivative ia matrix periodic operator, integer,  dimension(u-vector length + 1),         not modified
! From matvec_module:
!   amux -> multiplies a vector into a CSR matrix
!
! MODULE VARIABLES/ROUTINES:
! vecl -> total length of u-vector used for problem, integer, not modified
!******************************************************************************
! INPUTS:
! u   -> variable vector,         real(wp), dimension(u-vector length)
! eps -> Stiffness epsilon value, real(wp)
! OUTPUTS:
! dudt_D1     -> RHS vector involving just spacial derivatives, real(wp), dimension(u-vector length)
! dudt_source -> RHS vector involving just source terms         real(wp), dimension(u-vector length)
!******************************************************************************
      subroutine Broadwell_dUdt(u,dudt_D1,dudt_Source,ep)

!     u ordering:  flatted(  (neq,nx) ) 
!     dudt_D1:     flatted(  (neq,nx) ) 
      
      use SBP_Coef_Module, only: D1_per,jD1_per,iD1_per
      use matvec_module,   only: amux

      real(wp), dimension(:),    intent(in   ) :: u
      real(wp), dimension(vecl), intent(  out) :: dudt_D1,dudt_Source
      real(wp),                  intent(in   ) :: ep
           
      real(wp)                            :: epI
      real(wp), dimension(nx)             :: u1,u2,u3
      real(wp), dimension(nx)             :: t1

      u1(:) = u(1:vecl-2:3)
      u2(:) = u(2:vecl-1:3)
      u3(:) = u(3:vecl-0:3)

      dudt_D1(:)     = 0.0_wp
      dudt_Source(:) = 0.0_wp

      epI  = 1.0_wp / ep
!------------------------------------------------------------------------------     
      ! d (u2) / dx
      call amux(nx,u2,t1,D1_per,jD1_per,iD1_per)
      dudt_D1(1:vecl-2:3) = -t1(:)
      dudt_D1(3:vecl-0:3) = -t1(:)

      ! d (u3) / dx
      call amux(nx,u3,t1,D1_per,jD1_per,iD1_per)
      dudt_D1(2:vecl-1:3) = -t1(:)

      dudt_source(3:vecl-0:3) = epI * ( u1(:)*u1(:) + u2(:)*u2(:) - 2.0_wp*u1(:)*u3(:))

      end subroutine Broadwell_dUdt

!==============================================================================
!******************************************************************************
! Subroutine to set spatial Jacobian
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! SBP_COEF_MODULE.F90       *DEFINES CSR OPERATORS 
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From SBP_Coef_Module:
!    D1_per    -> First derivative  a matrix periodic operator, real(wp), dimension(4*num. points), not modified
!   jD1_per    -> First derivative ja matrix periodic operator, integer,  dimension(4*num. points), not modified
!   nnz_D1_per -> First derivative matrix number of non zeros,  integer,                            not modified
!
! MODULE VARIABLES/ROUTINES
! vecl -> total length of u-vector used for problem, integer, not modified
!******************************************************************************
! OUTPUTS:
! ia -> ia combined derivative matrix for output to main routine, integer,  dimension(u-vector length + 1)
! ja -> ja combined derivative matrix for output to main routine, integer,  dimension(4*num. points)
! a  ->  a combined derivative matrix for output to main routine, real(wp), dimension(4*num. points)
!******************************************************************************
      subroutine Broadwell_Build_Spatial_Jac(a,ja,ia)

      use unary_mod,         only: dperm
      use SBP_Coef_Module,   only: D1_per,jD1_per,nnz_D1_per

      integer,  dimension(vecl+1),     intent(  out) :: ia
      integer,  dimension(nnz_D1_per), intent(  out) :: ja
      real(wp), dimension(nnz_D1_per), intent(  out) ::  a
      
      integer,  dimension(vecl+1)                 :: iDeriv_comb
      integer,  dimension(neq*nnz_D1_per)         :: jDeriv_comb
      real(wp), dimension(neq*nnz_D1_per)         ::  Deriv_comb

      integer,  dimension(vecl)                   :: j_perm

      integer                                     :: i, nnz

!------------------------------------------------------------------------------     
! |0|D|0|
! |0|0|D|
! |0|D|0|

        nnz = nnz_D1_per

        do i = 1,vecl+1
          iDeriv_comb(i) = 1 + 4*(i-1)
        enddo 
        
        jDeriv_comb(      1:1*nnz) =  1*nx + jD1_per(1:nnz)
        jDeriv_comb(1*nnz+1:2*nnz) =  2*nx + jD1_per(1:nnz)
        jDeriv_comb(2*nnz+1:3*nnz) =  1*nx + jD1_per(1:nnz)

         Deriv_comb(      1:1*nnz) =       +  D1_per(1:nnz)
         Deriv_comb(1*nnz+1:2*nnz) =       +  D1_per(1:nnz)
         Deriv_comb(2*nnz+1:3*nnz) =       +  D1_per(1:nnz)

        j_perm(     1)=1
        j_perm(  nx+1)=2
        j_perm(2*nx+1)=3
        do i=2,vecl
          if (i==nx+1 .or. i==2*nx+1) cycle
          j_perm(i)=j_perm(i-1) + 3
        enddo

      ! permute D1's matrix   
        call dperm(vecl,Deriv_comb,jDeriv_comb,iDeriv_comb,a, &
     &             ja,ia,j_perm,j_perm,1)

      end subroutine Broadwell_Build_Spatial_Jac      

!==============================================================================
!******************************************************************************
! Subroutine to set source Jacobian
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
!
! MODULE VARIABLES/ROUTINES
! vecl -> total length of u-vector used for problem, integer, not modified
!******************************************************************************
! INPUTS:
! u  -> variable vector,                                      real(wp), dimension(u-vector length)
! ep -> stiffness constant (epsilon),                         real(wp)
! OUTPUTS:
! ia -> ia combined source matrix for output to main routine, integer,  dimension(u-vector length + 1)
! ja -> ja combined source matrix for output to main routine, integer,  dimension(u-vector length)
! a  ->  a combined source matrix for output to main routine, real(wp), dimension(u-vector length)
!******************************************************************************
      subroutine Broadwell_Build_Source_Jac(ep,u,a,ja,ia)

      real(wp),                    intent(in   ) :: ep
      real(wp), dimension(vecl),   intent(in   ) :: u
      
      integer,  dimension(vecl+1), intent(  out) :: ia
      integer,  dimension(vecl),   intent(  out) :: ja
      real(wp), dimension(vecl),   intent(  out) :: a
      
      integer                                    :: i,ii

      real(wp)                                   :: epI

!------------------------------------------------------------------------------     

      epI = 1.0_wp / ep
      
      ia(1:3)=1
      do i = 4,vecl+1,3 
        ia(i:i+2)=ia(i-3:i-1) + 3
      enddo
      
      do i = 1,nx
         ii = (i-1)*neq

         a(ii+1) = + epI * 2.0_wp*(u(ii+1)-u(ii+3))
         a(ii+2) = + epI * 2.0_wp*(u(ii+2)        )
         a(ii+3) = - epI * 2.0_wp*(u(ii+1)        )
      enddo

      do i = 1,vecl
        ja(i) =  i
      enddo

      end subroutine Broadwell_Build_Source_Jac      

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
! m   ->  a matrix for input Jacobian, real(wp), dimension(dependant on temporal_splitting, see main routine)
! jm  -> ja matrix for input Jacobian, integer,  dimension(dependant on temporal_splitting, see main routine)
! im  -> ia matrix for input Jacobian, intenger, dimension(u-vector length + 1)
!******************************************************************************
      subroutine Broadwell_Add_Diag_Jac(dt,akk,a,ja,ia)
      
      use unary_mod,         only: aplsca
      use Jacobian_CSR_Mod,  only: iaJac, jaJac,  aJac

      real(wp),               intent(in) :: dt,akk
      real(wp), dimension(:), intent(in) :: a
      integer,  dimension(:), intent(in) :: ja,ia
      
      integer,  dimension(vecl)          :: iw
      real(wp), dimension(size(a)+nx*2)  :: wrk
      integer,  dimension(size(a)+nx*2)  :: jwrk
      integer,  dimension(vecl+1)        :: iwrk   
      integer                            :: nnz
      
!------------------------------------------------------------------------------     
      nnz=size(a)
      
      wrk(:nnz)=-akk*dt*a(:)   
      wrk(nnz+1:)=0.0_wp
      jwrk(:nnz)=ja(:)
      jwrk(nnz+1:)=0
      iwrk(:)=ia(:)

      call aplsca(vecl,wrk,jwrk,iwrk,1.0_wp,iw)      

      iaJac=iwrk
      jaJac=jwrk
       aJac= wrk
      
      end subroutine Broadwell_Add_Diag_Jac      

!==============================================================================

      end module Broadwell_Mod
