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
      module Charney_DeVore10_mod

      use precision_vars,    only: wp, pi, eyeN
      use control_variables, only: temporal_splitting,probname,xjac,var_names,&
     &                             tol,dt_error_tol,uvec,uexact,programstep
      use eispack_module,    only: rg, qsortd


      implicit none; save  
  
      private
      public  :: Charney_DeVore10


!--------------------------------CONSTANTS-------------------------------------         
      real(wp), parameter :: sqrt2 = sqrt(2.0_wp)


!--------------------------------VARIABLES-------------------------------------         
      integer,  parameter :: neq  =   10 
      integer,  parameter :: vecl =   10

      real(wp), parameter :: CC  = +0.01_wp   !  Cep is governed from the main program

      !  DeStart, Acta1988
!     real(wp), parameter :: x1s = +2.50000_wp
!     real(wp), parameter :: x4s = -4.00000_wp
!     real(wp), parameter :: be  = +2.00_wp
!     real(wp), parameter :: ga  = +1.00_wp
!     real(wp), parameter :: b   = +1.60_wp

      !  Babaee, Acta1988
      real(wp), parameter :: x1s = +0.95_wp
      real(wp), parameter :: x4s = -0.76095_wp
      real(wp), parameter :: be  = +1.25_wp
      real(wp), parameter :: ga  = +0.2_wp
      real(wp), parameter :: b   = +0.5_wp

!     ============  Original definitions of coefficients  ===================
!     real(wp), parameter :: anm = ( (8*sqrt2*n*b)/(pi) * (m*m)/(4*m*m-1) * &
!                                    (n*n*b*b+(m*m-1)) / (n*n*b*b+m*m)  )

!     real(wp), parameter :: dnm = ( (64.0_wp*sqrt2*n*b)/(15*pi) * &
!                                    (n*n*b*b-(m*m-1))/(n*n*b*b+m*m)  )

!     real(wp), parameter :: epn = ( (16*sqrt2*n*b)/(5*pi) )

!     real(wp), parameter :: rnm = ((9*b/2) * ( ((n-2)*b)**2 - (m-2)**2 ) / (n*n*b*b+m*m) 

!     real(wp), parameter :: bnm = be * ( (n*b)/(n*n*b*b+m*m) )

!     real(wp), parameter :: gnm = ga * ( (4*m*m*m)/(4*m*m-1) *      &
!                                &   (sqrt2*n*b)/(pi*(n*n*b*b+m*m)) )

!     real(wp), parameter :: gnmS= ga * ( (4*m)/(4*m*m-1) * (sqrt2*n*b)/(pi) )

!     real(wp), parameter :: gnmP= ga * ( (3*b)/(4*(n*n*b*b+m*m)) )

!     ============  Original definitions of coefficients  ===================

      real(wp), parameter :: a11 = ( (8*sqrt2*1*b)/(pi) * (1*1)/(4*1*1-1) * &
                                     (1*1*b*b+(1*1-1)) / (1*1*b*b+1*1)  )
      real(wp), parameter :: a12 = ( (8*sqrt2*1*b)/(pi) * (2*2)/(4*2*2-1) * &
                                     (1*1*b*b+(2*2-1)) / (1*1*b*b+2*2)  )
      real(wp), parameter :: a21 = ( (8*sqrt2*2*b)/(pi) * (1*1)/(4*1*1-1) * &
                                     (2*2*b*b+(1*1-1)) / (2*2*b*b+1*1)  )
      real(wp), parameter :: a22 = ( (8*sqrt2*2*b)/(pi) * (2*2)/(4*2*2-1) * &
                                     (2*2*b*b+(2*2-1)) / (2*2*b*b+2*2)  )

      real(wp), parameter :: b11 = be * ( (1*b)/(1*1*b*b+1*1) )
      real(wp), parameter :: b12 = be * ( (1*b)/(1*1*b*b+2*2) )
      real(wp), parameter :: b21 = be * ( (2*b)/(2*2*b*b+1*1) )
      real(wp), parameter :: b22 = be * ( (2*b)/(2*2*b*b+2*2) ) 

      real(wp), parameter :: g11 = ga * ( (4*1*1*1)/(4*1*1-1) * (sqrt2*1*b)/(pi*(1*1*b*b+1*1)) )
      real(wp), parameter :: g12 = ga * ( (4*2*2*2)/(4*2*2-1) * (sqrt2*1*b)/(pi*(1*1*b*b+2*2)) )
      real(wp), parameter :: g21 = ga * ( (4*1*1*1)/(4*1*1-1) * (sqrt2*2*b)/(pi*(2*2*b*b+1*1)) )
      real(wp), parameter :: g22 = ga * ( (4*2*2*2)/(4*2*2-1) * (sqrt2*2*b)/(pi*(2*2*b*b+2*2)) )
                             
      real(wp), parameter :: d11 = ( (64.0_wp*sqrt2*1*b)/(15*pi) * &
                                     (1*1*b*b-(1*1-1))/(1*1*b*b+1*1)  )
      real(wp), parameter :: d12 = ( (64.0_wp*sqrt2*1*b)/(15*pi) * &
                                     (1*1*b*b-(2*2-1))/(1*1*b*b+2*2)  )
      real(wp), parameter :: d21 = ( (64.0_wp*sqrt2*2*b)/(15*pi) * &
                                     (2*2*b*b-(1*1-1))/(2*2*b*b+1*1)  )
      real(wp), parameter :: d22 = ( (64.0_wp*sqrt2*2*b)/(15*pi) * &
                                     (2*2*b*b-(2*2-1))/(2*2*b*b+2*2)  )

      real(wp), parameter :: r11 = (9*b/2) * ( ((1-2)*b)**2 - (1-2)**2 ) / (1*1*b*b+1*1) 
      real(wp), parameter :: r12 = (9*b/2) * ( ((1-2)*b)**2 - (2-2)**2 ) / (1*1*b*b+2*2) 
      real(wp), parameter :: r21 = (9*b/2) * ( ((2-2)*b)**2 - (1-2)**2 ) / (2*2*b*b+1*1) 
      real(wp), parameter :: r22 = (9*b/2) * ( ((2-2)*b)**2 - (2-2)**2 ) / (2*2*b*b+2*2) 

      real(wp), parameter :: g11S= ga * ( (4*1)/(4*1*1-1) * (sqrt2*1*b)/(pi) )
      real(wp), parameter :: g12S= ga * ( (4*2)/(4*2*2-1) * (sqrt2*1*b)/(pi) )

      real(wp), parameter :: g12P= ga * ( (3*b)/(4*(1*1*b*b+2*2)) )
      real(wp), parameter :: g21P= ga * ( (3*b)/(4*(2*2*b*b+1*1)) )

      real(wp), parameter :: ep1 = ( (16*sqrt2*1*b)/(5*pi) )
      real(wp), parameter :: ep2 = ( (16*sqrt2*2*b)/(5*pi) )

      real(wp), dimension(neq,neq) :: eye= reshape((/                  &
                                   &  +1,+0,+0,+0,+0,+0,+0,+0,+0,+0,   &
                                   &  +0,+1,+0,+0,+0,+0,+0,+0,+0,+0,   &
                                   &  +0,+0,+1,+0,+0,+0,+0,+0,+0,+0,   &
                                   &  +0,+0,+0,+1,+0,+0,+0,+0,+0,+0,   &
                                   &  +0,+0,+0,+0,+1,+0,+0,+0,+0,+0,   &
                                   &  +0,+0,+0,+0,+0,+1,+0,+0,+0,+0,   &
                                   &  +0,+0,+0,+0,+0,+0,+1,+0,+0,+0,   &
                                   &  +0,+0,+0,+0,+0,+0,+0,+1,+0,+0,   &
                                   &  +0,+0,+0,+0,+0,+0,+0,+0,+1,+0,   &
                                   &  +0,+0,+0,+0,+0,+0,+0,+0,+0,+1/), &
                                   &         (/10,10/) ) 

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
! Charney_DeVore10_dudt              -> Subroutine to build dudt (LHS)
! Charney_DeVore10_Build_Jac -> Subroutine to build spatial Jacobian
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

      subroutine Charney_DeVore10(nveclen,eq,Cep,dt,tfinal,iDT,resE_vec,resI_vec,OTD_RHS,akk)

      use control_variables, only: temporal_splitting,probname,Jac_case,     &
     &                             tol,dt_error_tol,uvec,uexact,programstep, &
     &                             var_names, OTD, OTDN, Ovec, wr, wi
      use OTD_Module,        only: Orthogonalize_SubSpace
      use Fourier_Coef_Module, only: Fourier_Coefficients

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
      real(wp), dimension(vecl),      intent(  out) :: resE_vec,resI_vec
      real(wp), dimension(vecl,OTDN), intent(  out) :: OTD_RHS
      real(wp), dimension(vecl)                :: dudt_E,dudt_I
      real(wp), dimension(vecl,vecl)           :: Jac_E, Jac_I, Z
      real(wp), dimension(vecl)                :: r
      
      integer                  :: i
      integer,  dimension(OTDN)                :: ind

      !Jacob vars
      real(wp), intent(in   )  :: akk

      real(wp), dimension(:,:), allocatable  :: Tmat

!------------------------------------------------------------------------------   
      
      Program_Step_Select: select case(programstep)

        !**Pre-initialization. Get problem name, vector length and grid**

        case('INITIALIZE_PROBLEM_INFORMATION')
          nveclen = vecl
          eq      = vecl
          probname='Charney_DeVore10'     
          Jac_case='DENSE'
          tol=1.0e-14_wp  
          dt_error_tol=5.0e-14_wp
          
          allocate(var_names(neq))
          var_names(:)=(/ 'Differential', 'Differential',  &
                          'Differential', 'Differential',  &
                          'Differential', 'Differential',  &
                          'Differential', 'Differential',  &
                          'Differential', 'Differential'/)
          
          call Fourier_Coefficients(9,Tmat)
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
          
          call random_number(r(:)) ;

          uvec(1)  = +0.14_wp*r(1) ;
          uvec(2)  = +2.00_wp*r(2) ;
          uvec(3)  = +0.00_wp*r(3) ;
          uvec(4)  = -1.91_wp*r(4) ;
          uvec(5)  = +1.00_wp*r(5) ;
          uvec(6)  = +0.20_wp*r(6) ;
          uvec(7)  = -0.10_wp*r(7) ;
          uvec(8)  = -0.20_wp*r(8) ;
          uvec(9)  = +0.30_wp*r(9) ;
          uvec(10) = -0.40_wp*r(10) ;

          call exact_Charney_DeVore10(uexact,Cep)                    !set exact solution at tfinal

          if(OTD .eqv. .true.) then

            call Charney_DeVore10_Jacobian_Eigenvalues(uvec,Cep,Z)

!           call random_number(Ovec(:,:))
            Ovec(:,:) = Z(:,1:OTDN)

            call Orthogonalize_SubSpace(vecl,OTDN,Ovec(:,:))

            call Charney_DeVore10_Build_OTD_Ritz(uvec,Ovec,Cep,OTDN,wr,wi)

          endif

          if(maxval(abs(eyeN(OTDN) - matmul(transpose(Ovec),Ovec))) >= tol) then
            write(*,*)'OTDN eqns are not orthogonal;  stopping'
            stop
          endif

        case('BUILD_RHS')

          call Charney_DeVore10_dUdt(uvec,Cep,dudt_E,dudt_I)        

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
          
          if(OTD .eqv. .true.) then
             call Charney_DeVore10_Build_OTD_RHS(uvec,Ovec,Cep,OTDN,OTD_RHS)
             OTD_RHS = dt * OTD_RHS
          endif

        case('BUILD_JACOBIAN')

          call Charney_DeVore10_Build_Jac(uvec,Cep,Jac_E,Jac_I)

          choose_Jac_type: select case(temporal_splitting)

            case('EXPLICIT') ! For fully implicit schemes

              xjac(:,:) = eye(:,:)

            case('IMEX') ! For IMEX schemes

              xjac(:,:) = eye(:,:) - akk*dt* (             + Jac_I(:,:))

            case('IMPLICIT') ! For fully implicit schemes

              xjac(:,:) = eye(:,:) - akk*dt* (+ Jac_E(:,:) + Jac_I(:,:))

          end select choose_Jac_type

        case('CALCULATE_RITZ_VALUES')

!         call Charney_DeVore10_Jacobian_Eigenvalues(uvec,Cep,Z)

!         call Charney_DeVore10_Build_OTD_Ritz(uvec,Ovec,Cep,OTDN,wr,wi)

      end select Program_Step_Select     

      end subroutine Charney_DeVore10

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

      subroutine exact_Charney_DeVore10(u,ep)

      real(wp),                  intent(in   ) :: ep
      real(wp), dimension(vecl), intent(  out) :: u

      integer                                  :: i
      real(wp), dimension(81,vecl+1)           :: ExactTot
      real(wp) :: diff

      u(:)=0.0_wp*ep
      
      !**Exact Solution** 
      open(unit=39,file='exact.Charney_DeVore10.data')
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

      end subroutine exact_Charney_DeVore10

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

      subroutine Charney_DeVore10_dUdt(u,Cep,dudt_E,dudt_I)

      real(wp), dimension(:), intent(in   ) :: u
      real(wp),               intent(in   ) :: Cep
      real(wp), dimension(:), intent(  out) :: dudt_E, dudt_I
           
      real(wp)                              :: x1,x2,x3,x4,x5, &
                                             & x6,x7,x8,x9,x10

        x1 = u(1); x2 = u(2); x3 = u(3); x4 = u(4); x5  = u(05); 
        x6 = u(6); x7 = u(7); x8 = u(8); x9 = u(9); x10 = u(10);
  
        dudt_E(01) =                   + g11S*x3              
        dudt_E(02) = - (a11*x1-b11)*x3                        - d11*x4*x6  - r11*(x5*x8 - x6*x7)
        dudt_E(03) = + (a11*x1-b11)*x2 - g11 *x1              + d11*x4*x5  + r11*(x5*x7 + x6*x8)

        dudt_E(04) =                   + g12S*x6              + ep1*(x2*x6 - x3*x5) + ep2*(x7*x10 - x8*x9)

        dudt_E(05) = - (a12*x1-b12)*x6                        - d12*x4*x3  + r12*(x2*x8 - x3*x7) + g12p*x8
        dudt_E(06) = + (a12*x1-b12)*x5 - g12 *x4              + d12*x4*x2  - r12*(x2*x7 + x3*x8) - g12p*x7

        dudt_E(07) = - (a21*x1-b21)*x8                        - d21*x4*x10 - r21*(x2*x6 + x3*x5) + g21p*x6
        dudt_E(08) = + (a21*x1-b21)*x7                        + d21*x4*x9  + r21*(x2*x5 - x3*x6) - g21p*x5

        dudt_E(09) = - (a22*x1-b22)*x10                       - d22*x4*x8
        dudt_E(10) = + (a22*x1-b22)*x9                        + d22*x4*x7


        dudt_I(01) =                           - Cep*(x1-x1S)
        dudt_I(02) =                           - Cep* x2    
        dudt_I(03) =                           - Cep* x3    
        dudt_I(04) =                           - Cep*(x4-x4S)
        dudt_I(05) =                           - Cep* x5    
        dudt_I(06) =                           - Cep* x6   
        dudt_I(07) =                           - Cep* x7   
        dudt_I(08) =                           - Cep* x8   
        dudt_I(09) =                           - Cep* x9   
        dudt_I(10) =                           - Cep* x10   


      end subroutine Charney_DeVore10_dUdt

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

      subroutine Charney_DeVore10_Build_Jac(u,Cep,Jac_E,Jac_I)

      real(wp), dimension(:  ), intent(in   ) :: u
      real(wp),                 intent(in   ) :: Cep
      real(wp), dimension(:,:), intent(  out) :: Jac_E, Jac_I
      
      real(wp)                              :: x1,x2,x3,x4,x5, &
                                             & x6,x7,x8,x9,x10

        x1 = u(1); x2 = u(2); x3 = u(3); x4 = u(4); x5  = u(05); 
        x6 = u(6); x7 = u(7); x8 = u(8); x9 = u(9); x10 = u(10);
  
        Jac_E(:,:) = 0.0_wp

        Jac_E(1,3) = g11S
        Jac_E(2,1) = -(a11*x3)
        Jac_E(2,3) = b11 - a11*x1
        Jac_E(2,4) = -(d11*x6)
        Jac_E(2,5) = -(r11*x8)
        Jac_E(2,6) = -(d11*x4) + r11*x7
        Jac_E(2,7) = r11*x6
        Jac_E(2,8) = -(r11*x5)
        Jac_E(3,1) = -g11 + a11*x2
        Jac_E(3,2) = -b11 + a11*x1
        Jac_E(3,4) = d11*x5
        Jac_E(3,5) = d11*x4 + r11*x7
        Jac_E(3,6) = r11*x8
        Jac_E(3,7) = r11*x5
        Jac_E(3,8) = r11*x6
        Jac_E(4,2) = ep1*x6
        Jac_E(4,3) = -(ep1*x5)
        Jac_E(4,5) = -(ep1*x3)
        Jac_E(4,6) = g12S + ep1*x2
        Jac_E(4,7) = ep2*x10
        Jac_E(4,8) = -(ep2*x9)
        Jac_E(4,9) = -(ep2*x8)
        Jac_E(4,10) = ep2*x7
        Jac_E(5,1) = -(a12*x6)
        Jac_E(5,2) = r12*x8
        Jac_E(5,3) = -(d12*x4) - r12*x7
        Jac_E(5,4) = -(d12*x3)
        Jac_E(5,6) = b12 - a12*x1
        Jac_E(5,7) = -(r12*x3)
        Jac_E(5,8) = g12p + r12*x2
        Jac_E(6,1) = a12*x5
        Jac_E(6,2) = d12*x4 - r12*x7
        Jac_E(6,3) = -(r12*x8)
        Jac_E(6,4) = -g12 + d12*x2
        Jac_E(6,5) = -b12 + a12*x1
        Jac_E(6,7) = -g12p - r12*x2
        Jac_E(6,8) = -(r12*x3)
        Jac_E(7,1) = -(a21*x8)
        Jac_E(7,2) = -(r21*x6)
        Jac_E(7,3) = -(r21*x5)
        Jac_E(7,4) = -(d21*x10)
        Jac_E(7,5) = -(r21*x3)
        Jac_E(7,6) = g21p - r21*x2
        Jac_E(7,8) = b21 - a21*x1
        Jac_E(7,10) = -(d21*x4)
        Jac_E(8,1) = a21*x7
        Jac_E(8,2) = r21*x5
        Jac_E(8,3) = -(r21*x6)
        Jac_E(8,4) = d21*x9
        Jac_E(8,5) = -g21p + r21*x2
        Jac_E(8,6) = -(r21*x3)
        Jac_E(8,7) = -b21 + a21*x1
        Jac_E(8,9) = d21*x4
        Jac_E(9,1) = -(a22*x10)
        Jac_E(9,4) = -(d22*x8)
        Jac_E(9,8) = -(d22*x4)
        Jac_E(9,10) = b22 - a22*x1
        Jac_E(10,1) = a22*x9
        Jac_E(10,4) = d22*x7
        Jac_E(10,7) = d22*x4
        Jac_E(10,9) = -b22 + a22*x1

        Jac_I(:,:) = 0.0_wp

        Jac_I(1,1) = -Cep
        Jac_I(2,2) = -Cep
        Jac_I(3,3) = -Cep
        Jac_I(4,4) = -Cep
        Jac_I(5,5) = -Cep
        Jac_I(6,6) = -Cep
        Jac_I(7,7) = -Cep
        Jac_I(8,8) = -Cep
        Jac_I(9,9) = -Cep
        Jac_I(10,10) = -Cep

      end subroutine Charney_DeVore10_Build_Jac      

!==============================================================================
! Subroutine to set spatial OTD RHS 
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

      subroutine Charney_DeVore10_Build_OTD_RHS(uvec,Ovec,Cep,OTDN,OTD_RHS)

      integer,                  intent(in   ) :: OTDN
      real(wp), dimension(:  ), intent(in   ) :: uvec
      real(wp), dimension(:,:), intent(in   ) :: Ovec
      real(wp),                 intent(in   ) :: Cep
      real(wp), dimension(:,:), intent(  out) :: OTD_RHS

      real(wp), dimension(neq,neq)            :: Jac_E, Jac_I, Jac
      real(wp), dimension(neq ,OTDN)          :: Jac_Ovec
      real(wp), dimension(OTDN,OTDN)          :: Ovec_Jac_Ovec
      
        call Charney_DeVore10_Build_Jac(uvec,Cep,Jac_E,Jac_I)
  
        Jac(:,:) = Jac_E(:,:) + Jac_I(:,:)

        Jac_Ovec(:,:) = matmul(Jac(:,:),Ovec(:,:))
        Ovec_Jac_Ovec(:,:) = matmul(transpose(Ovec),Jac_Ovec(:,:))

        OTD_RHS(:,:) = Jac_Ovec(:,:) - matmul(Ovec(:,:),Ovec_Jac_Ovec(:,:))

      end subroutine Charney_DeVore10_Build_OTD_RHS      

!==============================================================================

      subroutine Charney_DeVore10_Build_OTD_Ritz(uvec,Ovec,Cep,OTDN,wr,wi)

      integer,                  intent(in   ) :: OTDN
      real(wp), dimension(:  ), intent(in   ) :: uvec
      real(wp), dimension(:,:), intent(in   ) :: Ovec
      real(wp),                 intent(in   ) :: Cep
      real(wp), dimension(:  ), intent(  out) :: wr,wi

!     local storage for eigenvalue computation

      integer                        :: i

      real(wp), dimension(neq,neq)            :: Jac_E, Jac_I, Jac
      real(wp), dimension(OTDN,OTDN)          :: Ovec_Jac_Ovec, a
      
      integer                         :: nm, n, matz, ierr
      integer,  dimension(OTDN)       :: iv1, ind
      real(wp), dimension(OTDN)       :: fv1
      real(wp), dimension(OTDN,OTDN)  :: z

        call Charney_DeVore10_Build_Jac(uvec,Cep,Jac_E,Jac_I)
  
        Jac(:,:) = Jac_E(:,:) + Jac_I(:,:)

        Ovec_Jac_Ovec(:,:) = matmul(transpose(Ovec),matmul(Jac(:,:),Ovec(:,:)))

        nm = OTDN ; n = OTDN ; a(:,:) = Ovec_Jac_Ovec(:,:) ; matz = 0

        call rg(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr)

        call qsortd(wr,ind,nm)

        write(85,'(20(f8.3,1x))')(wr(ind(nm+1-i)),wi(ind(nm+1-i)),i=1,nm)

      end subroutine Charney_DeVore10_Build_OTD_Ritz      

!==============================================================================

      subroutine Charney_DeVore10_Jacobian_Eigenvalues(uvec,Cep,Z)

      real(wp), dimension(:  ), intent(in   ) :: uvec
      real(wp),                 intent(in   ) :: Cep
      real(wp), dimension(:,:), intent(  out) :: Z

!     local storage for eigenvalue computation


      integer                        :: i,j

      integer                        :: nm, n, matz, ierr
      integer,  dimension(neq)       :: iv1, ind

      real(wp), dimension(neq)       :: fv1, wr, wi
      real(wp), dimension(neq,neq)   :: Jac_E, Jac_I, Jac, a

        call Charney_DeVore10_Build_Jac(uvec,Cep,Jac_E,Jac_I)
  
        Jac(:,:) = Jac_E(:,:) + Jac_I(:,:)

        nm = neq ; n = neq ; a(:,:) = Jac(:,:) ; matz = neq

        call rg(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr)

        call qsortd(wr,ind,neq)

        write(90,'(20(f8.3,1x))')(wr(ind(neq+1-i)),wi(ind(neq+1-i)),i=1,neq)

        Jac_E(:,:) = Z(:,:)
        do i = 1,neq
          do j = 1,neq
            Z(i,j) = Jac_E(i,ind(neq+1-j))
          enddo
        enddo

      end subroutine Charney_DeVore10_Jacobian_Eigenvalues

!==============================================================================

      end module Charney_DeVore10_Mod
