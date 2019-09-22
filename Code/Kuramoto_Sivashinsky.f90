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
      module Kuramoto_Sivashinsky_Mod

      use precision_vars,    only: wp, pi, eyeN, half, third, twothird
      use control_variables, only: temporal_splitting,probname,xjac,var_names,&
     &                             tol,dt_error_tol,uvec,uexact,programstep
      use eispack_module,    only: rg, qsortd


      implicit none; save  
  
      private
      public  :: Kuramoto_Sivashinsky


!--------------------------------CONSTANTS-------------------------------------         
      real(wp), parameter :: sqrt2 = sqrt(2.0_wp)


!--------------------------------VARIABLES-------------------------------------         
      integer,  parameter :: neq  =   33 
      integer,  parameter :: vecl =   33 

!     real(wp), parameter :: xLen = 2 * pi
!     real(wp), parameter :: xLen = 100.0_wp
      real(wp), parameter :: xLen =  50.0_wp

      real(wp), dimension(vecl) :: x
      real(wp)                  :: dx

      logical :: update_Jac

      real(wp), dimension(vecl,vecl) :: Dmat1, Dmat2, Dmat4

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
! Kuramoto_Sivashinsky_dudt              -> Subroutine to build dudt (LHS)
! Kuramoto_Sivashinsky_Build_Jac -> Subroutine to build spatial Jacobian
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

      subroutine Kuramoto_Sivashinsky(nveclen,eq,Cep,dt,tfinal,iDT,resE_vec,resI_vec,akk)

      use control_variables,   only: temporal_splitting, temporal_splitting_OTD,&
                                     probname,Jac_case,                         &
     &                               tol,dt_error_tol,uvec,uexact,programstep,  &
     &                               var_names, OTD, OTDN, Ovec, wr, wi,        &
                                     dudtE_OTD, dudtI_OTD, resE_OTD, resI_OTD,  &
                                     Jac_OTD
      use OTD_Module,          only: Orthogonalize_SubSpace
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
      real(wp)                :: qq
      
      !RHS vars
      real(wp),                       intent(in   ) :: akk
      real(wp), dimension(vecl),      intent(  out) :: resE_vec,resI_vec
      real(wp), dimension(vecl)                     :: dudt_E,dudt_I
      real(wp), dimension(vecl,vecl)                :: Jac_E, Jac_I, Z
      real(wp), dimension(vecl)                     :: r
      
      !  local
      integer                                       :: i
      integer,  dimension(OTDN)                     :: ind


!------------------------------------------------------------------------------   
      
      Program_Step_Select: select case(programstep)

        !**Pre-initialization. Get problem name, vector length and grid**

        case('INITIALIZE_PROBLEM_INFORMATION')
          nveclen     = vecl
          eq          = vecl
          probname    = 'Kuramoto_S'     
          Jac_case    = 'DENSE'
          tol         = 1.0e-9_wp  
          dt_error_tol= 5.0e-14_wp
          
          allocate(var_names(neq))
          var_names(:)= 'Differential'
          
!         Fourier coefficients defined on interval 0 <= xi <= 2 pi  
          call Fourier_Coefficients(vecL,Dmat1,Dmat2,Dmat4)

!         scale coefficient to appropriate length
          qq = 2 * pi / xLen

          Dmat1 = Dmat1 * qq
          Dmat2 = Dmat2 * qq*qq
          Dmat4 = Dmat4 * qq*qq*qq*qq

        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')
          
          !Update=true for Jac/RHS updates
          update_Jac=.true. !reset every new epsilon/dt

          dt_max =   0001.000_wp                                  ! maximum dt
          tfinal =   0250.00_wp                                   ! final time   
          choose_dt: select case(temporal_splitting)
            case('EXPLICIT'); dt = dt_max/10**((iDT-1)/20.0_wp)   ! explicit timestep
            case default    ; dt = dt_max/10**((iDT-1)/20.0_wp)   ! implicit timestep      
          end select choose_dt
          
          call random_number(r(:)) ;

          uvec(:) = r(:)

          call exact_Kuramoto_Sivashinsky(uexact,Cep)                    !set exact solution at tfinal

        case('BUILD_RHS')

          call Kuramoto_Sivashinsky_dUdt(uvec,Cep,dudt_E,dudt_I)        

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

          call Kuramoto_Sivashinsky_Build_Jac(uvec,Cep,Jac_E,Jac_I)

          choose_Jac_type: select case(temporal_splitting)

            case('EXPLICIT') ! For fully implicit schemes

              xjac(:,:) = eyeN(vecl)

            case('IMEX') ! For IMEX schemes

              xjac(:,:) = eyeN(vecl) - akk*dt* (             + Jac_I(:,:))

            case('IMPLICIT') ! For fully implicit schemes

              xjac(:,:) = eyeN(vecl) - akk*dt* (+ Jac_E(:,:) + Jac_I(:,:))

            case('FIRK') ! For fully implicit schemes

              xjac(:,:) = (+ Jac_E(:,:) + Jac_I(:,:))

          end select choose_Jac_type

        case('SET_INITIAL_CONDITIONS_OTD')

          call Kuramoto_Sivashinsky_Jacobian_Eigenvalues(uvec,Cep,Z)

!         call random_number(Ovec(:,:))
          Ovec(:,:) = Z(:,1:OTDN)

          call Orthogonalize_SubSpace(vecl,OTDN,Ovec(:,:))

          call Kuramoto_Sivashinsky_Build_OTD_Ritz(uvec,Ovec,Cep,OTDN,wr,wi)

          if(maxval(abs(eyeN(OTDN) - matmul(transpose(Ovec),Ovec))) >= tol) then
            write(*,*)'OTDN eqns are not orthogonal;  stopping'
            stop
          endif

        case('BUILD_RHS_OTD')

          call Kuramoto_Sivashinsky_Build_OTD_RHS(uvec,Ovec,Cep,OTDN,dudtE_OTD,dudtI_OTD)

          choose_RHS_OTD_type: select case (Temporal_Splitting_OTD)
            case('EXPLICIT') ! For fully explicit schemes
              resE_OTD(:,:)= dt * (+ dudtE_OTD(:,:) + dudtI_OTD(:,:))
              resI_OTD(:,:)= dt * 0.0_wp
            case('IMEX') ! For IMEX schemes
              resE_OTD(:,:)= dt * (+ dudtE_OTD(:,:)                  ) 
              resI_OTD(:,:)= dt * (                  + dudtI_OTD(:,:))
            case('IMPLICIT') ! For fully implicit schemes
          end select choose_RHS_OTD_type

        case('BUILD_JACOBIAN_OTD')

          call Kuramoto_Sivashinsky_Build_Jac(uvec,Cep,Jac_E,Jac_I)
  
          Jac_OTD(:,:) = Jac_E(:,:) + Jac_I(:,:)

          choose_Jac_type_OTD: select case(temporal_splitting_OTD)

            case('EXPLICIT') ! For fully implicit schemes

              xjac(:,:) = eyeN(neq)

            case('IMEX') ! For IMEX schemes

              xjac(:,:) = eyeN(neq) - akk*dt* (             + Jac_OTD(:,:))

          end select choose_Jac_type_OTD

        case('CALCULATE_RITZ_VALUES')

!         call Kuramoto_Sivashinsky_Jacobian_Eigenvalues(uvec,Cep,Z)

!         call Kuramoto_Sivashinsky_Build_OTD_Ritz(uvec,Ovec,Cep,OTDN,wr,wi)

      end select Program_Step_Select     

      end subroutine Kuramoto_Sivashinsky

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

      subroutine exact_Kuramoto_Sivashinsky(u,ep)

      real(wp),                  intent(in   ) :: ep
      real(wp), dimension(vecl), intent(  out) :: u

      integer                                  :: i
      real(wp), dimension(81,vecl+1)           :: ExactTot
      real(wp) :: diff

      u(:)=0.0_wp*ep
      
      !**Exact Solution** 
      open(unit=39,file='./Exact_Data/exact.Kuramoto_Sivashinsky.data')

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

      end subroutine exact_Kuramoto_Sivashinsky

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
!==============================================================================

      subroutine grid()
      integer :: i
      
      do i=1,vecl
        x(i)= 2*pi*(xlen/2 + (xlen)*(i-1.0_wp)/(vecl-1.0_wp))
      enddo
      
      dx=x(2)-x(1)

      end subroutine grid

!==============================================================================

      subroutine Kuramoto_Sivashinsky_dUdt(u,Cep,dudt_E,dudt_I)

      real(wp), dimension(:), intent(in   ) :: u
      real(wp),               intent(in   ) :: Cep
      real(wp), dimension(:), intent(  out) :: dudt_E, dudt_I
           
      real(wp), dimension(vecl)             :: f
      
!--------------------------FLUX------------------------------------------------
      !     Canonical splitting for Inviscid flux in Burgers eqn
      !     f(u)_x = 2/3 (u u/2)_x + 1/3 u u_x

        f(:) = half*u(:)*u(:)

        dudt_E(:) = - twothird*(matmul(Dmat1,f)) - third*u(:)*matmul(Dmat1,u)

        dudt_I(:) = - matmul(Dmat2,u) - matmul(Dmat4,u) 

      end subroutine Kuramoto_Sivashinsky_dUdt

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

      subroutine Kuramoto_Sivashinsky_Build_Jac(u,Cep,Jac_E,Jac_I)

      real(wp), dimension(:  ), intent(in   ) :: u
      real(wp),                 intent(in   ) :: Cep
      real(wp), dimension(:,:), intent(  out) :: Jac_E, Jac_I

      real(wp), dimension(vecl,vecl)          :: DiagUv, DiagDU
      real(wp), dimension(vecl)               :: t1

      integer                                 :: i
      
!       jac = I  -  (- Dmat2 -Dmat4 -2/3*Dmat1*[u]-1/3*[u]*d1-1/3*[d1*u]) * akk * dt
 
        t1 = matmul(Dmat1,u) ;

        diagUv(:,:) = 0.0_wp ; diagDU(:,:) = 0.0_wp ;

        do i = 1,vecl
          diagUv(i,i) =  u(i)
          diagDU(i,i) = t1(i)
        enddo

        Jac_E(:,:) = - twothird* matmul(Dmat1,DiagUv)  &
                   & -    third*(matmul(diagUv,Dmat1) + diagDU(:,:))

        Jac_I(:,:) = - Dmat2(:,:) - Dmat4(:,:)


      end subroutine Kuramoto_Sivashinsky_Build_Jac      

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

      subroutine Kuramoto_Sivashinsky_Build_OTD_RHS(uvec,Ovec,Cep,OTDN,dudtE_OTD,dudtI_OTD)

      integer,                  intent(in   ) :: OTDN
      real(wp), dimension(:  ), intent(in   ) :: uvec
      real(wp), dimension(:,:), intent(in   ) :: Ovec
      real(wp),                 intent(in   ) :: Cep
      real(wp), dimension(:,:), intent(  out) :: dudtE_OTD, dudtI_OTD

      real(wp), dimension(neq,neq)            :: Jac_E, Jac_I, Jac
      real(wp), dimension(neq ,OTDN)          :: Jac_Ovec
      real(wp), dimension(OTDN,OTDN)          :: Ovec_Jac_Ovec
      
        call Kuramoto_Sivashinsky_Build_Jac(uvec,Cep,Jac_E,Jac_I)
  
        Jac(:,:) = Jac_E(:,:) + Jac_I(:,:)

        Jac_Ovec(:,:) = matmul(Jac(:,:),Ovec(:,:))
        Ovec_Jac_Ovec(:,:) = matmul(transpose(Ovec),Jac_Ovec(:,:))

        dudtI_OTD(:,:) = Jac_Ovec(:,:)
        dudtE_OTD(:,:) =               - matmul(Ovec(:,:),Ovec_Jac_Ovec(:,:))

      end subroutine Kuramoto_Sivashinsky_Build_OTD_RHS      

!==============================================================================

      subroutine Kuramoto_Sivashinsky_Build_OTD_Ritz(uvec,Ovec,Cep,OTDN,wr,wi)

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

        call Kuramoto_Sivashinsky_Build_Jac(uvec,Cep,Jac_E,Jac_I)
  
        Jac(:,:) = Jac_E(:,:) + Jac_I(:,:)

        Ovec_Jac_Ovec(:,:) = matmul(transpose(Ovec),matmul(Jac(:,:),Ovec(:,:)))

        nm = OTDN ; n = OTDN ; a(:,:) = Ovec_Jac_Ovec(:,:) ; matz = 0

        call rg(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr)

        call qsortd(wr,ind,nm)

        write(585,'(20(f8.3,1x))')(wr(ind(nm+1-i)),wi(ind(nm+1-i)),i=1,nm)

      end subroutine Kuramoto_Sivashinsky_Build_OTD_Ritz      

!==============================================================================

      subroutine Kuramoto_Sivashinsky_Jacobian_Eigenvalues(uvec,Cep,Z)

      real(wp), dimension(:  ), intent(in   ) :: uvec
      real(wp),                 intent(in   ) :: Cep
      real(wp), dimension(:,:), intent(  out) :: Z

!     local storage for eigenvalue computation


      integer                        :: i,j

      integer                        :: nm, n, matz, ierr
      integer,  dimension(neq)       :: iv1, ind

      real(wp), dimension(neq)       :: fv1, wr, wi
      real(wp), dimension(neq,neq)   :: Jac_E, Jac_I, Jac, a

        call Kuramoto_Sivashinsky_Build_Jac(uvec,Cep,Jac_E,Jac_I)
  
        Jac(:,:) = Jac_E(:,:) + Jac_I(:,:)

        nm = neq ; n = neq ; a(:,:) = Jac(:,:) ; matz = neq

        call rg(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr)

        call qsortd(wr,ind,neq)

        write(590,'(20(f8.3,1x))')(wr(ind(neq+1-i)),wi(ind(neq+1-i)),i=1,neq)

        Jac_E(:,:) = Z(:,:)
        do i = 1,neq
          do j = 1,neq
            Z(i,j) = Jac_E(i,ind(neq+1-j))
          enddo
        enddo

      end subroutine Kuramoto_Sivashinsky_Jacobian_Eigenvalues

!==============================================================================

      end module Kuramoto_Sivashinsky_Mod
