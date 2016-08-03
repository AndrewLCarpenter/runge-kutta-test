!******************************************************************************
! Module containing routines to describe Boscarino's linear system in "ON A 
! CLASS OF UNIFORMLY ACCURATE IMEX RUNGE–KUTTA SCHEMES AND APPLICATIONS TO 
! HYPERBOLIC SYSTEMS WITH RELAXATION" pg 11, equations 31a/b
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
! SBP_COEF_MODULE.F90           *DEFINES CSR OPERATORS 
! JACOBIAN_CSR_MOD.F90          *ALLOCATE AND STORE CSR JACOBIAN VARIABLES
! UNARY_MOD.F90                 *PERFORMS SPARSE MATRIX OPERATIONS
! MATVEC_MODULE.F90             *PERFORMS SPARSE MATRIX*VECTOR OPERATIONS
!******************************************************************************
      module Boscarino31_Mod 

      use precision_vars, only: wp,two,pi

      implicit none; save 
  
      private
      public  :: Boscarino31

!--------------------------------VARIABLES-------------------------------------         
      integer,  parameter :: vecl=256 !must be even  
      integer,  parameter :: neq=2
      real(wp), parameter :: a=0.5_wp
      real(wp), parameter :: xL=-1.0_wp,xR=1.0_wp      
   
      real(wp), dimension(vecl/2) :: x      
      real(wp)                    :: dx

      integer,  dimension(4*vecl) :: jDeriv_comb_p
      real(wp), dimension(4*vecl) :: Deriv_comb_p
      real(wp), dimension(5*vecl) :: wrk_Jac
      integer,  dimension(vecl+1) :: iSource_p,iDeriv_comb_p,iwrk_Jac
      integer,  dimension(5*vecl) :: jwrk_Jac
      
      real(wp), dimension(vecl)   :: Source_p
      integer,  dimension(vecl)   :: jSource_p
      logical :: update_RHS,update_Jac
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
!******************************************************************************

      subroutine Boscarino31(nveclen,eq,ep,dt,tfinal,iDT,resE_vec,resI_vec,akk)

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
!------------------------------------------------------------------------------   
      
      Program_Step_Select: select case(programstep)
        !**Pre-initialization. Get problem name, vector length and grid**
        case('INITIALIZE_PROBLEM_INFORMATION')
          nveclen = vecl
          eq = neq
          probname='Boscar_31'     
          Jac_case='SPARSE'
          tol=1.0e-12_wp  
          dt_error_tol=5.0e-14_wp
          
          allocate(var_names(neq))
          var_names(:)=(/'Differential', 'Algebraic   '/)

          call grid()           

        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')
          
          !Update=true for Jac/RHS updates
          update_RHS=.true.; update_Jac=.true. !reset every new epsilon/dt
          
          !Allocate derivative operators
          if(.not. allocated(D1_per))call Define_CSR_Operators(vecl/2,dx)   
          
          !Time information
          tfinal = 0.2_wp  ! final time   
          choose_dt: select case(temporal_splitting)
            case('EXPLICIT'); dt = 0.000025_wp*0.1_wp/10**((iDT-1)/20.0_wp) ! explicit timestep
            case default    ; dt = 0.2_wp/10**((iDT-1)/20.0_wp)             ! implicit timestep      
          end select choose_dt
          
          ! Set IC's
          uvec(1:vecl:2)=sin(two*pi*x(:))
          uvec(2:vecl:2)=a*sin(two*pi*x(:))+ep*(a**2-1)*two*pi*cos(two*pi*x(:))
         
          !set exact solution at tfinal
          call exact_Bosc(uexact,ep) 
              
        case('BUILD_RHS')
          call Bosc_dUdt(uvec,dudt_D1,dudt_Source,ep)        
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
          update_RHS=.false. !no need to update RHS until next epsilon
          
        case('BUILD_JACOBIAN')
          choose_Jac_type: select case(temporal_splitting)
            case('IMPLICIT')
              if (update_Jac) then
                nnz_Jac=5*vecl+vecl/2
                call Allocate_Jac_CSR_Storage(vecl,nnz_Jac)
             
                call aplb(vecl,vecl,1,Source_p,jSource_p,iSource_p,    &
     &                    Deriv_comb_p,jDeriv_comb_p,iDeriv_comb_p,    &
     &                    wrk_Jac,jwrk_Jac,iwrk_Jac,5*vecl,iw,ierr) 
                if (ierr/=0) then; print*,'Build Jac ierr=',ierr; stop; endif 
              endif
              call Bosc_Jac(dt,akk,wrk_Jac,jwrk_Jac,iwrk_Jac)       
                        
            case('IMEX')
              nnz_Jac=vecl+vecl/2
              call Allocate_Jac_CSR_Storage(vecl,nnz_Jac)
              call Bosc_Jac(dt,akk,Source_p,jSource_p,iSource_p)
              
          end select choose_Jac_type
          update_Jac=.false.  !no need to update matrix that forms Jacobian until next epsilon/dt
          
      end select Program_Step_Select     
      end subroutine Boscarino31
!==============================================================================
!==============================================================================
!==============================================================================
! PRODUCES GRID
      subroutine grid()
      integer :: i
      
      do i=1,vecl/2
        x(i)= xL + (xR-xL)*(i-1.0_wp)/(vecl/2)
      enddo
      
      dx=x(2)-x(1)

      end subroutine grid
!==============================================================================
! RETURNS VECTOR WITH EXACT SOLUTION
      subroutine exact_Bosc(u,eps)

      real(wp),                  intent(in   ) :: eps
      real(wp), dimension(vecl), intent(  out) :: u

      integer                                  :: i
      real(wp), dimension(81,vecl+1) :: ExactTot
      real(wp) :: diff

      u(:)=0.0_wp
      
      !**Exact Solution** 
      open(unit=39,file='exact.Bosc_31_256.data')
      rewind(39)
      do i=1,81
        read(39,*)ExactTot(i,1:vecl)
        ExactTot(i,vecl+1) = 1.0_wp/10**((i-1)/(10.0_wp))  !  used for 81 values of ep
      enddo
      do i=1,81
        diff = abs(ExactTot(i,vecl+1) - eps)
        if(diff <= 1.0e-10_wp)then
          u(:) = ExactTot(i,:vecl)
          exit
        endif
      enddo

      end subroutine exact_Bosc

!==============================================================================
! DEFINES RHS 
      subroutine Bosc_dUdt(u,dudt_D1,dudt_source,eps)
      
      use SBP_Coef_Module, only: D1_per,jD1_per,iD1_per
      use matvec_module,   only: amux
      use unary_mod,       only: aplb, dperm, csort

      real(wp), dimension(vecl), intent(in   ) :: u
      real(wp), dimension(vecl), intent(  out) :: dudt_D1,dudt_source
      real(wp),                  intent(in   ) :: eps
           
      real(wp), dimension(vecl/2)  :: LR_Source,LL_Source
      integer,  dimension(vecl/2)  :: jLR_Source,jLL_Source
      
      real(wp), dimension(2*vecl)  :: UR_D1,LL_D1
      integer,  dimension(2*vecl)  :: UR_jD1,LL_jD1
       
      real(wp), dimension(vecl)    :: Source
      integer,  dimension(vecl)    :: jSource,j_perm
        
      real(wp), dimension(4*vecl)  :: Deriv_comb
      integer,  dimension(4*vecl)  :: jDeriv_comb
            
      integer,  dimension(5*vecl)  :: iw
      
      integer,  dimension(vecl+1)  :: iLR_Source,iLL_Source,UR_iD1,LL_iD1
      integer,  dimension(vecl+1)  :: iDeriv_comb,iSource
           
      integer,  dimension(10*vecl) :: iwork           
           
      integer,dimension(2) :: ierr
      integer              :: i

      ierr = 0
!------------------------------------------------------------------------------     
      if (update_RHS) then
        
! Set lower right diagonal
! | | |
! | |x|
        LR_Source(:)=-1.0_wp/eps
        iLR_Source(:vecl/2)   = 1
        iLR_Source(vecl/2+1:) = (/ (i, i=1, vecl/2+1) /)
        jLR_Source(:)=(/ (i, i=vecl/2+1, vecl) /)

! Set lower left diagonal
! | | |
! |x| |
        LL_Source(:)=a/eps
        iLL_Source(:vecl/2)   = 1
        iLL_Source(vecl/2+1:) = (/ (i, i=1, vecl/2+1) /)
        jLL_Source(:)=(/ (i, i=1, vecl/2) /)      
      
! Set upper right D1
! | |x|
! | | |      
        UR_D1=(-1.0_wp)*D1_per
        UR_iD1(:vecl/2+1)=iD1_per
        UR_iD1(vecl/2+2:)=iD1_per(vecl/2+1)
        UR_jD1=jD1_per(:)+vecl/2
      
! Set lower left D1
! | | |
! |x| |      
        LL_D1=(-1.0_wp)*D1_per
        LL_iD1(:vecl/2)=1
        LL_iD1(vecl/2+1:)=iD1_per
        LL_jD1=jD1_per

! Add all matricies
      ! combine both diagonals
      ! | | |
      ! |x|x|
        call aplb(vecl,vecl,1,LL_Source,jLL_Source,iLL_Source,LR_Source, &
     &            jLR_Source,iLR_Source,Source,jSource,iSource,vecl,iw,ierr(1))
      
      ! combine both D1's
      ! | |x|
      ! |x| |
        call aplb(vecl,vecl,1,LL_D1,LL_jD1,LL_iD1,UR_D1,UR_jD1,UR_iD1,    &
     &            Deriv_comb,jDeriv_comb,iDeriv_comb,4*vecl,iw,ierr(2))   
          
! permute matrix into:
! |x| |
! | |x|
        j_perm(1)=1
        j_perm(vecl/2+1)=2
        do i=2,vecl
          if (i==vecl/2+1) cycle
          j_perm(i)=j_perm(i-1) + 2
        enddo
        
      ! permute source's matrix    
        call dperm(vecl,Source,jSource,iSource,Source_p,jSource_p,iSource_p, &
     &             j_perm,j_perm,1)
        call csort(vecl,Source_p,jSource_p,iSource_p,iwork,.true.)  
   
      ! permute D1's matrix   
        call dperm(vecl,Deriv_comb,jDeriv_comb,iDeriv_comb,Deriv_comb_p,      &
     &             jDeriv_comb_p,iDeriv_comb_p,j_perm,j_perm,1)
        call csort(vecl,Deriv_comb_p,jDeriv_comb_p,iDeriv_comb_p,iwork,.true.)           
   
      endif
      
      ! get dudt     
      call amux(vecl,u,dudt_source,Source_p,    jSource_p,    iSource_p    )
      call amux(vecl,u,dudt_D1,    Deriv_comb_p,jDeriv_comb_p,iDeriv_comb_p)                 

      if (sum(ierr)/=0) then ! Catch errors
        print*,'Error building dudt'
        write(*,*)'ierr',ierr
        stop
      endif
  
      end subroutine Bosc_dUdt
!==============================================================================
      subroutine Bosc_Jac(dt,akk,a,ja,ia)
      
      use unary_mod,         only: aplsca
      use Jacobian_CSR_Mod,  only: iaJac, jaJac,  aJac
     
      real(wp),                  intent(in) :: dt,akk
      real(wp), dimension(:),    intent(in) :: a
      integer,  dimension(:),    intent(in) :: ja,ia
      
      integer,  dimension(vecl)             :: iw
      real(wp), dimension(size(a)+vecl/2)   :: wrk
      integer,  dimension(size(a)+vecl/2)   :: jwrk
      integer,  dimension(vecl+1)           :: iwrk      
      
      integer                               :: nnz

      nnz = size(a)
      
      wrk(:nnz)=-akk*dt*a(:)   
      wrk(nnz+1:)=0.0_wp
      jwrk(:nnz)=ja(:)
      jwrk(nnz+1:)=0
      iwrk(:)=ia(:)

      call aplsca(vecl,wrk,jwrk,iwrk,1.0_wp,iw)      

      iaJac=iwrk
      jaJac=jwrk
      aJac = wrk
      
      end subroutine Bosc_Jac      
!==============================================================================
      end module Boscarino31_mod
