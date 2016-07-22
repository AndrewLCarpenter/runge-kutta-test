!*****************************************************************************
! Module containing routines to describe Boscarino's linear system in "ON A 
! CLASS OF UNIFORMLY ACCURATE IMEX RUNGE–KUTTA SCHEMES AND APPLICATIONS TO 
! HYPERBOLIC SYSTEMS WITH RELAXATION" pg 11, equations 31a/b. 31a was modified 
! to be of the form "u_t+v_x=(b*u+c*v)/epsilon"
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
!******************************************************************************

      module Boscarino31_Mod 

      use precision_vars, only: wp,two,pi

      implicit none; save
  
      private
      public  :: Boscarino31

!--------------------------------VARIABLES-------------------------------------         
      integer,  parameter :: vecl=256 !must be even  
      real(wp), parameter :: a=0.5_wp
      real(wp), parameter :: xL=-1.0_wp,xR=1.0_wp      
   
      real(wp), dimension(vecl/2) :: x      
      real(wp)                    :: dx

      integer,  dimension(:), allocatable :: jap_mat 
      real(wp), dimension(:), allocatable ::  ap_mat
      integer,  dimension(vecl+1)         :: iap_mat
      real(wp)                            :: ep_store=0.0_wp
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

      subroutine Boscarino31(nveclen,ep,dt,tfinal,iDT,resE_vec,resI_vec,akk)

      use control_variables, only: temporal_splitting,probname,Jac_case, &
     &                             tol,dt_error_tol,uvec,uexact,programstep
      use SBP_Coef_Module,   only: Define_CSR_Operators,D1_per
      use unary_mod,         only: aplsca
      use Jacobian_CSR_Mod,  only: Allocate_CSR_Storage,  iaJac, jaJac,  aJac

!-----------------------VARIABLES----------------------------------------------
      !INIT vars     
      real(wp), intent(in   ) :: ep
      real(wp), intent(inout) :: dt
      integer,  intent(  out) :: nveclen
      real(wp), intent(  out) :: tfinal
      integer,  intent(in   ) :: iDT

      integer                    :: i
      
      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE_vec,resI_vec
      
      !Jacob vars
      real(wp), intent(in   ) :: akk
      
      integer,  dimension(vecl)                :: iw
      real(wp), dimension(size(ap_mat)+vecl/2) :: wrk
      integer,  dimension(size(ap_mat)+vecl/2) :: jwrk
      integer,  dimension(vecl+1)              :: iwrk
      integer                                  :: nnz_Jac
      
!------------------------------------------------------------------------------   
      
      Program_Step_Select: select case(programstep)
        !**Pre-initialization. Get problem name, vector length and grid**
        case('INITIALIZE_PROBLEM_INFORMATION')
          nveclen = vecl
          probname='Boscar_31'     
          Jac_case='SPARSE'
          tol=1.0e-12_wp  
          dt_error_tol=5.0e-14_wp
          call grid()           

        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')

          !Allocate derivative operators
          if(.not. allocated(D1_per)) call Define_CSR_Operators(vecl/2,dx)   

          !Time information
          tfinal = 0.2_wp  ! final time   
          choose_dt: select case(temporal_splitting)
            case('EXPLICIT'); dt = 0.000025_wp*0.1_wp/10**((iDT-1)/20.0_wp) ! explicit timestep
            case default    ; dt = 0.2_wp/10**((iDT-1)/20.0_wp)             ! implicit timestep      
          end select choose_dt
          
          ! Set IC's
          uvec(1:vecl:2)=sin(two*pi*x(:))
          uvec(2:vecl:2)=a*sin(two*pi*x(:))+ep*(a**2 - 1)*two*pi*cos(two*pi*x(:))
         
          !set exact solution at tfinal
          call exact_Bosc(uexact,ep) 
              
        case('BUILD_RHS')
          choose_RHS_type: select case (Temporal_Splitting)       
            case('EXPLICIT')
                call Bosc_dUdt(uvec,resE_vec,ep,dt)
                resI_vec(:)=0.0_wp                         
            case('IMPLICIT')
                 call Bosc_dUdt(uvec,resI_vec,ep,dt)
                 resE_vec(:)=0.0_wp
            case('IMEX')
          end select choose_RHS_type
       
        case('BUILD_JACOBIAN')
          choose_Jac_type: select case(temporal_splitting)
            case('IMPLICIT')
              nnz_Jac=size(ap_mat)+vecl/2
              call Allocate_CSR_Storage(vecl,nnz_Jac)
                            
              wrk=-akk*dt*ap_mat   
              jwrk(:size(jap_mat))=jap_mat
              jwrk(size(jap_mat)+1:)=0
              iwrk=iap_mat

              call aplsca(vecl,wrk,jwrk,iwrk,1.0_wp,iw)      

              iaJac=iwrk
              jaJac=jwrk
              aJac = wrk
            case('IMEX')
            
          end select choose_Jac_type
          
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
      open(unit=39,file='exact.Bosc_31_512.data')
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
      subroutine Bosc_dUdt(u,dudt,eps,dt)
      
      use SBP_Coef_Module, only: D1_per,jD1_per,iD1_per
      use matvec_module,   only: amux
      use unary_mod,       only: aplb, dperm, csort

      real(wp), dimension(vecl), intent(in   ) :: u
      real(wp), dimension(vecl), intent(  out) :: dudt
      real(wp),                  intent(in   ) :: eps, dt
           
      real(wp), dimension(vecl/2)                :: LR_Diag,LL_Diag
      integer,  dimension(vecl/2)                :: jLR_Diag,jLL_Diag
      
      real(wp), dimension(size( D1_per))         :: UR_D1,LL_D1
      integer,  dimension(size(jD1_per))         :: UR_jD1,LL_jD1
       
      real(wp), dimension(vecl)                  :: Diags
      integer,  dimension(vecl)                  :: jDiags,j_perm
        
      real(wp), dimension(2*size( D1_per))       :: Deriv_comb
      integer,  dimension(2*size( D1_per))       :: jDeriv_comb
            
      real(wp), dimension(2*size(D1_per)+vecl)   :: a_mat
      integer,  dimension(2*size(D1_per)+vecl)   :: iw,ja_mat
      
      integer,  dimension(vecl+1)                :: iLR_Diag,iLL_Diag,UR_iD1,LL_iD1
      integer,  dimension(vecl+1)                :: iDiags,iDeriv_comb,ia_mat
           
      integer,  dimension(4*size(D1_per)+2*vecl) :: iwork           
           
      integer,dimension(3) :: ierr
      integer              :: i

!------------------------------------------------------------------------------
      if (.not. allocated(ap_mat)) then
        allocate(ap_mat(2*size(D1_per)+vecl),jap_mat(2*size(D1_per)+vecl))
      endif

      if (.not. abs(eps-ep_store)<=1.0e-12_wp) then !epsilon value has changed
        ep_store=eps
        
! Set lower right diagaonal
! | | |
! | |x|
        LR_Diag(:)=-1.0_wp/eps
        iLR_Diag(:vecl/2)   = 1
        iLR_Diag(vecl/2+1:) = (/ (i, i=1, vecl/2+1) /)
        jLR_Diag(:)=(/ (i, i=vecl/2+1, vecl) /)

! Set lower left diagonal
! | | |
! |x| |
        LL_Diag(:)=a/eps
        iLL_Diag(:vecl/2)   = 1
        iLL_Diag(vecl/2+1:) = (/ (i, i=1, vecl/2+1) /)
        jLL_Diag(:)=(/ (i, i=1, vecl/2) /)      
      
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
        call aplb(vecl,vecl,1,LL_Diag,jLL_Diag,iLL_Diag,LR_Diag,jLR_Diag, &
     &            iLR_Diag,Diags,jDiags,iDiags,vecl,iw,ierr(1))
      
      ! combine both D1's
      ! | |x|
      ! |x| |
        call aplb(vecl,vecl,1,LL_D1,LL_jD1,LL_iD1,UR_D1,UR_jD1,UR_iD1,    &
     &            Deriv_comb,jDeriv_comb,iDeriv_comb,2*size(D1_per),iw,ierr(2))   
       
      ! combine diagonals and D1's
      ! | |x|
      ! |x|x|
        call aplb(vecl,vecl,1,Diags,jDiags,iDiags,Deriv_comb,jDeriv_comb, &
     &          iDeriv_comb,a_mat,ja_mat,ia_mat,2*size(D1_per)+vecl,iw,ierr(3)) 
  
! permiate matrix into:
! |x| |
! | |x|
        j_perm(1)=1
        j_perm(vecl/2+1)=2
        do i=2,vecl
          if (i==vecl/2+1) cycle
          j_perm(i)=j_perm(i-1) + 2
        enddo

        call dperm(vecl,a_mat,ja_mat,ia_mat,ap_mat,jap_mat,iap_mat,j_perm,j_perm,1)
        call csort(vecl,ap_mat,jap_mat,iap_mat,iwork,.true.)     
      endif
      
 ! get dudt     
      call amux(vecl,u,dudt,ap_mat,jap_mat,iap_mat)
   
      dudt(:)=dt*dudt(:)      
      
      if (sum(ierr)/=0) then ! Catch errors
        print*,'Error building dudt'
        stop
      endif
  
      end subroutine Bosc_dUdt
      
!==============================================================================
      end module Boscarino31_mod
