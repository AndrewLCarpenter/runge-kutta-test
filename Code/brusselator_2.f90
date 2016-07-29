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
! u'=1-4*u+u*v*v+ep*d^2u/dx^2
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
      
      integer,  parameter :: vecl=8     ! total uvec length == 2*vecl
      real(wp), parameter :: xL=0.0_wp, xR=1.0_wp
      real(wp), dimension(vecl) :: x
      real(wp) ::dx
      
      integer,  dimension(10*vecl) :: jDeriv2_p
      real(wp), dimension(10*vecl) :: Deriv2_p
      real(wp), dimension(5*vecl) :: wrk_Jac
      integer,  dimension(vecl*2+1) :: iSource_p,iDeriv2_p,iwrk_Jac
      integer,  dimension(5*vecl) :: jwrk_Jac
      
      real(wp), dimension(4*vecl)   :: Source_p
      integer,  dimension(4*vecl)   :: jSource_p
      
      contains
      
!==============================================================================      
      subroutine Brusselator(nveclen,ep,dt,tfinal,iDT,resE_vec,resI_vec,akk)

      use control_variables, only: temporal_splitting,probname,Jac_case, &
     &                             tol,dt_error_tol,uvec,uexact,programstep
      use SBP_Coef_Module,   only: Define_CSR_Operators,D2_per,nnz_D2_per
      use unary_mod,         only: aplb
      use Jacobian_CSR_Mod,  only: Allocate_Jac_CSR_Storage
!-----------------------VARIABLES----------------------------------------------

      !INIT vars
      real(wp), intent(in   ) :: ep
      real(wp), intent(inout) :: dt
      integer,  intent(  out) :: nveclen
      real(wp), intent(  out) :: tfinal
      integer,  intent(in   ) :: iDT

!      real(wp), dimension(81,vecl+1) :: ExactTot
!      real(wp)                       :: diff
!      integer                        :: i

      !RHS vars
      real(wp), dimension(vecl*2), intent(  out) :: resE_vec,resI_vec
      real(wp), dimension(vecl*2)                :: dudt_Dx,dudt_Source
            
      !Jacob vars
      real(wp), intent(in   ) :: akk
      integer :: ierr
      integer :: nnz_Jac
      integer, dimension(vecl*2) :: iw 
!------------------------------------------------------------------------------
      
      Program_Step_Select: select case(programstep)
        !**Pre-initialization. Get problem name and vector length**
        case('INITIALIZE_PROBLEM_INFORMATION')
          nvecLen = 2*vecl
          probname='Brusselat'     
          tol=9.9e-11_wp  
          dt_error_tol=1.0e-13_wp
          Jac_case='SPARSE'
          
          call grid()

        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')
          
          !Allocate derivative operators
          if(.not. allocated(D2_per))call Define_CSR_Operators(vecl,dx) 
           
          !Time information
          dt = 0.25_wp/10**((iDT-1)/20.0_wp) ! timestep   
          tfinal = 1.0_wp                   ! final time
 
          uexact(:)=0.0_wp         
          !**Exact Solution** 
!          open(unit=39,file='exact.brusselator4.data')
!          rewind(39)
!          do i=1,81
!            read(39,*)ExactTot(i,1),ExactTot(i,2)
!            ExactTot(i,3) = 1.0_wp/10**((i-1)/(10.0_wp))  !  used for 81 values of ep
!          enddo
!          do i=1,81
!            diff = abs(ExactTot(i,3) - ep)
!            if(diff.le.1.0e-10_wp)then
!              uexact(:) = ExactTot(i,:vecl)
!              exit
!            endif
!          enddo

          !**IC**
          uvec(1:2*vecl:2) = 1.0_wp+sin(two*pi*x(:))
          uvec(2:2*vecl:2) = 3.0_wp
      
          !**Exact Solution**
          uexact(:)=0.0_wp
          
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

              nnz_Jac=12*vecl
              call Allocate_Jac_CSR_Storage(vecl*2,nnz_Jac)
             
!              call aplb(vecl*2,vecl*2,1,Source_p,jSource_p,iSource_p,    &
!     &                  Deriv2_p,jDeriv2_p,iDeriv2_p,    &
!     &                  wrk_Jac,jwrk_Jac,iwrk_Jac,12*vecl,iw,ierr) 
!              if (ierr/=0) then; print*,'Build Jac ierr=',ierr; stop; endif 

              call Build_Jac(uvec,dt,akk,Deriv2_p,jDeriv2_p,iDeriv2_p)       
                        
            case('IMEX')
              nnz_Jac=vecl+vecl/2
              call Allocate_Jac_CSR_Storage(vecl*2,nnz_Jac)
              call Build_Jac(uvec,dt,akk,Source_p,jSource_p,iSource_p)
              
          end select choose_Jac_type
        
      end select Program_Step_select
      end subroutine Brusselator
!==============================================================================
!==============================================================================
!==============================================================================
! PRODUCES GRID
      subroutine grid()
      integer :: i
      
      do i=1,vecl
        x(i)= xL + (xR-xL)*(i-1.0_wp)/(vecl-1.0_wp)
      enddo

      dx=x(2)-x(1)      
      
      end subroutine grid
!==============================================================================
! DEFINES RHS 
      subroutine Bruss_dUdt(vec,dudt_deriv2,dudt_source,eps)
      
      use SBP_Coef_Module, only: D2_per,jD2_per,iD2_per
      use unary_mod,       only: aplb,dperm
      use matvec_module,   only: amux

      real(wp), dimension(:),      intent(in   ) :: vec
      real(wp), dimension(vecl*2), intent(  out) :: dudt_deriv2,dudt_source
      real(wp),                    intent(in   ) :: eps

      !real(wp), dimension(vecl) :: dudt_s,dvdt_s,d2udx2,d2vdx2,u,v
      
      integer :: i
      integer, dimension(4) :: ierr=0
      integer, dimension(vecl*2) :: iw,j_perm
      
      real(wp), dimension(vecl) :: u,v
!      integer,  dimension(vecl*2+1) :: iLR_Source,iLL_Source
!      integer,  dimension(vecl*2+1) :: iUR_Source,iUL_Source
      integer,  dimension(vecl*2+1) :: iUL_D2,iLR_D2
            
!      real(wp), dimension(vecl) :: LR_Source,LL_Source,UR_Source,UL_Source
!      integer,  dimension(vecl) :: jLR_Source,jLL_Source,jUR_Source,jUL_Source
      
      real(wp), dimension(5*vecl) :: UL_D2,LR_D2
      integer,  dimension(5*vecl) :: jUL_D2,jLR_D2
      
!      real(wp), dimension(2*vecl) :: U_Source,L_Source
!      integer,  dimension(2*vecl) :: jU_Source,jL_Source
!      integer,  dimension(vecl*2+1) :: iU_Source,iL_Source
      
!      real(wp), dimension(4*vecl) :: Source
!      integer,  dimension(4*vecl) :: jSource
!      integer,  dimension(vecl*2+1) :: iSource
      
      real(wp), dimension(10*vecl) :: Deriv2
      integer,  dimension(10*vecl) :: jDeriv2
      integer,  dimension(vecl*2+1) :: iDeriv2
      
      
      u=vec(1:vecl*2:2); v=vec(2:vecl*2:2)

! Set upper left D2
! |x| |
! | | |      
        UL_D2(:)=eps*D2_per(:)
        iUL_D2(:vecl)=iD2_per(:)
        iUL_D2(vecl+1:)=iD2_per(vecl+1)
        jUL_D2(:)=jD2_per(:)
         
! Set lower right D2
! | | |
! | |x|      
        LR_D2(:)=eps*D2_per(:)
        iLR_D2(:vecl)=1
        iLR_D2(vecl+1:)=iD2_per(:)
        jLR_D2(:)=jD2_per(:)+vecl

! Add all matricies   
! combine both D2's
! |x| |
! | |x|
        call aplb(vecl*2,vecl*2,1,LR_D2,jLR_D2,iLR_D2,UL_D2,jUL_D2,iUL_D2,    &
     &            Deriv2,jDeriv2,iDeriv2,vecl*10,iw,ierr(4))   
      
! permute matrix into:
! |x| |
! | |x|
        j_perm(1)=1
        j_perm(vecl+1)=2
        do i=2,vecl*2
          if (i==vecl+1) cycle
          j_perm(i)=j_perm(i-1) + 2
        enddo
 
      ! permute source's matrix    
!        call dperm(vecl*2,Source,jSource,iSource,Source_p,jSource_p,iSource_p,&
!     &             j_perm,j_perm,1)
      ! permute D2's matrix   
        call dperm(vecl*2,Deriv2,jDeriv2,iDeriv2,Deriv2_p,      &
     &             jDeriv2_p,iDeriv2_p,j_perm,j_perm,1)

      ! get dudt     
!      call amux(vecl*2,vec,dudt_source,Source_p,jSource_p,iSource_p)
      call amux(vecl*2,vec,dudt_Deriv2,Deriv2_p,jDeriv2_p,iDeriv2_p)   
      
      dudt_Source(1:2*vecl:2)=1.0_wp-4.0_wp*u(:)+u(:)*v(:)*v(:)
      dudt_Source(2:2*vecl:2)=       3.0_wp*u(:)-u(:)*u(:)*v(:)

      if (sum(ierr)/=0) then ! Catch errors
        print*,'Error building dudt'
        write(*,*)'ierr',ierr
        stop
      endif
    
      end subroutine Bruss_dUdt
!==============================================================================
!  CREATES JACOBIAN
      subroutine Build_Jac(vec,dt,akk,a,ja,ia)

  !    use SBP_Coef_Module,  only: nnz_D2_per,D2_per!,Pinv,D1_per,jD1_per,iD1_per,nnz_D1_per,jD2_per,iD2_per       
      use matvec_module,    only: amux
      use unary_mod,        only: aplb,aplsca,dperm!,amudia,diamua,apldia
      use Jacobian_CSR_Mod, only: iaJac,jaJac,aJac
      
      real(wp), dimension(:), intent(in) :: vec
      real(wp),               intent(in) :: dt,akk
      real(wp), dimension(:), optional, intent(in) :: a
      integer,  dimension(:), optional, intent(in) :: ja,ia
      
      real(wp), dimension(vecl) :: u,v
      integer :: ierr=0
      integer, dimension(2) :: jArray
      integer,dimension(4*vecl)::jSource
      real(wp), dimension(4*vecl)::Source
      integer,dimension(1+vecl*2)::iSource
      integer :: i,nnz
      
      integer,dimension(:),allocatable::jwrk
      integer,dimension(vecl*2+1)::iwrk
      real(wp), dimension(:),allocatable::wrk
      
      integer,dimension(vecl*2)::iw

          integer, dimension(vecl*2)::j_perm
       
       
      u=vec(1:vecl*2:2); v=vec(2:vecl*2:2)
      
      jArray=(/0,vecl/)
      
 !     if(present(ia)) then
 !       nnz=12*vecl
 !     else
 !       nnz=4*vecl
 !     endif
     
! Set Source Jacobian
      Source(1:vecl*2:2)=-4.0_wp+v(:)*v(:)
      Source(2:vecl*2:2)=2.0_wp*u(:)*v(:)
      Source(vecl*2+1:vecl*4:2)=3.0_wp-2.0_wp*u(:)*v(:)
      Source(vecl*2+2:vecl*4:2)=-u(:)*u(:)
      iSource(:) = (/ (i, i=1, vecl*4+1,2) /)
      do i=1,vecl
        jSource(2*i-1:2*i)=jArray(:)+i
      enddo
      jSource(vecl*2+1:)=jSource(:vecl*2)
    
      j_perm(1)=1
      j_perm(vecl+1)=2
      do i=2,vecl*2
        if (i==vecl+1) cycle
        j_perm(i)=j_perm(i-1) + 2
      enddo

      ! permute source's matrix    
      call dperm(vecl*2,Source,jSource,iSource,Source_p,jSource_p,iSource_p,&
     &             j_perm,j_perm,1)

      if(present(ia)) then
        nnz=12*vecl
        allocate(wrk(nnz),jwrk(nnz))
        call aplb(vecl*2,vecl*2,1,Source_p,jSource_p,iSource_p,a,ja,ia,    &
     &            wrk,jwrk,iwrk,nnz,iw,ierr)
      else
        nnz=4*vecl
        allocate(wrk(nnz),jwrk(nnz))
        wrk(:)=Source_p(:)
        jwrk(:)=jSource_p(:)
        iwrk(:)=iSource_p(:)
      endif

      wrk(:)=-akk*dt*wrk(:)   
      call aplsca(vecl*2,wrk,jwrk,iwrk,1.0_wp,iw)      

      iaJac=iwrk
      jaJac=jwrk
      aJac = wrk

      deallocate(wrk,jwrk)
      
!-----------------------------------------------------------------------------  
     ! if (sum(ierr)/=0) then ! Catch errors
      if (ierr/=0) then ! Catch errors
        print*,'Error building Jacobian'
        print*,'ierr=',ierr
        stop
      endif

      end subroutine Build_Jac
      
      end module Brusselator_mod
