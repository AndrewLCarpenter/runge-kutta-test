!******************************************************************************
! Module containing global variables for test_cases.f90 as well as allocation
! and deallocation routines
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
! SBP_COEF_MODULE.F90       *DEFINES CSR DERIVATIVE OPERATORS
!******************************************************************************

      module control_variables
      
      use precision_vars, only : wp
      
      implicit none; save
      
      private
      public :: Temporal_splitting, probname, Jac_case, isamp,jmax,jactual
      public :: tol,dt_error_tol
      public :: uvec,uexact,b,usum,uveco
      public :: errvec,errvecT,tmpvec
      public :: resE,resI,error,errorP,xjac,errorL2
      public :: b1save,b1Psave,b1L2save,ustage,predvec 
      public :: allocate_vars,deallocate_vars
      public :: programstep,var_names
!--------------------------VARIABLES-------------------------------------------
      character(len=80), parameter :: Temporal_Splitting = 'IMEX' 
!      character(len=80), parameter :: Temporal_Splitting = 'IMPLICIT'  
!      character(len=80), parameter :: Temporal_Splitting = 'EXPLICIT'
      character(len=9)             :: probname
      character(len=6)             :: Jac_case='DENSE' !default value
      character(len=80)            :: programstep
      integer, parameter           :: isamp=71   
      integer, parameter           :: jmax=81     
      integer, parameter           :: jactual=81     
      
      real(wp) :: tol,dt_error_tol

      real(wp), dimension(:),   ALLOCATABLE :: uvec,uexact,b,usum,uveco
      real(wp), dimension(:),   allocatable :: errvec,errvecT,tmpvec
      real(wp), dimension(:,:), allocatable :: resE,resI,error,errorP,xjac
      real(wp), dimension(:,:), allocatable :: b1save,b1Psave,ustage,predvec
      real(wp), dimension(:,:), allocatable :: errorL2,b1L2save
      character(len=12), dimension(:), allocatable :: var_names
!------------------------------------------------------------------------------
 
      contains

!==============================================================================
!******************************************************************************
! Subroutine to allocate global variables
!******************************************************************************
! MODULE VARIABLES:
! uvec     -> Array containing variables       
! uexact   -> Array containing exact solution to variables              
! resE     -> Explicit RHS vector
! resI     -> Implicit RHS vector
! error    -> solution error
! errorP   -> solution error, predicted
! b        -> convergence rate
! errorL2  -> solution error, L2 norm
! usum     -> Array containing summation from test_cases
! xjac     -> matrix containing dense Jacobian
! ustage   -> stage value predictor
! predvec  -> stage value predictor
! uveco    -> newton iteration / usum storage (?)
! errvec   -> temporary error storage
! errvecT  -> error storage for each time step
! tmpvec   -> (uvec-uexact)
! b1save   -> convergence rate storage
! b1Psave  -> convergence rate storage, predicted
! b1L2save -> convergence rate storage, L2 norm
! isamp    -> number of dt's
! jmax     -> number of epsilon values
!******************************************************************************
! INPUTS:
! nveclen -> u-vector length,          integer
! neq     -> number of equations,      integer
! is      -> maximum number of stages, integer
!******************************************************************************
      subroutine allocate_vars(nveclen,neq,is)
      
      integer, intent(in) :: nveclen,neq,is
    
      !**ALLOCATE VARIABLES**
      !problemsub
      AllOCATE(uvec(nveclen),uexact(nveclen),resE(nveclen,is),resI(nveclen,is))
   
      !data outs
      ALLOCATE(error(isamp,nveclen),errorP(isamp,nveclen),b(nveclen*2+neq))
      ALLOCATE(errorL2(isamp,neq))

      !Newton_iteration
      ALLOCATE(usum(nveclen),xjac(nveclen,nveclen))       

      !internal
      ALLOCATE(ustage(nveclen,is),predvec(nveclen,is),uveco(nveclen))
      ALLOCATE(errvec(nveclen),errvecT(nveclen),tmpvec(nveclen))
      ALLOCATE(b1save(jmax,nveclen),b1Psave(jmax,nveclen),b1L2save(jmax,neq))

      end subroutine allocate_vars
     
!==============================================================================
!******************************************************************************
! Subroutine to deallocate global variables at end of problem loop
!******************************************************************************
! SBP_COEF_MODULE.F90       *DEFINES CSR OPERATORS 
! JACOBIAN_CSR_MOD.F90      *ALLOCATE AND STORE CSR JACOBIAN VARIABLES
!******************************************************************************
! GLOBAL VARIABLES:
! From SBP_Coef_Module:
!   Pmat,Pinv,iD1,jD1,D1,iD2,jD2,D2,D1_per,iD1_per,jD1_per,D2_per,jD2_per,iD2_per
! From Jacobian_CSR_Mod:
!   iaJac,jaJac,aJac,jUJac,jLUJac,aLUJac,iw
! 
! MODULE VARIABLES:
! uvec     -> Array containing variables       
! uexact   -> Array containing exact solution to variables              
! resE     -> Explicit RHS vector
! resI     -> Implicit RHS vector
! error    -> solution error
! errorP   -> solution error, predicted
! b        -> convergence rate
! errorL2  -> solution error, L2 norm
! usum     -> Array containing summation from test_cases
! xjac     -> matrix containing dense Jacobian
! ustage   -> stage value predictor
! predvec  -> stage value predictor
! uveco    -> newton iteration / usum storage (?)
! errvec   -> temporary error storage
! errvecT  -> error storage for each time step
! tmpvec   -> (uvec-uexact)
! b1save   -> convergence rate storage
! b1Psave  -> convergence rate storage, predicted
! b1L2save -> convergence rate storage, L2 norm
!******************************************************************************
! INPUTS:
! nveclen -> u-vector length,          integer
! neq     -> number of equations,      integer
! is      -> maximum number of stages, integer
!******************************************************************************
      subroutine deallocate_vars()
      
      use SBP_Coef_Module,  only: Pmat,Pinv,iD1,jD1,D1,iD2,jD2,D2,D1_per,&
     &                            iD1_per,jD1_per,D2_per,jD2_per,iD2_per
      use Jacobian_CSR_Mod, only: iaJac,jaJac,aJac,jUJac,jLUJac,aLUJac,iw
      
      logical :: open_logical
      integer :: io_int,i
      
      !**DEALLOCATE VARIABLES**
      !problemsub
      DEAllOCATE(uvec,uexact,resE,resI,var_names)
      
      !data out
      DEALLOCATE(error,errorP,b,errorL2)
      
      !Newton_iteration
      DEALLOCATE(usum,xjac)   
         
      !internal
      DEALLOCATE(ustage,predvec,uveco)
      DEALLOCATE(errvec,errvecT,tmpvec,b1save,b1Psave,b1L2save)
      
      !D CSR
      if(allocated(D1)) deallocate(Pmat,Pinv,jD1,jD2,iD1,iD2,D1,D2,D1_per,&
     &                              iD1_per,jD1_per,D2_per,jD2_per,iD2_per) 
     
      !Jacobian CSR
      if(allocated(iaJac)) deallocate(iaJac,jaJac,aJac,jUJac,jLUJac,aLUJac,iw)  
      
      !close all open files
     ! do i=10,900
       ! inquire(unit=i, opened=open_logical)!, iostat=io_int)
       ! if (open_logical) print*,'io_int',io_int
       ! printstop
       ! enddo                               
      
      end subroutine deallocate_vars
      
!==============================================================================
      
      end module control_variables      
