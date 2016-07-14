!******************************************************************************
! Module containing global variables for test_cases.f90 as well as allocation
! and deallocation routines
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
!******************************************************************************

      module control_variables
      
      use precision_vars, only : wp
      
      implicit none; save
      
      private
      public :: Temporal_splitting, probname, Jac_case, isamp,jmax,jactual
      public :: tol,dt_error_tol
      public :: uvec,uexact,b,usum,uveco
      public :: errvec,errvecT,tmpvec
      public :: resE,resI,error,errorP,xjac
      public :: b1save,b1Psave,ustage,predvec 
      public :: allocate_vars,deallocate_vars
!--------------------------VARIABLES-------------------------------------------
!      character(len=80) :: Temporal_Splitting = 'IMEX'
      character(len=80)  :: Temporal_Splitting = 'IMPLICIT'
!      character(len=80) :: Temporal_Splitting = 'EXPLICIT'
      character(len=9)   :: probname
      character(len=6)   :: Jac_case='DENSE' !default value
      integer, parameter :: isamp=71   
      integer, parameter :: jmax=81     
      integer, parameter :: jactual=81     
      
      real(wp) :: tol,dt_error_tol

      real(wp), dimension(:),   ALLOCATABLE :: uvec,uexact,b,usum,uveco
      real(wp), dimension(:),   allocatable :: errvec,errvecT,tmpvec
      real(wp), dimension(:,:), allocatable :: resE,resI,error,errorP,xjac
      real(wp), dimension(:,:), allocatable :: b1save,b1Psave,ustage,predvec
      
!------------------------------------------------------------------------------
 
      contains

!==============================================================================
! ALLOCATE GLOBAL VARIABLES
      subroutine allocate_vars(nveclen,is)
      
      integer, intent(in) :: nveclen,is
      
      !**ALLOCATE VARIABLES**
      !problemsub
      AllOCATE(uvec(nveclen),uexact(nveclen),resE(nveclen,is),resI(nveclen,is))
   
      !data outs
      ALLOCATE(error(isamp,nveclen),errorP(isamp,nveclen),b(nveclen*2))

      !Newton_iteration
      ALLOCATE(usum(nveclen),xjac(nveclen,nveclen))       

      !internal
      ALLOCATE(ustage(nveclen,is),predvec(nveclen,is),uveco(nveclen))
      ALLOCATE(errvec(nveclen),errvecT(nveclen),tmpvec(nveclen))
      ALLOCATE(b1save(jmax,nveclen),b1Psave(jmax,nveclen))
       
       end subroutine allocate_vars
     
!==============================================================================
! DEALLOCATE GLOBAL VARIABLES
      subroutine deallocate_vars()
      
      use SBP_Coef_Module, only: Pmat,Pinv,iD1,jD1,D1,iD2,jD2,D2
      
      !**DEALLOCATE VARIABLES**
      !problemsub
      DEAllOCATE(uvec,uexact,resE,resI)
      
      !data out
      DEALLOCATE(error,errorP,b)
      
      !Newton_iteration
      DEALLOCATE(usum,xjac)   
         
      !internal
      DEALLOCATE(ustage,predvec,uveco)
      DEALLOCATE(errvec,errvecT,tmpvec,b1save,b1Psave)
      
      !CSR
      if(allocated(D1)) deallocate(Pmat,Pinv,jD1,jD2,iD1,iD2,D1,D2) 
      
      end subroutine deallocate_vars
      
!==============================================================================
      
      end module control_variables      
