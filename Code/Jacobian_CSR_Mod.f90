! Module to allocate and store Jacobian CSR variables 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
!******************************************************************************
      module Jacobian_CSR_Mod
      
      use precision_vars,  only: wp
      
      implicit none; save
      
      public :: iaJac,jaJac,aJac,jUJac,jLUJac,aLUJac,iw
      public :: Allocate_Jac_CSR_Storage
      private
      
      integer,  allocatable, dimension(:) ::  iaJac,jaJac,jUJac,jLUJac,iw 
      real(wp), allocatable, dimension(:) ::  aJac,aLUJac     
      
      contains
      
!==============================================================================
!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Boscarino problem 
!******************************************************************************
! MODULE VARIABLES:
! iaJac  -> ia matrix for global storage of Jacobian, integer,  dimension(nJac+1), allocated
! jaJac  -> ja matrix for global storage of Jacobian, integer,  dimension(nnz),    allocated
!  aJac  -> a matrix for global storage of Jacobian,  real(wp), dimension(nnz),    allocated
! iw     -> work column matrix for Jacobian ops.,     integer,  dimension(nJac),   allocated
!  jUJac -> matrix for storage of LU Jacobian,        integer,  dimension(nJac),   allocated & init
! jLUJac -> matrix for storage of LU Jacobian,        integer,  dimension(nnz*4),  allocated & init
! aLUJac -> matrix for storage of LU Jacobian,        real(wp), dimension(nnz*4),  allocated & init
!******************************************************************************
! INPUTS:
! nJac   -> row/column size of the Jacobian matrx     integer
! nnz    -> number of non zero terms in the Jacobian  integer
!******************************************************************************
      subroutine Allocate_Jac_CSR_Storage(nJac,nnz)
      
      integer, intent(in) :: nJac !nJac=nveclen
      integer, intent(in) :: nnz 
      
      if(.not. allocated(iaJac)) then
      
        allocate(iaJac(nJac+1))
        allocate(jaJac(nnz))
        allocate( aJac(nnz))

        allocate(jUJac(nJac))
        allocate(jLUJac(nnz*4))
        allocate(aLUJac(nnz*4))

        allocate(iw(nJac))
        
        jLUJac(:)=0.0_wp ; aLUJac(:)=0.0_wp; aJac=0.0_wp
        
      endif

      end subroutine Allocate_Jac_CSR_Storage

!==============================================================================
      end module Jacobian_CSR_Mod
