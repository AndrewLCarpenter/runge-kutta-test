      module Jacobian_CSR_Mod
      
      use precision_vars
      use SBP_Coef_Module, only:nnz_D2
      
      implicit none
      
      public :: iaJac,jaJac,aJac,jUJac,jLUJac,aLUJac,iw,Allocate_CSR_Storage
      public :: Jacobian_CSR
      private
      
      real(wp),    parameter              ::   toljac = 1.0e-13_wp

      integer                             ::   nJac
      integer                             :: nnzJac
      logical                             :: allo_test=.false.

      integer,  allocatable, dimension(:) ::  iaJac
      integer,  allocatable, dimension(:) ::  jaJac
      real(wp), allocatable, dimension(:) ::   aJac

      integer,  allocatable, dimension(:) ::  jUJac
      integer,  allocatable, dimension(:) ::  jLUJac
      real(wp), allocatable, dimension(:) ::  aLUJac
!      integer                             ::   ierr
      integer,  allocatable, dimension(:) ::  iw         
      
      contains
      
!==============================================================================
      subroutine Allocate_CSR_Storage(nJac)

      integer, intent(in)  :: nJac !nJac=nveclen
      if(.not.allo_test) then

        allocate(iaJac(nJac+1))
        allocate(jaJac(nnz_D2))
        allocate( aJac(nnz_D2))

        allocate(jUJac(nJac))
        allocate(jLUJac(nnz_D2))
        allocate(aLUJac(nnz_D2))

        allocate(iw(nJac))
        
      endif
      allo_test=.true.
      end subroutine

!==============================================================================
      subroutine Jacobian_CSR(nJac,xjac)
 
      integer,                  intent(in   ) :: nJac
      real(wp), dimension(:,:), intent(in   ) :: xjac
      integer                                 :: i,j,icnt,jcnt

!      print*,'calling2'
      nnz_D2=nJac**2
      call Allocate_CSR_Storage(nJac)
 !     print*,'called2'
 !     print*,size(iaJac),size(jaJac),size(aJac)
!     U_t = F(U);  Jac = \frac{\partial F(U)}{\partial U};  xjac = I - akk dt Jac

      ! Initialize CSR
      iaJac(:) = 0
      iaJac(1) = 1
      jaJac(:) = 0 
      aJac(:) = 0.0_wp
      
      ! Store dense matrix into CSR format
      icnt = 0   
      do i = 1,nJac
        jcnt = 0   
        do j = 1,nJac
          if(abs(xjac(i,j)) >= tolJac) then
            icnt = icnt + 1 
            jcnt = jcnt + 1 
             
            jaJac(icnt) = j
             aJac(icnt) = xjac(i,j)
          endif
        enddo
        iaJac(i+1) = iaJac(i) + jcnt
      enddo

      end subroutine Jacobian_CSR
!==============================================================================
      end module Jacobian_CSR_Mod
