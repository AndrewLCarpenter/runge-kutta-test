      module CSR_Variables

      use precision_vars

      implicit none

      real(wp),    parameter              ::   toljac = 1.0e-13_wp

      integer                             ::   nJac
      integer                             :: nnzJac
      integer                             ::   ierr

      integer,  allocatable, dimension(:) ::  iaJac
      integer,  allocatable, dimension(:) ::  jaJac
      real(wp), allocatable, dimension(:) ::   aJac

      integer,  allocatable, dimension(:) ::   juJac
      integer,  allocatable, dimension(:) ::  jLUJac
      real(wp), allocatable, dimension(:) ::  aLUJac

      integer,  allocatable, dimension(:) ::  iw

      end module CSR_Variables

