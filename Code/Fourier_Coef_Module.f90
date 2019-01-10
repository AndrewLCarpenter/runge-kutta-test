!******************************************************************************
! Module containing routines to create Fourier derivative operators
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90      *DEFINES PRECISION FOR ALL VARIABLES
!******************************************************************************
      module Fourier_Coef_Module

      use precision_vars, only: wp
      use blas_module,    only: toeplitz

      implicit none

      real(wp),  parameter                      :: tol=1.0e-12_wp
      real(wp), dimension(:,:), allocatable     :: dmat

      contains

!==============================================================================

      subroutine Fourier_Coefficients(n,dmat)
      
      integer,                               intent(in   )   :: n
      real(wp), dimension(:,:), allocatable, intent(  out)   :: dmat

      real(wp), dimension(:), allocatable       :: Ftmp
      real(wp), dimension(:), allocatable       :: nrow, ncol

      integer                                   :: i, j, nmid
!     storage for derivative matrices

      allocate(dmat(n,n))        ;  dmat(:,:) = 0.0_wp ;
      allocate(Ftmp((n-1)/2))    ;  Ftmp(:)   = 0.0_wp ;
      allocate(nrow(n),ncol(n))  ;  nrow(:)   = 0.0_wp ; ncol(:) = 0.0_wp ;
      
      select case(n)

        case(09)
          Ftmp = reshape((/+1.46190220008155e+00_wp,-7.77861913430207e-01_wp,+5.77350269189626e-01_wp, &
                         & -5.07713305942873e-01_wp/), (/4/)) ;
        case(17)
          Ftmp = reshape((/+2.72109557587590e+00_wp,-1.38411497565444e+00_wp,+9.49789992441034e-01_wp, &
                         & -7.42174904521600e-01_wp,+6.26552889972156e-01_wp,-5.58557309966455e-01_wp, &
                         & +5.19844738543910e-01_wp,-5.02142049457839e-01_wp/), (/8/)) ;  
        case(33)
          Ftmp = reshape((/+5.26005483312601e+00_wp,-2.64199055435430e+00_wp,+1.77473276644212e+00_wp, &
                         & -1.34530672222908e+00_wp,+1.09116337944630e+00_wp,-9.24828432956910e-01_wp, &
                         & +8.08853403647791e-01_wp,-7.24554712838073e-01_wp,+6.61594815223310e-01_wp, &
                         & -6.13816303736287e-01_wp,+5.77350269189631e-01_wp,-5.49672837535945e-01_wp, &
                         & +5.29100070706893e-01_wp,-5.14503033610254e-01_wp,+5.05141613269013e-01_wp, &
                         & -5.00566972596134e-01_wp/), (/16/)) ;
        case(65)
          Ftmp = reshape((/+1.03491000818112e+01_wp,-5.18059980635193e+00_wp,+3.46046985214936e+00_wp, &
                         & -2.60244908613881e+00_wp,+2.08929073443020e+00_wp,-1.74858386669814e+00_wp, &
                         & +1.50644221065292e+00_wp,-1.32592638435841e+00_wp,+1.18651800096810e+00_wp, &
                         & -1.07590933716754e+00_wp,+9.86271390601783e-01_wp,-9.12387935110402e-01_wp, &
                         & +8.50650808352034e-01_wp,-7.98485818517622e-01_wp,+7.54008317746921e-01_wp, &
                         & -7.15807966591665e-01_wp,+6.82809497124985e-01_wp,-6.54179924027914e-01_wp, &
                         & +6.29265096120166e-01_wp,-6.07545324666924e-01_wp,+5.88603735735435e-01_wp, &
                         & -5.72103305080104e-01_wp,+5.57769941082280e-01_wp,-5.45379860948039e-01_wp, &
                         & +5.34750068687560e-01_wp,-5.25731112119135e-01_wp,+5.18201542434539e-01_wp, &
                         & -5.12063667578899e-01_wp,+5.07240307091642e-01_wp,-5.03672338432850e-01_wp, &
                         & +5.01316884412492e-01_wp,-5.00146035600253e-01_wp/), (/32/)) ;

        case default

        write(*,*)'Not a valid dimension for Fourier method'
        write(*,*)'Stopping'

      end select

      nmid = (n-1)/2 ;
      do i = 1,nmid
         nrow(i+1) = +Ftmp(i) ; nrow(n+1-i) = -Ftmp(i) ;
      enddo
      ncol(:) = -nrow(:)

      call toeplitz(n,ncol,nrow,dmat)

      end subroutine Fourier_Coefficients

! ======================================================================================

      end module Fourier_Coef_Module
