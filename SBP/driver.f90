      program driver

      use SBP_Coef_Module,    only: D1_242

      implicit none

      integer,   parameter                      :: wp = 8

      integer,   parameter                      ::   n     = 100
      integer,   parameter                      ::   order = 242

      real(wp),  dimension(n  )                 ::  Pmat, Pinv

!     CSR storage for derivative matrices
      integer                                   :: nnz_D1
      integer,   dimension(n+1)                 :: iD1
      integer,   dimension(:), allocatable      :: jD1
      real(wp),  dimension(:), allocatable      ::  D1

      integer                                   :: i,j

      continue

      if(order == 242) then
        nnz_D1 = 28+4*(n-8)
        allocate(jD1(nnz_D1),D1(nnz_D1))
      endif
      call D1_242(n,nnz_D1,iD1,jD1,D1,Pmat,Pinv)

      stop
      end program
