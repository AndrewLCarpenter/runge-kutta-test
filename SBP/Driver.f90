      program driver

      use SBP_Coef_Module,    only: D1_242, D2_242

      implicit none

      integer,   parameter                      :: wp = 8

      integer,   parameter                      ::   n     =  10
      integer,   parameter                      ::   order = 242

      real(wp),  dimension(n  )                 ::  Pmat, Pinv

!     CSR storage for derivative matrices
      integer                                   :: nnz_D1, nnz_D2
      integer,   dimension(n+1)                 :: iD1, iD2
      integer,   dimension(:), allocatable      :: jD1, jD2
      real(wp),  dimension(:), allocatable      ::  D1,  D2

      integer                                   :: i,j

      continue

      if(order == 242) then
        nnz_D1 = 28+4*(n-8)
        allocate(jD1(nnz_D1),D1(nnz_D1))
        nnz_D2 = 38+5*(n-8)
        allocate(jD2(nnz_D2),D2(nnz_D2))
      endif
     
      call D1_242(n,nnz_D1,iD1,jD1,D1,Pmat,Pinv)
      call D2_242(n,nnz_D2,iD2,jD2,D2,Pmat,Pinv)

      stop
      end program
