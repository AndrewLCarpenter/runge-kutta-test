      program driver

      use SBP_Coef_Module,    only: D1_242

      implicit none

      integer,   parameter                      :: wp = 8

      integer,   parameter                      ::   n = 9
      integer,   parameter                      :: nnz = 10+4*n

      integer,   dimension(n+1)                 :: ia
      integer,   dimension(nnz)                 :: ja

      real(wp),  dimension(nnz)                 ::  a

      real(wp),  dimension(n  )                 ::  Pmat, Pinv

      integer                                   :: i,j

      continue

      call D1_242(n,nnz,ia,ja,a,Pmat,Pinv)

      write(*,*)'ia'
      write(*,*) ia 
      write(*,*)'ja'
      write(*,*) ja 
      write(*,*)' a'
      write(*,*)  a 

      stop
      end program
