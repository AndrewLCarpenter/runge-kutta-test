!******************************************************************************
! Module to perform QR decomposition
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
!******************************************************************************

      module QR_Module

      use precision_vars, only : wp

      implicit none

      public    ::  qrdcmp, qrsolv, qrsolv_MRHS
      private
      
      real(wp), parameter   :: tolerance = 1.0e-17_wp
      contains

!==============================================================================

      SUBROUTINE qrdcmp(a,n,np,c,d,sing)

      INTEGER,                     intent(in   ) :: n,np
      REAL(wp),  dimension(np,np), intent(inout) :: a
      REAL(wp),  dimension(n),     intent(  out) :: c,d
      LOGICAL,                     intent(  out) :: sing

      !  Constructs the QR decomposition of a(1:n,1:n) , with physical dimension np .  The upper
      !  triangular  matrix R is returned in  the upper triangle  of a , except for the diagonal elements
      !  of R which are returned in d(1:n) .  The orthogonal matrix Q is represented as a product of
      !  n - 1 Householder matrices Q 1 ::: Q n - 1 ,where Q j = 1 - u j ?  u j =c j .The i th component
      !  of u j is zero for i =1 ;:::;j - 1 while the nonzero components are returned in a(i,j) for
      !  i = j;:::;n .  sing returns  as  true  if  singularity  is  encountered during  the decomposition,
      !  but  the  decomposition  is  still  completed  in  this  case.

      INTEGER                                    :: i,j,k
      REAL(wp)                                   :: scl,sigma,wrk,tau

      sing=.false.
      do k=1,n-1
        scl=0.0_wp
        do i=k,n
          scl=max(scl,abs(a(i,k)))
        enddo
        if(abs(scl) <= tolerance)then
          sing=.true.
          c(k)=0.0_wp
          d(k)=0.0_wp
        else
          do i=k,n
            a(i,k)=a(i,k)/scl
          enddo
          wrk=0.0_wp
          do i=k,n
            wrk=wrk+a(i,k)**2
          enddo
          sigma=sign(sqrt(wrk),a(k,k))
          a(k,k)=a(k,k)+sigma
          c(k)=sigma*a(k,k)
          d(k)=-scl*sigma
          do j=k+1,n
            wrk=0.0_wp
            do i=k,n
              wrk=wrk+a(i,k)*a(i,j)
            enddo
            tau=wrk/c(k)
            do i=k,n
              a(i,j)=a(i,j)-tau*a(i,k)
            enddo
          enddo
        endif
      enddo
      d(n)=a(n,n)
      if(abs(d(n)) <= tolerance)sing=.true.
      
      END subroutine qrdcmp

!==============================================================================

      SUBROUTINE qrsolv(a,n,np,c,d,b)

      implicit none

      INTEGER,                     intent(in   ) :: n,np
      REAL(wp),  dimension(np,np), intent(inout) :: a
      REAL(wp),  dimension(n),     intent(in   ) :: c,d
      REAL(wp),  dimension(n),     intent(inout) :: b

!     USES rsolv
!      Solves the set of n linear equations A x = b ,where a is a matrix with physical dimension np.
!      a , c ,and d are input  as  the output of the routine qrdcmp and  are not modi ed.  b(1:n)
!      is input as the right-hand side vector, and is overwritten with the solution vector on output.

      INTEGER        ::  i,j
      REAL(wp)       ::  wrk,tau

      do j=1,n-1                            !  Form Q^T . b
        wrk=0.0_wp
        do i=j,n
          wrk=wrk+a(i,j)*b(i)
        enddo
        tau=wrk/c(j)
        do i=j,n
          b(i)=b(i)-tau*a(i,j)
        enddo
      enddo
      call rsolv(a,n,np,d,b)                   !  Solve R x = Q^T b
      
      END subroutine qrsolv

!==============================================================================

      SUBROUTINE rsolv(a,n,np,d,b)

      INTEGER,                     intent(in   ) :: n,np
      REAL(wp),  dimension(np,np), intent(inout) :: a
      REAL(wp),  dimension(n),     intent(in   ) :: d
      REAL(wp),  dimension(n),     intent(inout) :: b

!     Solves the set of n linear  equations R x = b ,where R is an upper triangular  matrix stored
!     in a and d .  a and d are input  as  the output  of the routine qrdcmp and  are not modified.
!     b(1:n) is  input  as  the  right-hand  side  vector, and  is  overwritten  with the  solution  vector
!     on  output.

      INTEGER                                    :: i,j
      REAL(wp)                                   :: wrk

      b(n)=b(n)/d(n)
      do i=n-1,1,-1
        wrk=0.0_wp
        do j=i+1,n
          wrk=wrk+a(i,j)*b(j)
        enddo
        b(i)=(b(i)-wrk)/d(i)
      enddo
      
      end subroutine rsolv

!==============================================================================

      SUBROUTINE qrsolv_MRHS(a,n,np,m,c,d,b)

      implicit none

      INTEGER,                     intent(in   ) :: n,np,m
      REAL(wp),  dimension(np,np), intent(inout) :: a
      REAL(wp),  dimension(n),     intent(in   ) :: c,d
      REAL(wp),  dimension(n,m),   intent(inout) :: b

!     USES rsolv
!      Solves the set of n linear equations A x = b ,where a is a matrix with physical dimension np.
!      a , c ,and d are input  as  the output of the routine qrdcmp and  are not modi ed.  b(1:n,k)
!      is input as the right-hand side vector, and is overwritten with the solution vector on output.

      INTEGER                      ::  i,j,k
      REAL(wp)                     ::  wrk,tau

      do k=1,m                                !  loop over RHS vectors
        do j=1,n-1                            !  Form Q^T . b
          wrk = 0.0_wp
          do i=j,n
            wrk=wrk+a(i,j)*b(i,k)
          enddo
          tau=wrk/c(j)
          do i=j,n
            b(i,k)=b(i,k)-tau*a(i,j)
          enddo
        enddo
        call rsolv(a,n,np,d,b(:,k))                   !  Solve R x = Q^T b
      enddo
      
      END subroutine qrsolv_MRHS

!==============================================================================

      end module QR_Module
