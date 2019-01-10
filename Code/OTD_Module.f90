! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
! MATVEC_MODULE.F90             *PERFORMS SPARSE MATRIX*VECTOR OPERATIONS
!******************************************************************************

      module OTD_Module

      use precision_vars,    only: wp

      implicit none
  
      private
      public  :: Orthogonalize_SubSpace

      contains

!==============================================================================

      subroutine Orthogonalize_SubSpace(vecL,OTDN,V)

!-----------------------VARIABLES----------------------------------------------
      !INIT vars     
      integer,                        intent(in   ) :: vecL,OTDN

      real(wp), dimension(vecL,OTDN), intent(inout) :: V
      
      integer                              :: i,j
      real(wp)                             :: t1, rij

      continue

      t1 = sqrt(dot_product(V(:,1),V(:,1))) ; V(:,1) = V(:,1) / t1

      if(OTDN == 1) return

      do i = 2,OTDN

        do j = 1,i-1

           rij = dot_product(V(:,i),V(:,j))

           V(:,i) = V(:,i) - V(:,j) * rij

        enddo
    
        t1 = sqrt(dot_product(V(:,i),V(:,i))) ; V(:,i) = V(:,i) / t1

      enddo

      end subroutine Orthogonalize_SubSpace      

!==============================================================================

      end module OTD_Module      

