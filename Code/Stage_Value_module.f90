!******************************************************************************
! Module to perform stage value prediction on uvec
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90         *ONTAINS VARIABLES USED IN THE PROGRAM
! RUNGE_KUTTA.F90               *CONTAINS RK CONSTANTS
!******************************************************************************

      module Stage_value_module
  
      use precision_vars
      use control_variables
      use runge_kutta

      implicit none
    
      public :: Stage_Value_Predictor,xnorm
      private
          
      contains
      
!==============================================================================      
      
      subroutine Stage_Value_Predictor(ipred,L,ktime)
      
      integer, intent(in  ) :: ipred,L,ktime
      
      if (ipred==2) then
        call Stage_Value_Predictor2(L,ktime)
      else 
        call Stage_Value_Predictor1(L)
      endif   
      
      end subroutine Stage_Value_Predictor
      
!==============================================================================
      
      subroutine Stage_Value_Predictor1(L)
      
!  The SVP routine is called after the newton iteration has converged.  (stage L)
!  Thus, it is predicting the starting guesses of the next stage (L+1)

      integer,                    intent(in   ) :: L
    !  real(wp),   dimension(:,:), intent(in   ) :: ustage
     ! real(wp),   dimension(:,:), intent(inout) :: predvec
      !real(wp),   dimension(is,4),         intent(in   ) :: bD
      !real(wp),   dimension(is,is),        intent(  out) :: alpha

      !real(wp),   dimension(is)                          :: bint


!     real(wp),   dimension(is,is,0:is),   intent(  out) :: svpB
!     real(wp),   dimension(is,is),        intent(  out) :: al3N
!     real(wp),   dimension(is,is),        intent(  out) :: al3D
!     real(wp),   dimension(is,is),        intent(  out) :: al4N
!     real(wp),   dimension(is,is),        intent(  out) :: al4D

      integer :: j,M

      M = L+1

      if(L <= ns-1) then
        predvec(:,M) = ustage(:,L) ! previous guess as starter
      else
        predvec(:,2) = ustage(:,6)
      endif

!     if(L .ge. 2)then
!       predvec(:,M)   = ustage(:,1)                 ! put predict into start guess
!       do j = 2,L
!         predvec(:,M) = predvec(:,M) + alpha(M,j)*(ustage(:,j)-ustage(:,1))
!       enddo
!     endif

!     if(L == 2)then
!       bint(1) =  9.9518675264213746_wp
!       bint(2) =  4.8366852488953721_wp
!       bint(3) =-24.163405114569394_wp
!       bint(4) = 14.152132944153401_wp
!       bint(5) =  0.94399768676237158_wp
!       predvec(:,M) = ustage(:,3) + bint(1)*(ustage(:,3)-ustage(:,3))&
!                    &             + bint(2)*(ustage(:,4)-ustage(:,3))&
!                    &             + bint(3)*(ustage(:,5)-ustage(:,3))&
!                    &             + bint(4)*(ustage(:,6)-ustage(:,3))&
!                    &             + bint(5)*(ustage(:,2)-ustage(:,3))
!     endif

      end subroutine Stage_Value_Predictor1
 
!==============================================================================

      subroutine Stage_Value_Predictor2(L,ktime)
      
      integer,                  intent(in  )  :: L,ktime
      !real(wp), dimension(:,:), intent(in   ) :: ustage
      !real(wp), dimension(:),   intent(inout) :: uvec
   !   real(wp), dimension(:),   intent(in   ) :: uveco
      !real(wp), dimension(:,:), intent(in   ) :: alpha      
            
      integer :: j
      real(wp), dimension(5) :: bint
  
        uvec(:)  = uveco(:)          ! put predict into start guess
        do j = 2,L-1
          uvec(:) = uvec(:) + alpha(L,j)*(ustage(:,j)-uvec(:))
        enddo
    !    print*,'uvec',uvec
       ! print*,'alpha',alpha(:,:)
     !   print*,'ustage',ustage(:,:)
      !  stop
                
        if(L == 2 .and.ktime/=1)then
          bint(1) =  9.9518675264213746_wp 
          bint(2) =  4.8366852488953721_wp 
          bint(3) =-24.163405114569394_wp 
          bint(4) = 14.152132944153401_wp
          bint(5) =  0.94399768676237158_wp
          uvec(:) = ustage(:,3) +bint(1)*(ustage(:,3)-ustage(:,3))&
     &                          +bint(2)*(ustage(:,4)-ustage(:,3))&
     &                          +bint(3)*(ustage(:,5)-ustage(:,3))&
     &                          +bint(4)*(ustage(:,6)-ustage(:,3))&
     &                          +bint(5)*(ustage(:,2)-ustage(:,3))
        endif    

      end subroutine Stage_Value_Predictor2

!==============================================================================
!  COMPUTES NORM OF TWO VECTORS
      function xnorm(vec1,vec2)
      
      real(wp)                           :: xnorm
      real(wp), dimension(:), intent(in) :: vec1,vec2
      
      xnorm = sqrt(dot_product(vec1(:)-vec2(:),vec1(:)-vec2(:))/  &
     &                          dot_product(vec1(:),vec1(:)))
      
      return
      end function xnorm
     
!==============================================================================
      end module Stage_value_module
