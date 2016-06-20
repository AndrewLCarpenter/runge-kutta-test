      subroutine Stage_Value_Predictor(L,nrk,ustage,predvec,  &
                       & bD,alpha)!,svpB,al3N,al3D,al4N,al4D)
      
!  The SVP routine is called after the newton iteration has converged.  (stage L)
!  Thus, it is predicting the starting guesses of the next stage (L+1)
      use precision_vars

      implicit none

      integer,   parameter                               :: is=9
      integer,   parameter                               :: ivarlen=4

      integer,                             intent(in   ) :: L, nrk
      real(wp),   dimension(ivarlen,is),   intent(in   ) :: ustage
      real(wp),   dimension(ivarlen,is),   intent(inout) :: predvec
      real(wp),   dimension(is,4),         intent(in   ) :: bD
      real(wp),   dimension(is,is),        intent(  out) :: alpha

      real(wp),   dimension(is)                          :: bint


!     real(wp),   dimension(is,is,0:is),   intent(  out) :: svpB
!     real(wp),   dimension(is,is),        intent(  out) :: al3N
!     real(wp),   dimension(is,is),        intent(  out) :: al3D
!     real(wp),   dimension(is,is),        intent(  out) :: al4N
!     real(wp),   dimension(is,is),        intent(  out) :: al4D
!     character(len=25)                                  :: casename

      integer                                            :: j,M

      M = L+1

      if(L <= nrk-1) then
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

      end subroutine Stage_Value_Predictor
