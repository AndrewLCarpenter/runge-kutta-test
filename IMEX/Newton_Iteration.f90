      subroutine Newton_Iteration(uvec,iprob,Temporal_Splitting,L,ep,dt,nveclen,iDT,&
     & resE,resI,aI,usum,icount,k)

      use precision_vars

      implicit none

      integer, parameter :: is=9

      real(wp), dimension(nveclen),    intent(inout) :: uvec
      integer,                         intent(in   ) :: iprob
      character(80),                   intent(in   ) :: Temporal_Splitting
      integer,                         intent(in   ) :: L
      real(wp),                        intent(in   ) :: ep,dt
      integer,                         intent(in   ) :: nveclen,iDT
      real(wp), dimension(nveclen,is), intent(inout) :: resE,resI
      real(wp),                        intent(in   ) :: aI!from aI(L,L)
      real(wp), dimension(nveclen),    intent(in   ) :: usum
      integer,                         intent(inout) :: icount,k
      
      character(len=9)                     :: probname
      integer                              :: i,j  !Do loop variables
      integer                              :: programStep !input to problemsub
      real(wp), dimension(nveclen)         :: uveciter !Storage of uvec
      real(wp), dimension(nveclen)         :: uexact !input to problemsub (not needed)
      real(wp)                             :: tfinal !output of problemsub (not needed)
      real(wp), dimension(nveclen)         :: Rnewton !Newton 
      real(wp), dimension(nveclen,nveclen) :: xjac,xjacinv !Jacobians
      real(wp)                             :: tmp !temp variable for accuracy
!============================================================
      
      do k = 1,20
        icount = icount + 1

        uveciter(:) = uvec(:) !store old uvec

        programStep=2
        call problemsub(iprob,programStep,probname,nveclen,Temporal_Splitting,uvec,ep,uexact,dt,&
     &                   tfinal,iDT,resE(1,L),resI(1,L),aI,xjac)

        Rnewton(:) =uvec(:)-aI*resI(:,L)-usum(:)
         
        !**GET INVERSE JACOBIAN**
        programStep=3
        call problemsub(iprob,programStep,probname,nveclen,Temporal_Splitting,uvec,ep,uexact,dt,&
     &                   tfinal,iDT,resE(1,L),resI(1,L),aI,xjac)

!!      call JACOB2(nveclen,xjac,xjacinv)
!!      call Jacobian(uvec,xjac,dt,ep,aI,iprob,nveclen,sigma,rho,beta)
        call Invert_Jacobian(nveclen,xjac,xjacinv)

!        call Jacobian_CSR(uvec,dt,ep,aI,iprob,nveclen,sigma,rho,beta)
!        call ilu0(nveclen,aJac,jaJac,aluJac,jluJac,iw,ierr)
!        Backsolve to get solution

       !**want to use this but there is a different trunkation/round off error compared to the explicit do loop
       ! uvec(:)=uvec(:)-matmul(xjacinv,Rnewton)
        do i = 1,nvecLen
          do j = 1,nvecLen
            uvec(i) = uvec(i) - xjacinv(i,j)*Rnewton(j)     !u^n+1=u^n-J^-1*F
          enddo
       enddo

       tmp = sum(abs(uvec(:)-uveciter(:))) !check accuracy of zeros
       if(tmp.lt.1.0e-12_wp) exit
      enddo

      return
      end subroutine 
