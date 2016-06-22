      subroutine Newton_Iteration(uvec,iprob,L,ep,dt,nveclen,iDT,&
     & resE,resI,aI,usum,icount,k)

      use precision_vars

      implicit none

      integer, parameter :: ivarlen=4,is=9

      real(wp), dimension(ivarlen),    intent(inout) :: uvec
      integer,                         intent(in   ) :: iprob,L
      real(wp),                        intent(in   ) :: ep,dt
      integer,                         intent(in   ) :: nveclen,iDT
      real(wp), dimension(ivarlen,is), intent(inout) :: resE,resI
      real(wp),                        intent(in   ) :: aI!from aI(L,L)
      real(wp), dimension(ivarlen),    intent(in   ) :: usum
      integer,                         intent(inout) :: icount,k
      

      integer                              :: ival,i,j  !Do loop variables
      integer                              :: programStep !input to problemsub
      real(wp), dimension(ivarlen)         :: uveciter !Storage of uvec
      real(wp), dimension(ivarlen)         :: uexact !input to problemsub (not needed)
      real(wp)                             :: tfinal !output of problemsub (not needed)
      real(wp), dimension(ivarlen)         :: Rnewton !Newton 
      real(wp), dimension(ivarlen,ivarlen) :: xjac,xjacinv !Jacobians
      real(wp)                             :: tmp !temp variable for accuracy
!============================================================
      
      Rnewton(:)=0.0_wp
      do k = 1,20

        icount = icount + 1

        uveciter(:) = uvec(:) !store old uvec

        programStep=2
        call problemsub(iprob,programStep,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE(1,L),resI(1,L))!,aI,xjac)
        !call RHS(uvec,resI(1,L),resE(1,L),dt,ep,iprob,sigma,rho,beta) !get res() for newton iteration

        Rnewton(:) =uvec(:)-aI*resI(:,L)-usum(:)
           
        !**GET INVERSE JACOBIAN**
        programStep=3
        call problemsub(iprob,programStep,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE(1,L),resI(1,L),aI,xjac)
!!        call JACOB2(nveclen,xjac,xjacinv)
!!      call Jacobian(uvec,xjac,dt,ep,aI,iprob,nveclen,sigma,rho,beta)
        call Invert_Jacobian(nveclen,xjac,xjacinv)
        !call JACOB(uvec,xjacinv,dt,ep,aI,iprob,nvecLen,sigma,rho,beta)

!        call Jacobian_CSR(uvec,dt,ep,aI,iprob,nveclen,sigma,rho,beta)
!        call ilu0(nveclen,aJac,jaJac,aluJac,jluJac,iw,ierr)
!        Backsolve to get solution

        do i = 1,nvecLen
          do j = 1,nvecLen
            uvec(i) = uvec(i) - xjacinv(i,j)*Rnewton(j)     !u^n+1=u^n-J^-1*F
          enddo
       enddo

       tmp = 0.0_wp
       do ival = 1,nvecLen
         tmp = tmp + abs(uvec(ival)-uveciter(ival)) !check accuracy of zeros
       enddo
       if(tmp.lt.1.0e-12_wp) exit

      enddo

      if (.false.) then
      !print*,'newton',resE!(1,L)
      write(*,*),uvec
      write(*,*),iprob
      write(*,*)L
      write(*,*)ep
      write(*,*)dt
      write(*,*)nveclen
      write(*,*)iDT
      write(*,*)resE
      write(*,*)resI
      write(*,*)aI
      write(*,*)usum
      write(*,*)'icount',icount
      write(*,*)ival
      write(*,*)i              
      write(*,*)j
      write(*,*)k
      write(*,*)programStep
      write(*,*)uveciter
      write(*,*)uexact
      write(*,*)tfinal
      write(*,*)Rnewton
      write(*,*)xjac
      write(*,*)xjacinv
      write(*,*)tmp
      endif
      return
      end subroutine 
