      subroutine Newton_Iteration(uvec,iprob,L,ep,dt,nveclen,iDT,time,&
     & resE,resI,aI,usum,icount,k)

      use precision_vars
      use QR_Module
      use control_variables
     ! use line_search_mod

      implicit none

      integer, parameter :: is=9

      real(wp), dimension(nveclen),    intent(inout) :: uvec
      integer,                         intent(in   ) :: iprob
      integer,                         intent(in   ) :: L
      real(wp),                        intent(in   ) :: ep,dt
      integer,                         intent(in   ) :: nveclen,iDT
      real(wp), dimension(nveclen,is), intent(inout) :: resE,resI
      real(wp),                        intent(in   ) :: aI!from aI(L,L)
      real(wp), dimension(nveclen),    intent(in   ) :: usum
      integer,                         intent(inout) :: icount,k
      real(wp),                        intent(in   ) :: time
      
      integer                              :: i,j  !Do loop variables
      integer                              :: programStep !input to problemsub
      integer                              :: ierr,typ
      real(wp), dimension(nveclen)         :: uveciter !Storage of uvec
      real(wp), dimension(nveclen)         :: uexact !input to problemsub (not needed)
      real(wp)                             :: tfinal !output of problemsub (not needed)
      real(wp), dimension(nveclen)         :: Rnewton !Newton 
      real(wp), dimension(nveclen,nveclen) :: xjac,xjacinv !Jacobians
      real(wp)                             :: tmp !temp variable for accuracy
      logical :: snglr,check
      !line search variables
      real(wp) :: al,rnorm,rnormt
      real(wp), dimension(nveclen) :: dxi
      real(wp), dimension(nveclen) :: wrk_c,wrk_d
      real, dimension(nveclen) :: uvectemp
      
      ierr = 0 
!============================================================
!set typ to 1 for line search
!            2 for newton
!            3 for QR

      typ = 2
!----------------LINE SEARCH----------------------------------------------------
      if (typ==1) then !true for linesearch
        !uvectemp=sngl(uvec)
        !call newt(uvectemp,nveclen,check)  
        !uvec=dble(uvectemp)
        do k = 1,150
          icount = icount + 1
  
          programStep=2
          call problemsub(iprob,programStep,nveclen,&
     &                 uvec,ep,uexact,dt,tfinal,iDT,time,resE(1,L),resI(1,L),aI,xjac)
     
          Rnewton(:) = uvec(:)-aI*resI(:,L)-usum(:)
          rnorm = sqrt(dot_product(Rnewton(:),Rnewton(:)))
  
          programStep=3
          call problemsub(iprob,programStep,nveclen,&
     &                 uvec,ep,uexact,dt,tfinal,iDT,time,resE(1,L),resI(1,L),aI,xjac)

          call Invert_Jacobian(nveclen,xjac,xjacinv)
  
          dxi = MatMul(xjacinv,Rnewton(:))
  
          al = 1.0_wp
          do j = 1,10    !   under-relax the value of the parameter alpha
  
            uveciter(:) = uvec(:) - al*dxi
  
            programStep=2
            call problemsub(iprob,programStep,nveclen,&
     &                 uveciter,ep,uexact,dt,tfinal,iDT,time,resE(1,L),resI(1,L),aI,xjac)

            Rnewton(:) =uveciter(:)-aI*resI(:,L)-usum(:)
            rnormt = sqrt(dot_product(Rnewton(:),Rnewton(:)))
  
             if((rnormt >= rnorm) .and. (rnorm >= 1.0e-4_wp)) then
               al = 0.5_wp * al
             else
               exit
             endif
          enddo

          uvec(:) = uvec(:) - al*dxi

          rnorm = rnormt
         if (k >= 125) print*,'L',L,'k',k,'tmp',rnorm,'j',j
          if(rnorm <= 1.0e-9_wp) then
            ierr = 0
            return
          endif
        enddo
!-------------------NEWTON AND QR----------------------------------------------
      else      

        do k = 1,20
          icount = icount + 1

          uveciter(:) = uvec(:) !store old uvec

          programStep=2
          call problemsub(iprob,programStep,nveclen,uvec,ep,uexact,dt,&
     &                     tfinal,iDT,time,resE(1,L),resI(1,L),aI,xjac)
  
          Rnewton(:) =uvec(:)-aI*resI(:,L)-usum(:)

          !**GET INVERSE JACOBIAN**
          programStep=3
          call problemsub(iprob,programStep,nveclen,uvec,ep,uexact,dt,&
     &                     tfinal,iDT,time,resE(1,L),resI(1,L),aI,xjac)
   
           if (typ==2) then !**NEWTON ITERATION**
             call Invert_Jacobian(nveclen,xjac,xjacinv)
!----------------------QR------------------------------------------------------             
           elseif (typ==3) then
             call qrdcmp(xjac,nveclen,nveclen,wrk_c,wrk_d,snglr)
         
             if(.not.snglr) then
               call qrsolv(xjac,nveclen,nveclen,wrk_c,wrk_d,rnewton)
             else
               write(*,*)'matrix is singular in Newton Iteration: Stopping'
               stop
             endif
             uvec(:)=uvec(:)-rnewton(:)  
            endif
!---------------------END QR---------------------------------------------------            
!          Backsolve to get solution

          if (typ==2) then !**NEWTON ITERATION**
          !**want to use this but there is a different truncation/round off error compared to the explicit do loop
          ! uvec(:)=uvec(:)-matmul(xjacinv,Rnewton)
           do i = 1,nvecLen
             do j = 1,nvecLen
               uvec(i) = uvec(i) - xjacinv(i,j)*Rnewton(j)     !u^n+1=u^n-J^-1*F
             enddo
           enddo
         endif !**NEWTON ITERATION**
         
         tmp = sum(abs(uvec(:)-uveciter(:))) !check accuracy of zeros
         !tmp=sqrt(dot_product(rnewton,rnewton))
         if (k>=15) print*,'tmp',tmp,'L',L,'k',k,'time',time
         !if (tmp/=tmp)  print*, 'hit nan', uveciter
         
         if(tmp.lt.1.0e-11_wp) then
           ierr = 0
           return 
         endif
        enddo
      endif
!----------------END NEWTON AND QR---------------------------------------------      
!     if(ierr == 0) write(*,*)'Newton Iteration did not converge to machine precision'

      end subroutine Newton_Iteration
