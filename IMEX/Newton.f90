      module Newton
            
      use precision_vars
      use QR_Module
      use control_variables
      use Runge_Kutta
      
      public :: Newton_Iteration
      private      
      
      contains
!==============================================================================      
      subroutine Newton_Iteration(uvec,iprob,L,ep,dt,nveclen,iDT,time,&
     & resE,resI,aI,usum,icount,k)

      implicit none

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
      real(wp), dimension(nveclen,nveclen) :: xjac,xjacinv !Jacobian

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
!----------------LINE SEARCH---------------------------------------------------
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

          xjacinv=Invert_Jacobian(xjac)
  
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
             !print*,xjac

             xjacinv=Invert_Jacobian(xjac)
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
           exit 
         endif
        enddo

      endif
!----------------END NEWTON AND QR---------------------------------------------      
!     if(ierr == 0) write(*,*)'Newton Iteration did not converge to machine precision'

      end subroutine Newton_Iteration
!==============================================================================

      function Invert_Jacobian(xjac)

      implicit none
      
      real(wp), dimension(:,:), intent(in   ) :: xjac
      integer :: mat_size
      !real(wp), dimension(:,:), allocatable, intent(  out) :: xinv
      real(wp), dimension(size(xjac(:,1)),size(xjac(1,:))) :: xinv
      real(wp), dimension(size(xjac(:,1)),size(xjac(1,:))) :: Invert_Jacobian
      real(wp) :: x11,x12,x13,x14
      real(wp) :: x21,x22,x23,x24
      real(wp) :: x31,x32,x33,x34
      real(wp) :: x41,x42,x43,x44
      real(wp) :: det,detI   

      mat_size=size(xjac(:,1))
     
      if(mat_size==2)then
        det = (xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1))
        xinv(1,1) =  xjac(2,2)/det
        xinv(1,2) = -xjac(1,2)/det
        xinv(2,1) = -xjac(2,1)/det
        xinv(2,2) =  xjac(1,1)/det
      elseif(mat_size==3)then

        x11 = xjac(1,1)
        x12 = xjac(1,2)
        x13 = xjac(1,3)
        x21 = xjac(2,1)
        x22 = xjac(2,2)
        x23 = xjac(2,3)
        x31 = xjac(3,1)
        x32 = xjac(3,2)
        x33 = xjac(3,3)

        det = - x13*x22*x31 + x12*x23*x31 +  x13*x21*x32& 
     &        - x11*x23*x32 - x12*x21*x33 +  x11*x22*x33

        detI = 1./det

        xinv(1,1) =  (- x23*x32 + x22*x33) * detI
        xinv(1,2) =  (+ x13*x32 - x12*x33) * detI
        xinv(1,3) =  (- x13*x22 + x12*x23) * detI
        xinv(2,1) =  (+ x23*x31 - x21*x33) * detI
        xinv(2,2) =  (- x13*x31 + x11*x33) * detI
        xinv(2,3) =  (+ x13*x21 - x11*x23) * detI
        xinv(3,1) =  (- x22*x31 + x21*x32) * detI
        xinv(3,2) =  (+ x12*x31 - x11*x32) * detI
        xinv(3,3) =  (- x12*x21 + x11*x22) * detI

      elseif(mat_size==4)then

        x11 = xjac(1,1)
        x12 = xjac(1,2)
        x13 = xjac(1,3)
        x14 = xjac(1,4)
        x21 = xjac(2,1)
        x22 = xjac(2,2)
        x23 = xjac(2,3)
        x24 = xjac(2,4)
        x31 = xjac(3,1)
        x32 = xjac(3,2)
        x33 = xjac(3,3)
        x34 = xjac(3,4)
        x41 = xjac(4,1)
        x42 = xjac(4,2)
        x43 = xjac(4,3)
        x44 = xjac(4,4)

        det =& 
     &  (x14*x23*x32*x41 - x13*x24*x32*x41 - &
     &   x14*x22*x33*x41 + x12*x24*x33*x41 + &
     &   x13*x22*x34*x41 - x12*x23*x34*x41 - &
     &   x14*x23*x31*x42 + x13*x24*x31*x42 + &
     &   x14*x21*x33*x42 - x11*x24*x33*x42 - &
     &   x13*x21*x34*x42 + x11*x23*x34*x42 + &
     &   x14*x22*x31*x43 - x12*x24*x31*x43 - &
     &   x14*x21*x32*x43 + x11*x24*x32*x43 + &
     &   x12*x21*x34*x43 - x11*x22*x34*x43 - &
     &   x13*x22*x31*x44 + x12*x23*x31*x44 + &
     &   x13*x21*x32*x44 - x11*x23*x32*x44 - &
     &   x12*x21*x33*x44 + x11*x22*x33*x44)

        detI = 1./det

        xinv(1,1) = (&
     & -(x24*x33*x42) + x23*x34*x42 + x24*x32*x43 - &
     &   x22*x34*x43  - x23*x32*x44 + x22*x33*x44  ) * detI
        xinv(1,2) = (&
     &   x14*x33*x42  - x13*x34*x42 - x14*x32*x43 + &
     &   x12*x34*x43  + x13*x32*x44 - x12*x33*x44  ) * detI
        xinv(1,3) = (&
     & -(x14*x23*x42) + x13*x24*x42 + x14*x22*x43 - &
     &   x12*x24*x43  - x13*x22*x44 + x12*x23*x44  ) * detI
        xinv(1,4) = (&
     &   x14*x23*x32  - x13*x24*x32 - x14*x22*x33 + &
     &   x12*x24*x33  + x13*x22*x34 - x12*x23*x34  ) * detI
        xinv(2,1) = (&
     &   x24*x33*x41  - x23*x34*x41 - x24*x31*x43 + &
     &   x21*x34*x43  + x23*x31*x44 - x21*x33*x44  ) * detI
        xinv(2,2) = (&
     & -(x14*x33*x41) + x13*x34*x41 + x14*x31*x43 - &
     &   x11*x34*x43  - x13*x31*x44 + x11*x33*x44  ) * detI
        xinv(2,3) = (&
     &   x14*x23*x41  - x13*x24*x41 - x14*x21*x43 + &
     &   x11*x24*x43  + x13*x21*x44 - x11*x23*x44  ) * detI
        xinv(2,4) = (&
     & -(x14*x23*x31) + x13*x24*x31 + x14*x21*x33 - &
     &   x11*x24*x33  - x13*x21*x34 + x11*x23*x34  ) * detI
        xinv(3,1) = (&
     & -(x24*x32*x41) + x22*x34*x41 + x24*x31*x42 - &
     &   x21*x34*x42  - x22*x31*x44 + x21*x32*x44  ) * detI
        xinv(3,2) = (&
     &   x14*x32*x41  - x12*x34*x41 - x14*x31*x42 + &
     &   x11*x34*x42  + x12*x31*x44 - x11*x32*x44  ) * detI
        xinv(3,3) = (&
     & -(x14*x22*x41) + x12*x24*x41 + x14*x21*x42 - &
     &   x11*x24*x42  - x12*x21*x44 + x11*x22*x44  ) * detI
        xinv(3,4) = (&
     &   x14*x22*x31  - x12*x24*x31 - x14*x21*x32 + &
     &   x11*x24*x32  + x12*x21*x34 - x11*x22*x34  ) * detI
        xinv(4,1) = (&
     &   x23*x32*x41  - x22*x33*x41 - x23*x31*x42 + &
     &   x21*x33*x42  + x22*x31*x43 - x21*x32*x43  ) * detI
        xinv(4,2) = (&
     & -(x13*x32*x41) + x12*x33*x41 + x13*x31*x42 - &
     &   x11*x33*x42  - x12*x31*x43 + x11*x32*x43  ) * detI
        xinv(4,3) = (&
     &   x13*x22*x41  - x12*x23*x41 - x13*x21*x42 + &
     &   x11*x23*x42  + x12*x21*x43 - x11*x22*x43  ) * detI
        xinv(4,4) = (&
     & -(x13*x22*x31) + x12*x23*x31 + x13*x21*x32 - &
     &   x11*x23*x32  - x12*x21*x33 + x11*x22*x33  ) * detI
      endif
      Invert_Jacobian=xinv
      return
      end function Invert_Jacobian
      end module Newton
