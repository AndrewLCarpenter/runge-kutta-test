      module Newton
            
      use precision_vars
      use QR_Module
      use control_variables
      use Runge_Kutta
      
      implicit none
      
      public :: Newton_Iteration
      private      
      
      contains
!==============================================================================      
      subroutine Newton_Iteration(uvec,iprob,L,ep,dt,nveclen,iDT,time,&
     & resE,resI,aI,usum,icount,k)


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
!////////////////////////////
!set typ to 1 for line search
!            2 for newton
!            3 for QR
!/////////////////////////////     
      typ = 2
!----------------LINE SEARCH---------------------------------------------------
      if (typ==1) then
        call newt_line_search(uvec,iprob,L,ep,dt,nveclen,time,resE,resI, &
     &                            aI,usum,icount,k)
!-------------------NEWTON AND QR----------------------------------------------
      else      

        do k = 1,20
          icount = icount + 1

          uveciter(:) = uvec(:) !store old uvec

          call Build_Rnewton(Rnewton,uvec,usum,ep,dt,time,aI,iprob,L)
 
          !**GET INVERSE JACOBIAN**
!!          call Build_Jac()
          programStep=3
          call problemsub(iprob,programStep,nveclen,uvec,ep,uexact,dt,&
     &                     tfinal,iDT,time,resE(1,L),resI(1,L),aI,xjac)
   
           if (typ==2) then !**NEWTON ITERATION**
             !print*,xjac

             xjacinv=Mat_invert(xjac)
!----------------------QR------------------------------------------------------             
           elseif (typ==3) then
             call QR_decomp(xjac,nveclen,rnewton,uvec)
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
         
         if(tmp.lt.1.0e-12_wp) then
           ierr = 0
           exit 
         endif
        enddo

      endif
!----------------END NEWTON AND QR---------------------------------------------      
!     if(ierr == 0) write(*,*)'Newton Iteration did not converge to machine precision'

      end subroutine Newton_Iteration
!==============================================================================
      
      subroutine Build_Rnewton(Rnewton,uvec,usum,ep,dt,time,aI,iprob,L)
            
      real(wp), dimension(:), intent(  out) :: Rnewton
      real(wp), dimension(:), intent(in   ) :: uvec,usum
      real(wp),               intent(in   ) :: ep,dt,time,aI
      integer,                intent(in   ) :: iprob,L
      
      real(wp), dimension(size(uvec),size(uvec)) :: xjac
      integer :: iDT, programStep,nveclen
      real(wp) :: tfinal      
      real(wp), dimension(size(uvec)) :: uexact
      real(wp), dimension(size(uvec),is) :: resE,resI
      
      programStep=2
      call problemsub(iprob,programStep,nveclen,uvec,ep,uexact,dt,&
     &                     tfinal,iDT,time,resE(1,L),resI(1,L),aI,xjac)
  
      Rnewton(:) = uvec(:)-aI*resI(:,L)-usum(:)
      
      end subroutine Build_Rnewton

!==============================================================================
!  PERFORMS LINE SEARCH AND NEWTON ITERATION 
      subroutine newt_line_search(uvec,iprob,L,ep,dt,nveclen,time,resE,resI, &
     &                            aI,usum,icount,k)
      
      real(wp), dimension(nveclen),    intent(inout) :: uvec
      integer,                         intent(in   ) :: iprob
      integer,                         intent(in   ) :: L
      real(wp),                        intent(in   ) :: ep,dt
      integer,                         intent(in   ) :: nveclen
      real(wp),                        intent(in   ) :: time
      real(wp), dimension(nveclen,is), intent(inout) :: resE,resI
      real(wp),                        intent(in   ) :: aI!from aI(L,L)
      real(wp), dimension(nveclen),    intent(in   ) :: usum
      integer,                         intent(inout) :: icount,k
      
      integer                              :: j  !Do loop variables
      integer                              :: programStep !input to problemsub
      integer                              :: ierr,iDT
      real(wp), dimension(nveclen)         :: uveciter !Storage of uvec
      real(wp), dimension(nveclen)         :: uexact !input to problemsub (not needed)
      real(wp)                             :: tfinal !output of problemsub (not needed)
      real(wp), dimension(nveclen)         :: Rnewton !Newton 
      real(wp), dimension(nveclen,nveclen) :: xjac,xjacinv !Jacobian

      real(wp)                             :: tmp !temp variable for accuracy
      real(wp) :: al,rnorm,rnormt
      real(wp), dimension(nveclen) :: dxi
      
      do k = 1,150
        icount = icount + 1

        programStep=2
        call problemsub(iprob,programStep,nveclen,&
      &          uvec,ep,uexact,dt,tfinal,iDT,time,resE(1,L),resI(1,L),aI,xjac)

        Rnewton(:) = uvec(:)-aI*resI(:,L)-usum(:)
        rnorm = sqrt(dot_product(Rnewton(:),Rnewton(:)))

        programStep=3
        call problemsub(iprob,programStep,nveclen,&
      &          uvec,ep,uexact,dt,tfinal,iDT,time,resE(1,L),resI(1,L),aI,xjac)

        xjacinv=Mat_invert(xjac)

        dxi = MatMul(xjacinv,Rnewton(:))

        al = 1.0_wp
        do j = 1,10    !   under-relax the value of the parameter alpha
          uveciter(:) = uvec(:) - al*dxi

          programStep=2
          call problemsub(iprob,programStep,nveclen,&
     &       uveciter,ep,uexact,dt,tfinal,iDT,time,resE(1,L),resI(1,L),aI,xjac)

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
      end subroutine newt_line_search
      
!==============================================================================
!  PEFORMS QR DECOMPOSITION USING QR       
      subroutine QR_decomp(xjac,nveclen,rnewton,uvec)
      
      real(wp), dimension(:,:), intent(inout) :: xjac
      integer,                  intent(in   ) :: nveclen
      real(wp), dimension(:),   intent(inout) :: rnewton,uvec
      
      real(wp), dimension(nveclen) :: wrk_c,wrk_d
      logical                :: snglr
      
      call qrdcmp(xjac,nveclen,nveclen,wrk_c,wrk_d,snglr)
      if (snglr) then
        write(*,*)'matrix is singular in Newton Iteration: Stopping'
        stop
      else
        call qrsolv(xjac,nveclen,nveclen,wrk_c,wrk_d,rnewton)
      endif
      uvec(:)=uvec(:)-rnewton(:) 
      end subroutine QR_decomp       
      
!==============================================================================
!  INVERTS MATRIX OF SIZE 2X2 TO 4X4
      function Mat_invert(mat)

      real(wp), dimension(:,:), intent(in   ) :: mat
      integer :: mat_size
      !real(wp), dimension(:,:), allocatable, intent(  out) :: xinv
      real(wp), dimension(size(mat(:,1)),size(mat(1,:))) :: xinv
      real(wp), dimension(size(mat(:,1)),size(mat(1,:))) :: Mat_invert
      real(wp) :: x11,x12,x13,x14
      real(wp) :: x21,x22,x23,x24
      real(wp) :: x31,x32,x33,x34
      real(wp) :: x41,x42,x43,x44
      real(wp) :: det,detI   

      mat_size=size(mat(:,1))
     
      if(mat_size==2)then
        det = (mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1))
        xinv(1,1) =  mat(2,2)/det
        xinv(1,2) = -mat(1,2)/det
        xinv(2,1) = -mat(2,1)/det
        xinv(2,2) =  mat(1,1)/det
      elseif(mat_size==3)then

        x11 = mat(1,1)
        x12 = mat(1,2)
        x13 = mat(1,3)
        x21 = mat(2,1)
        x22 = mat(2,2)
        x23 = mat(2,3)
        x31 = mat(3,1)
        x32 = mat(3,2)
        x33 = mat(3,3)

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

        x11 = mat(1,1)
        x12 = mat(1,2)
        x13 = mat(1,3)
        x14 = mat(1,4)
        x21 = mat(2,1)
        x22 = mat(2,2)
        x23 = mat(2,3)
        x24 = mat(2,4)
        x31 = mat(3,1)
        x32 = mat(3,2)
        x33 = mat(3,3)
        x34 = mat(3,4)
        x41 = mat(4,1)
        x42 = mat(4,2)
        x43 = mat(4,3)
        x44 = mat(4,4)

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
      Mat_invert=xinv
      return
      end function Mat_invert
!==============================================================================      
      end module Newton
