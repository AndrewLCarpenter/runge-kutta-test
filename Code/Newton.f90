!*****************************************************************************
! Module to perform Newton Iteration in one of three modes: QR decomp, line 
! search, or normal Newton Iteration. 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
! QR_MODULE.F90                 *CONTAINS QR ROUTINES
! CONTROL_VARIABLES.F90         *ONTAINS VARIABLES USED IN THE PROGRAM
! RUNGE_KUTTA.F90               *CONTAINS RK CONSTANTS
! JACOBIAN_CSR_MOD.F90          *CONTAINS CSR JACOBIAN VARIABLES
! ILUT_MODULE.F90               *PERFORMS LU FACTORIZATION AND SOLVING
!******************************************************************************

      module Newton
            
      use precision_vars, only:wp

      implicit none; save
      
      public :: Newton_Iteration
      private      
      
      contains
      
!==============================================================================      
!  PERFORMS NEWTON ITERATION
      subroutine Newton_Iteration(iprob,L,ep,dt,nveclen,time,aI,icount,k)
      
      use control_variables, only: temporal_splitting,jac_case,uvec,usum,xjac
      use Jacobian_CSR_Mod,  only: Jacobian_CSR

      integer, parameter :: iter_max=20
      logical, parameter :: Line_search=.false. !set to TRUE for line search
      logical, parameter :: QR=.false. !set to TRUE for QR factorization

      integer,  intent(in   ) :: iprob,L
      real(wp), intent(in   ) :: ep,dt
      integer,  intent(in   ) :: nveclen
      real(wp), intent(in   ) :: time,aI
      integer,  intent(inout) :: icount
      integer,  intent(  out) :: k
      
      integer                              :: i,j  !Do loop variables
      real(wp), dimension(nveclen)         :: Rnewton,uveciter!Newton 
      real(wp), dimension(nveclen,nveclen) :: xjacinv !Jacobian

!------------------------------------------------------------------------------
      select case(temporal_splitting)
        case('EXPLICIT')
          do k = 1,iter_max
            icount = icount + 1
            uveciter(:) = uvec(:)
            uvec(:)=usum(:)
            if(check_exit(uveciter,k)) exit            
          enddo   
        case default !IMEX or IMPLICIT
          select case(Jac_case)
            case('SPARSE')
              do k = 1,iter_max
                icount = icount + 1
                uveciter(:) = uvec(:) !store old uvec

                Rnewton=Build_Rnewton(ep,dt,time,aI,iprob,L)
                call Build_Jac(ep,dt,time,aI,iprob,L)
              
                call LU_solver(Rnewton)
                
                if(check_exit(uveciter,k)) exit
                            
              enddo  
            case('DENSE')
              if (Line_search) then
                call newt_line_search(iprob,L,ep,dt,nveclen,time,aI,icount,k)    
              else
              
                do k = 1,iter_max
                  icount = icount + 1
                  uveciter(:) = uvec(:) !store old uvec   
                              
                  Rnewton=Build_Rnewton(ep,dt,time,aI,iprob,L)
                  call Build_Jac(ep,dt,time,aI,iprob,L)
                  
                  if (nveclen>4 .and. .not. QR) then !No explicit inverse
                    call Jacobian_CSR(nveclen,xjac) !convert dense xjac to csr
                    call LU_solver(Rnewton)
                    
                  elseif (nveclen>4 .and. QR) then
                    call QR_decomp(xjac,nveclen,rnewton) 
                     
                  elseif (nveclen<=4) then !Explicit inverse
                    xjacinv=Mat_invert(xjac)                  
                    do i = 1,nvecLen
                      do j = 1,nvecLen
                        uvec(i) = uvec(i) - xjacinv(i,j)*Rnewton(j)     !u^n+1=u^n-J^-1*F
                      enddo
                    enddo  
                  endif
                  
                  if(check_exit(uveciter,k)) exit
                  
                enddo
                
              endif                
          end select    
      end select
      end subroutine Newton_Iteration
      
!==============================================================================
! EXITS NEWTON ITERATION
      function check_exit(uveciter,k)
      
      use control_variables, only: tol,uvec
      
      real(wp), dimension(:), intent(in) :: uveciter
      integer,                intent(in) :: k
      real(wp)                           :: tmp
      logical                            :: check_exit
      
      check_exit=.false.
      tmp = sum(abs(uvec(:)-uveciter(:))) !check accuracy of zeros         
      if (k>=15) write(*,*),'tmp',tmp,'k',k
      if (tmp/=tmp) then
        print*,'stopping NaN k=',k
        stop
      endif      
      if(tmp < tol) check_exit=.true. 
      
      end function check_exit
!==============================================================================
!  CREATES RNEWTON FOR NEWTON ITERATION      
      function Build_Rnewton(ep,dt,time,aI,iprob,L)
      
      use control_variables, only: uvec,usum,resI
      use problemsub_mod,    only: problemsub
            
      real(wp), dimension(size(uvec))       :: Build_Rnewton
      real(wp),               intent(in   ) :: ep,dt,time,aI
      integer,                intent(in   ) :: iprob,L

      integer  :: iDT, programStep,nveclen
      real(wp) :: tfinal,dt_in  
      
      dt_in=dt
      
      programStep=2
      call problemsub(iprob,programStep,nveclen,ep,dt_in,tfinal,iDT,time,aI,L)
      Build_Rnewton(:) = uvec(:)-aI*resI(:,L)-usum(:)

      end function Build_Rnewton
!==============================================================================
!  CREATES JACOBIAN FOR NEWTON ITERATION      
      subroutine Build_Jac(ep,dt,time,aI,iprob,L)
      
      use problemsub_mod,    only: problemsub     
       
      real(wp), intent(in) :: ep,dt,time,aI
      integer,  intent(in) :: iprob,L
    
      integer  :: nveclen
      real(wp) :: tfinal,dt_in
      integer  :: programstep,iDT
      
      dt_in=dt
  
      programStep=3
      call problemsub(iprob,programStep,nveclen,ep,dt_in,tfinal,iDT,time,aI,L)
      
      end subroutine Build_Jac
      
!==============================================================================
! PERFORMS LU DECOMPOSITION AND SOLVING
      subroutine LU_solver(Rnewton)
      
      use control_variables, only: uvec      
      use Jacobian_CSR_Mod,  only: iaJac,jaJac,aJac,aLUJac,jLUJac,jUJac
      use ilut_module,       only: lusol,ilutp
      
      real(wp), dimension(:)               :: Rnewton
      integer                              :: nveclen,ierr=0
      real(wp), dimension(size(Rnewton))   :: r_wrk
      real(wp), dimension(size(Rnewton))   :: w
      integer,  dimension(2*size(Rnewton)) :: jw,iperm
      
      integer :: i,k
    
      nveclen=size(Rnewton)   

      call ilutp(nveclen,aJac,jaJac,iaJac,nveclen,1e-13_wp,0.1_wp,nveclen, &
     &           aLUJac,jLUJac,jUJac,size(alujac),w,jw,iperm,ierr)

      call lusol(nveclen,Rnewton,r_wrk,aLUJac,jLUJac,jUJac)
      uvec(:)=uvec(:)-r_wrk(:)        

      if (ierr/=0) then 
        print*,'Error in LU_solver!'
        print*,'ierr=',ierr
        stop
      endif
      end subroutine LU_solver 

!==============================================================================
!  PERFORMS LINE SEARCH AND NEWTON ITERATION 
      subroutine newt_line_search(iprob,L,ep,dt,nveclen,time,aI,icount,k)
      
      use control_variables, only: uvec,xjac
      
      integer,  intent(in   ) :: iprob,L
      real(wp), intent(in   ) :: ep,dt
      integer,  intent(in   ) :: nveclen
      real(wp), intent(in   ) :: time,aI
      integer,  intent(inout) :: icount,k
      
      integer                              :: j  !Do loop variables
      integer                              :: programStep !input to problemsub
      integer                              :: ierr,iDT
      real(wp)                             :: tfinal !output of problemsub (not needed)
      real(wp), dimension(nveclen)         :: Rnewton,ustor !Newton 
      real(wp), dimension(nveclen,nveclen) :: xjacinv !Jacobianxjac,

      real(wp) :: tmp !temp variable for accuracy
      real(wp) :: al,rnorm,rnormt
      real(wp), dimension(nveclen) :: dxi
      
      do k = 1,150
        icount = icount + 1
   
        Rnewton=Build_Rnewton(ep,dt,time,aI,iprob,L)
    
        rnorm = sqrt(dot_product(Rnewton(:),Rnewton(:)))

        call Build_Jac(ep,dt,time,aI,iprob,L)

        xjacinv=Mat_invert(xjac)

        dxi = MatMul(xjacinv,Rnewton(:))

        al = 1.0_wp
        do j = 1,10    !   under-relax the value of the parameter alpha
        
          ustor(:)=uvec(:)!becuause uvec is global, it needs to be temp
                          !stored so that calculations are done correctly
          uvec(:) = uvec(:) - al*dxi
             
          Rnewton=Build_Rnewton(ep,dt,time,aI,iprob,L)
    
          rnormt = sqrt(dot_product(Rnewton(:),Rnewton(:)))

          uvec(:)=ustor(:)!reset uvec from temp storage

           if((rnormt >= rnorm) .and. (rnorm >= 1.0e-4_wp)) then
             al = 0.5_wp * al
           else
             exit
           endif
        enddo

        uvec(:) = uvec(:) - al*dxi

        rnorm = rnormt
       if (k >= 140) print*,'L',L,'k',k,'tmp',rnorm,'j',j
        if(rnorm <= 1.0e-9_wp) then
          ierr = 0
          return
        endif
      enddo      
      end subroutine newt_line_search
      
!==============================================================================
!  PEFORMS QR DECOMPOSITION USING QR MODULE       
      subroutine QR_decomp(mat,nveclen,Rnewton)
      
      use control_variables, only: uvec
      use QR_Module, only:qrdcmp,qrsolv
      
      real(wp), dimension(:,:), intent(inout) :: mat
      integer,                  intent(in   ) :: nveclen
      real(wp), dimension(:),   intent(inout) :: Rnewton
      
      real(wp), dimension(nveclen) :: wrk_c,wrk_d
      logical                      :: snglr
      
      call qrdcmp(mat,nveclen,nveclen,wrk_c,wrk_d,snglr)
      if (snglr) then
        write(*,*)'matrix is singular in Newton Iteration: Stopping'
        stop
      else
        call qrsolv(mat,nveclen,nveclen,wrk_c,wrk_d,Rnewton)
      endif
      uvec(:)=uvec(:)-Rnewton(:) 
      end subroutine QR_decomp       
      
!==============================================================================
!  INVERTS MATRIX OF SIZE 2X2 TO 4X4
      function Mat_invert(mat)

      real(wp), dimension(:,:), intent(in   ) :: mat
      integer  :: mat_size
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

        detI = 1.0_wp/det

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

        detI = 1.0_wp/det

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
