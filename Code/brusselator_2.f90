!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate Brusselator
!
! x'=a+x*x*y-b*x-x
! y'=b*x-x*x*y
!
! (du/dt)=(df/dx df/dy)*(u)
! (dv/dt) (dg/dx dg/dy) (v)
!
! J=(b-1 a)
!   (-b -a)
!
! u'=1+u*u*v+ep*d^2u/dx^2
! v'=3*u-u*u*v+ep*d^2v/dx^2
! where:
! u(x,t), v(x,t)
! u(0,t)=u(1,t)=1       v(0,t)=v(1,t)=3
! u(x,0)=1+sin(2*pi*x)  v(x,0)=3
! 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
      module Brusselator_mod
      
      use precision_vars,    only: wp,two,pi
      
      implicit none; save            
      
      private
      public :: Brusselator
      
      integer,  parameter :: vecl=16     ! total uvec length == 2*vecl
      real(wp), parameter :: xL=0.0_wp, xR=1.0_wp
      real(wp), dimension(vecl) :: x
      real(wp) ::dx
      real(wp), dimension(4), parameter :: d1vec0= (/-24.0_wp/17.0_wp,  &
                                                  &  +59.0_wp/34.0_wp,  &
                                                  &   -4.0_wp/17.0_wp,  &
                                                  &   -3.0_wp/34.0_wp/)
                                                                      
      real(wp), dimension(4), parameter :: d1vec1= (/+ 3.0_wp/34.0_wp,  &
                                                  &  + 4.0_wp/17.0_wp,  &
                                                  &  -59.0_wp/34.0_wp,  &
                                                  &  +24.0_wp/17.0_wp/)
      
      contains
      
!==============================================================================      
      subroutine Brusselator(nveclen,ep,dt,tfinal,iDT,resE_vec,resI_vec,akk)

      use control_variables, only: temporal_splitting,probname,Jac_case, &
     &                             tol,dt_error_tol,uvec,uexact,programstep
      use SBP_Coef_Module,   only: Define_CSR_Operators,D2,nnz_D2
      use Jacobian_CSR_Mod,  only: Allocate_Jac_CSR_Storage
!-----------------------VARIABLES----------------------------------------------

      !INIT vars
      real(wp), intent(in   ) :: ep
      real(wp), intent(inout) :: dt
      integer,  intent(  out) :: nveclen
      real(wp), intent(  out) :: tfinal
      integer,  intent(in   ) :: iDT

!      real(wp), dimension(81,vecl+1) :: ExactTot
!      real(wp)                       :: diff
      integer                        :: i

      !RHS vars
      real(wp), dimension(vecl*2), intent(  out) :: resE_vec,resI_vec
      real(wp), dimension(vecl*2)                :: dudt_Dx,dudt_Source
            
      !Jacob vars
      real(wp), intent(in   ) :: akk
!------------------------------------------------------------------------------
      
      Program_Step_Select: select case(programstep)
        !**Pre-initialization. Get problem name and vector length**
        case('INITIALIZE_PROBLEM_INFORMATION')
          nvecLen = 2*vecl
          probname='Brusselat'     
          tol=9.9e-11_wp  
          dt_error_tol=1.0e-13_wp
          Jac_case='SPARSE'
          
          call grid()

        !**Initialization of problem information**        
        case('SET_INITIAL_CONDITIONS')
          
          !Allocate derivative operators
          if(.not. allocated(D2))call Define_CSR_Operators(vecl,dx) 
           
          !Time information
          dt = 0.25_wp/10**((iDT-1)/20.0_wp) ! timestep   
          tfinal = 1.0_wp                   ! final time
 
          uexact(:)=0.0_wp         
          !**Exact Solution** 
!          open(unit=39,file='exact.brusselator4.data')
!          rewind(39)
!          do i=1,81
!            read(39,*)ExactTot(i,1),ExactTot(i,2)
!            ExactTot(i,3) = 1.0_wp/10**((i-1)/(10.0_wp))  !  used for 81 values of ep
!          enddo
!          do i=1,81
!            diff = abs(ExactTot(i,3) - ep)
!            if(diff.le.1.0e-10_wp)then
!              uexact(:) = ExactTot(i,:vecl)
!              exit
!            endif
!          enddo

          !**IC**
          uvec(1:2*vecl:2) = 1.0_wp+sin(two*pi*x(:))
          uvec(2:2*vecl:2) = 3.0_wp
      
          !**Exact Solution**
          uexact(:)=0.0_wp
          
        case('BUILD_RHS')
          call Bruss_dudt(uvec,dudt_Dx,dudt_Source,ep)
          choose_RHS_type: select case (Temporal_Splitting)
            case('IMPLICIT')
              resE_vec(:)=0.0_wp
              resI_vec(:)=dt*(dudt_Dx(:)+dudt_Source(:))
            case('IMEX')
              resE_vec(:)=dt*dudt_Dx(:)
              resI_vec(:)=dt*dudt_Source(:)
          end select choose_RHS_type
       
        case('BUILD_JACOBIAN')
          choose_Jac_type: select case(temporal_splitting)
            case('IMEX') ! For IMEX schemes
!              xjac(1,1) = 1.0_wp-akk*dt*(2.0_wp*uvec(1)*uvec(2)*ep)
!              xjac(1,2) = 0.0_wp-akk*dt*(       uvec(1)*uvec(1)*ep)
!              xjac(2,1) = 0.0_wp-akk*dt*(2.0_wp*uvec(1)*uvec(2)*ep)
!              xjac(2,2) = 1.0_wp-akk*dt*(      -uvec(1)*uvec(1)*ep)
            case('IMPLICIT') ! For fully implicit schemes
              call Allocate_Jac_CSR_Storage(vecl*2,nnz_D2*2)
              call Build_Jac(uvec,ep,dt,akk)
!              xjac(1,1) = 1.0_wp-akk*dt*(   2.0_wp*uvec(1)*uvec(2)/ep-bb-1.0_wp)
!              xjac(1,2) = 0.0_wp-akk*dt*(          uvec(1)*uvec(1)/ep)
!              xjac(2,1) = 0.0_wp-akk*dt*(bb-2.0_wp*uvec(1)*uvec(2)/ep)
!              xjac(2,2) = 1.0_wp-akk*dt*(         -uvec(1)*uvec(1)/ep)
          end select choose_Jac_type
        
      end select Program_Step_select
      end subroutine Brusselator
!==============================================================================
!==============================================================================
!==============================================================================
! PRODUCES GRID
      subroutine grid()
      integer :: i
      
      do i=1,vecl
        x(i)= xL + (xR-xL)*(i-1.0_wp)/(vecl-1.0_wp)
      enddo

      dx=x(2)-x(1)      
      
      end subroutine grid
!==============================================================================
! DEFINES RHS 
      subroutine Bruss_dUdt(vec,dudt_d_comb,dudt_s_comb,eps)
      
      use SBP_Coef_Module, only: Pinv,D2,jD2,iD2 
      use matvec_module,   only: amux

      real(wp), dimension(vecl),   intent(in   ) :: vec
      real(wp), dimension(vecl*2), intent(  out) :: dudt_d_comb,dudt_s_comb
      real(wp),                    intent(in   ) :: eps

      real(wp), dimension(vecl) :: gsat_u,gsat_v,dudt_s,dvdt_s,d2udx2,d2vdx2,u,v
     ! real(wp), dimension(vecl)                :: f, df, dfv, gsat, wrk    
     ! real(wp)                                 :: uL,uR,du0,du1,a0,a1,g0,g1
      
      gsat_u = 0.0_wp ; gsat_v = 0.0_wp !; df  = 0.0_wp ; dfv = 0.0_wp ;

      u=vec(1:vecl*2:2); v=vec(2:vecl*2:2)
      
      !Boundary conditions---> buggy
!      gsat_u(1   )=(-1.0_wp)*Pinv(1   )*(u(1   )-1.0_wp)
!      gsat_u(vecl)=( 1.0_wp)*Pinv(vecl)*(u(vecl)-1.0_wp)
!      gsat_v(1   )=(-1.0_wp)*Pinv(1   )*(v(1   )-3.0_wp)
!      gsat_v(vecl)=( 1.0_wp)*Pinv(vecl)*(v(vecl)-3.0_wp)
      
      dudt_s(:)=1.0_wp+u(:)*v(:)*v(:)-4.0_wp*u(:)+gsat_u(:)
      dvdt_s(:)=      -u(:)*u(:)*v(:)+3.0_wp*u(:)+gsat_v(:)
            
      call amux(vecl,u,d2udx2,D2,jD2,iD2)
      call amux(vecl,v,d2vdx2,D2,jD2,iD2)          
              
      dudt_d_comb(1:vecl*2:2)=eps*d2udx2(:)
      dudt_d_comb(2:vecl*2:2)=eps*d2vdx2(:)
      dudt_s_comb(1:vecl*2:2)=dudt_s(:)
      dudt_s_comb(2:vecl*2:2)=dvdt_s(:)

!--------------Inflow Boundary Condition---------------------------------------

!      uR = u(1)

!      a0 = third*(uR + sqrt(uR**2+1.0e-28_wp))    !a0 = third*(uR + abs(uR))
!      g0 = a0 * NL_Burg_exactsolution(x(1),time,eps)                &
!         - eps* NL_Burg_exact_derivative(x(1),time,eps)

!      du0 = dot_product(d1vec0(1:4),u(1:4)) / dx

!      gsat(1) = gsat(1) + sig0 * Pinv(1) * (a0*uR - eps*du0 - g0) 

!---------------Outflow Boundary Condition-------------------------------------

!      uL = u(vecl)

!      a1 = third*(uL - sqrt(uL**2+1.0e-28_wp)) !a1 = third*(uL - abs(uL))
!      g1 = a1 * NL_Burg_exactsolution(x(vecl),time,eps)                &
!         - eps* NL_Burg_exact_derivative(x(vecl),time,eps)
!      du1 = dot_product(d1vec1(1:4),u(vecl-3:vecl)) / dx

!      gsat(vecl) = gsat(vecl) + sig1 * Pinv(vecl) * (a1*uL-eps*du1-g1)

!--------------------Sum all terms---------------------------------------------

!      dudt(:) = dt*(eps*dfv(:) - df(:) + gsat(:))

      end subroutine Bruss_dUdt
!==============================================================================
!  CREATES JACOBIAN
      subroutine Build_Jac(vec,eps,dt,akk)

      use SBP_Coef_Module,  only: Pinv,D1,D2,jD1,jD2,iD1,iD2,nnz_D2,nnz_D1       
      use matvec_module,    only: amux
      use unary_mod,        only: aplb,aplsca,amudia,diamua,apldia
      use Jacobian_CSR_Mod, only: iaJac,jaJac,aJac
      
      real(wp), dimension(vecl*2), intent(in) :: vec
      real(wp),                  intent(in) :: eps,dt,akk    
        
!     integer,  dimension(vecl+1) :: iwrk1,iwrk2,iwrk3,iJac
!      integer,  dimension(nnz_D2) :: jwrk1,jwrk2,jwrk3,jJac
!      real(wp), dimension(nnz_D2) :: wrk1,wrk2,wrk3,Jac,wrk4,eps_d2      
      
      integer,  dimension(vecl)  :: iw
!      real(wp), dimension(vecl)  :: diag
      integer,  dimension(3)     :: ierr
!      real(wp)                   :: uL,a1_d,uR,a0_d,a0,a1
      
      
      
      real(wp), dimension(vecl) :: D1x4,v,u,v_sqr,u_sqr
      integer :: counter, i
      integer, dimension(vecl+1)::iaJv,iaJu,iwrk1,iwrk2,iwrk3
      
      real(wp), dimension(nnz_D1) :: wrk3,wrk2
      integer, dimension(nnz_D1) :: jwrk3,jwrk2
      
      real(wp), dimension(nnz_D2) ::wrk1,aJu,aJv,eps_d2
      integer,  dimension(nnz_D2) ::jwrk1,jaJu,jaJv
      
      integer, dimension(vecl*2+1) ::iwrkJ
      integer, dimension(nnz_D2*2) ::jwrkJ
      real(wp), dimension(nnz_D2*2)::wrkJ
  
      u=vec(1:vecl*2:2); v=vec(2:vecl*2:2)
! Ju=v^2*D1-4*D1+eps*D2
! D1x4=-4*D1
! eps_d2=eps*D2
! wrk1=D1x4+eps_d2
! Jv=-u^2*D1    +eps*D2
     
      eps_d2=eps*D2
      D1x4=(-4.0_wp)*D1
      
      call aplb(vecl,vecl,1,eps_d2,jD2,iD2,D1x4,jD1,iD1,wrk1,jwrk1,iwrk1,nnz_D2,iw,ierr(1))

      v_sqr(:)=v(:)*v(:)
      
      call diamua(vecl,1,D1,jD1,iD1,v_sqr,wrk2,jwrk2,iwrk2)!nnz_D1
      
      call aplb(vecl,vecl,1,wrk2,jwrk2,iwrk2,wrk1,jwrk1,iwrk1,aJu,jaJu,iaJu,nnz_D2,iw,ierr(2))
      
      u_sqr=(-1.0_wp)*u(:)*u(:)
      
      call diamua(vecl,1,D1,jD1,iD1,u_sqr,wrk3,jwrk3,iwrk3)!nnz_D1
      
      call aplb(vecl,vecl,1,eps_d2,jD2,iD2,wrk3,jwrk3,iwrk3,aJv,jaJv,iaJv,nnz_D2,iw,ierr(3))
      
      !boundary conditions -- buggy
!      aJu(1   )=aJu(1       )+(-1.0_wp)*Pinv(1   )
!      aJu(nnz_D2)=aJu(nnz_D2)+( 1.0_wp)*Pinv(vecl)
!      aJv(1   )=aJv(1       )+(-1.0_wp)*Pinv(1   )
!      aJv(nnz_D2)=aJv(nnz_D2)+( 1.0_wp)*Pinv(vecl)
     
     
      iwrkJ(1)=1
      counter=1
      do i=2,2*vecl,2
        iwrkJ(i)=iwrkJ(i-1)+iaJu(counter+1)-iaJu(counter)
        iwrkJ(i+1)=iwrkJ(i)+iaJv(counter+1)-iaJv(counter)
        counter=counter+1
      enddo   
      
      counter=1
      do i=1,2*vecl-1,2
        wrkJ( iwrkJ(i  ):iwrkJ(i+1)-1)= aJu(iaJu(counter):iaJu(counter+1)-1)
        wrkJ( iwrkJ(i+1):iwrkJ(i+2)-1)= aJv(iaJv(counter):iaJv(counter+1)-1)
        jwrkJ(iwrkJ(i  ):iwrkJ(i+1)-1)=jaJu(iaJu(counter):iaJu(counter+1)-1)*2-1
        jwrkJ(iwrkJ(i+1):iwrkJ(i+2)-1)=jaJv(iaJv(counter):iaJv(counter+1)-1)*2
        counter=counter+1
      enddo
      
      wrkJ=-akk*dt*wrkJ                                                       
      call aplsca(vecl,wrkJ,jwrkJ,iwrkJ,1.0_wp,iw)                          

      iaJac=iwrkJ
      jaJac=jwrkJ
      aJac=wrkJ
     
!---------------------dgsat/dt-------------------------------------------------
!      uL=u(vecl)
!      a1 = third*(uL - sqrt(uL**2+1.0e-28_wp))
!      a1_d=third*(1-uL/sqrt(uL**2+1.0e-28_wp))
!   
!      uR=u(1)
!      a0 = third*(uR + sqrt(uR**2+1.0e-28_wp))
!      a0_d=third*(1+uR/sqrt(uR**2+1.0e-28_wp))

!******Steps to make xjac******
!       jac=eps*d2-2/3*d1*u-1/3*u*d1-1/3*diag(d1*u)+dgsat/du
!1      wrk1=-2/3*d1*u
!2      wrk2=-1/3*u*d1
!3      wrk3=wrk1+wrk2
!4      eps_d2=eps*d2
!5      jac=eps_d2+wrk3-1/3*diag(d1*u)+dgsat/du
!6      wrk4=-akk*dt*jac
!
!       xjac=I-akk*dt*jac=-akk*dt*jac+(1)*I
!7      xjac=wrk4+(1)*I
! 
!8      xjac-> wrk4, iJac, jJac
!******************************
! --------------Make xjac -----------------------------------------------------
!      call amudia(vecl,1,D1,jD1,iD1,u,wrk1,jwrk1,iwrk1)                   !1
!      wrk1(:)=-twothird*wrk1(:)                                              !1
      
!      call diamua(vecl,1,D1,jD1,iD1,u,wrk2,jwrk2,iwrk2)                   !2
!      wrk2(:)=-third*wrk2(:)                                                 !2
      
!      call aplb(vecl,vecl,1,wrk1,jwrk1,iwrk1,wrk2,jwrk2,iwrk2,wrk3,&   !3
!     &          jwrk3,iwrk3,nnz_D2,iw,ierr(1))                               !3

!      eps_d2=eps*D2(:)                                                       !4
!      call aplb(vecl,vecl,1,eps_d2,jD2,iD2,wrk3,jwrk3,iwrk3,Jac,&     !5a
!     &          jJac,iJac,nnz_D2,iw,ierr(2))                                !5a

!      call amux(vecl,-third*u,diag,D1,jD1,iD1) !First-derivative of u vec!5b
!      call apldia(vecl,0,Jac,jJac,iJac,diag,Jac,jJac,iJac,iw)            !5b

     ! R - 5c
!      Jac(1)=Jac(1)+sig0*Pinv(1)*(a0+a0_d*(uR- &
!     &       NL_Burg_exactsolution(x(1),time,eps))) 
!      Jac(1:4)=Jac(1:4)+sig0*Pinv(1)*(-eps)*d1vec0(:)/dx
      
      ! L - 5c
!      Jac(nnz_D2)=Jac(nnz_D2)+sig1*Pinv(vecl)*(a1+a1_d*(uL- &
!     &            NL_Burg_exactsolution(x(vecl),time,eps))) 
!      Jac(nnz_D2-3:nnz_D2)=Jac(nnz_D2-3:nnz_D2)+ &
!     &                     sig1*Pinv(vecl)*(-eps)*d1vec1(:)/dx
     
!      wrk4=-akk*dt*Jac                                                       !6
!      call aplsca(vecl,wrk4,jJac,iJac,1.0_wp,iw)                          !7

      !8
!      iaJac=iJac
!      jaJac=jJac
!      aJac=wrk4
!-----------------------------------------------------------------------------  
      if (sum(ierr)/=0) then ! Catch errors
        print*,'Error building Jacobian'
        stop
      endif
   
      end subroutine Build_Jac
      
      end module Brusselator_mod
