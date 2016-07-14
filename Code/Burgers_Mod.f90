!*****************************************************************************
! Module containing routines to describe Burger's equation 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
! SBP_COEF_MODULE.F90           *DEFINES CSR OPERATORS 
! UNARY_MOD.F90                 *PERFORMS SPARSE MATRIX OPERATIONS
! MATVEC_MODULE.F90             *PERFORMS SPARSE MATRIX*VECTOR OPERATIONS
!******************************************************************************

      module Burgers_Module 

      use precision_vars, only: third, half, twothird, wp

      implicit none; save
  
      private
      public    ::  Burgers

!--------------------------------VARIABLES-------------------------------------         
      real(wp)  :: dx
      character(120),         parameter :: exact_solution = 'tanh'
!      character(120),         parameter :: exact_solution = 'Olver'
      real(wp),               parameter :: sig0 = -1.0_wp
      real(wp),               parameter :: sig1 = +1.0_wp
      
      integer,  parameter    :: vecl=16     
      real(wp), dimension(vecl) :: x      
      real(wp), parameter :: xL=0.0_wp,xR=1.0_wp

      real(wp), dimension(4), parameter :: d1vec0= (/-24.0_wp/17.0_wp,  &
                                                  &  +59.0_wp/34.0_wp,  &
                                                  &   -4.0_wp/17.0_wp,  &
                                                  &   -3.0_wp/34.0_wp/)
                                                                      
      real(wp), dimension(4), parameter :: d1vec1= (/+ 3.0_wp/34.0_wp,  &
                                                  &  + 4.0_wp/17.0_wp,  &
                                                  &  -59.0_wp/34.0_wp,  &
                                                  &  +24.0_wp/17.0_wp/)
!------------------------------------------------------------------------------      
      contains    
!==============================================================================
!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Burgers problem 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
! BURGERS_MOD.F90           *CONTAINS ROUTINES TO BUILD BURGER
! JACOBIAN_CSR_MOD.F90      *CONTAINS CSR JACOBIAN VARIABLES
!******************************************************************************

      subroutine Burgers(programStep,nveclen,ep,dt, &
     &                   tfinal,iDT,time,resE_vec,resI_vec,akk)

      use control_variables, only: temporal_splitting,probname,Jac_case, &
     &                             tol,dt_error_tol,uvec,uexact
      use SBP_Coef_Module,   only: Define_CSR_Operators,nnz_D2,D1
      use Jacobian_CSR_Mod,  only: Allocate_CSR_Storage

!-----------------------VARIABLES----------------------------------------------
      integer, intent(in   ) :: programStep

      !INIT vars     
      real(wp), intent(in   ) :: ep
      real(wp), intent(inout) :: dt
      integer,  intent(  out) :: nveclen
      real(wp), intent(  out) :: tfinal
      real(wp), intent(in   ) :: time
      integer,  intent(in   ) :: iDT

      real(wp)                :: tinitial
      
      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE_vec,resI_vec
      
      !Jacob vars
      real(wp), intent(in   ) :: akk
!------------------------------------------------------------------------------

      !**Pre-initialization. Get problem name and vector length**
      if (programStep==-1) then
        nveclen = vecl
        probname='Burgers  '     
        tol=1.0e-12_wp  
        dt_error_tol=5.0e-14_wp
        jac_case='SPARSE'
        call grid()           
      !**Initialization of problem information**        
      elseif (programStep==0) then

        if(.not. allocated(D1)) call Define_CSR_Operators(vecl,dx)   

        tinitial=0.0_wp  ! initial time
        tfinal = 0.5_wp  ! final time           
!        dt = 0.00005_wp*0.1_wp/10**((iDT-1)/20.0_wp) ! timestep  explicit 
        dt = 0.5_wp*0.1_wp/10**((iDT-1)/20.0_wp) ! timestep  implcit
        call exact_Burg(uvec,ep,tinitial) !set initial conditions   
        call exact_Burg(uexact,ep,tfinal) !set exact solution at tfinal
              
      !**RHS and Jacobian**
      elseif (programStep>=1) then

        select case (Temporal_Splitting)
        
          case('EXPLICIT')
            !**RHS**
            if (programStep==1 .or.programStep==2) then
              call Burgers_dUdt(uvec,resE_vec,time,ep,dt)
              resI_vec(:)=0.0_wp
              
            !**Jacobian**              
            elseif (programStep==3) then
              
            endif
            
          case('IMPLICIT') ! For fully implicit schemes
            !**RHS**
            if (programStep==1 .or.programStep==2) then
              call Burgers_dUdt(uvec,resI_vec,time,ep,dt)
              resE_vec(:)=0.0_wp

            !**Jacobian**
            elseif (programStep==3) then

              call Allocate_CSR_Storage(vecl,nnz_D2)
              call Build_Jac(uvec,ep,dt,akk,time)          

            endif
            
          case default ! To catch invald inputs
            write(*,*)'Invaild case entered. Enter "EXPLICIT", or "IMPLICIT"'
            write(*,*)'Exiting'
            stop
            
        end select
        
      endif
      
      end subroutine Burgers
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

      return
      end subroutine grid
!==============================================================================
! RETURNS VECTOR WITH EXACT SOLUTION
      subroutine exact_Burg(u,eps,time)

      real(wp),                  intent(in   ) :: eps,time
      real(wp), dimension(vecl), intent(  out) :: u
      integer                                  :: i

      do i = 1,vecl
        u(i) =  NL_Burg_exactsolution(x(i),time,eps)
      enddo

      return
      end subroutine exact_Burg

!==============================================================================
! DEFINES RHS 
      subroutine Burgers_dUdt(u,dudt,time,eps,dt)
      
      use SBP_Coef_Module, only: Pinv,D1,D2,jD1,jD2,iD1,iD2 
      use matvec_module,   only: amux

      real(wp), dimension(vecl), intent(in   ) ::  u
      real(wp), dimension(vecl), intent(  out) :: dudt
      real(wp),                  intent(in   ) :: time, eps, dt

      real(wp), dimension(vecl)                :: f, df, dfv, gsat, wrk    
      integer                                  :: i
      real(wp)                                 :: uL,uR,du0,du1,a0,a1,g0,g1
      
      gsat = 0.0_wp ; df  = 0.0_wp ; dfv = 0.0_wp ;

!--------------------------FLUX------------------------------------------------
      !     Canonical splitting for Inviscid flux in Burgers eqn
      !     f(u)_x = 2/3 (u u/2)_x + 1/3 u u_x

      f(:) = half*u(:)*u(:)

      call amux(vecl,f,wrk,D1,jD1,iD1)

      df(:) = df(:) + twothird*wrk(:)

      call amux(vecl,u,wrk,D1,jD1,iD1)
      df(:) = df(:) + third*u(:)*wrk(:)

      !  Viscous Flux in Burgers eqn

      call amux(vecl,u,dfv,D2,jD2,iD2)

!--------------Inflow Boundary Condition---------------------------------------

      uR = u(1)

      a0 = third*(uR + sqrt(uR**2+1.0e-28_wp))    !a0 = third*(uR + abs(uR))
      g0 = a0 * NL_Burg_exactsolution(x(1),time,eps)                &
         - eps* NL_Burg_exact_derivative(x(1),time,eps)

      du0 = dot_product(d1vec0(1:4),u(1:4)) / dx

      gsat(1) = gsat(1) + sig0 * Pinv(1) * (a0*uR - eps*du0 - g0) 

!---------------Outflow Boundary Condition-------------------------------------

      uL = u(vecl)

      a1 = third*(uL - sqrt(uL**2+1.0e-28_wp)) !a1 = third*(uL - abs(uL))
      g1 = a1 * NL_Burg_exactsolution(x(vecl),time,eps)                &
         - eps* NL_Burg_exact_derivative(x(vecl),time,eps)
      du1 = dot_product(d1vec1(1:4),u(vecl-3:vecl)) / dx

      gsat(vecl) = gsat(vecl) + sig1 * Pinv(vecl) * (a1*uL-eps*du1-g1)

!--------------------Sum all terms---------------------------------------------

      dudt(:) = dt*(eps*dfv(:) - df(:) + gsat(:))

      end subroutine Burgers_dUdt
      
!==============================================================================
! CREATES JACOBIAN
      subroutine Build_Jac(u,eps,dt,akk,time)

      use SBP_Coef_Module,  only: Pinv,D1,D2,jD1,jD2,iD1,iD2       
      use matvec_module,    only: amux
      use unary_mod,        only: aplb,aplsca,amudia,diamua,apldia
      use Jacobian_CSR_Mod, only: iaJac,jaJac,aJac
      
      real(wp), dimension(vecl), intent(in) :: u
      real(wp),                  intent(in) :: eps,dt,akk,time    
        
      integer,  dimension(size(iaJac))  :: iwrk1,iwrk2,iwrk3,iJac
      integer,  dimension(size(jaJac))  :: jwrk1,jwrk2,jwrk3,jJac
      real(wp), dimension(size( aJac))  :: wrk1,wrk2,wrk3,Jac,wrk4,eps_d2      
      
      integer,  dimension(vecl)  :: iw
      real(wp), dimension(vecl)  :: diag
      integer,  dimension(2)     :: ierr
      real(wp)                   :: uL,a1_d,uR,a0_d,a0,a1
      integer                    :: nnz_max
!------------------------------------------------------------------------------ 
      nnz_max=size(aJac)
     

!---------------------dgsat/dt-------------------------------------------------
      uL=u(vecl)
      a1 = third*(uL - sqrt(uL**2+1.0e-28_wp))
      a1_d=third*(1-uL/sqrt(uL**2+1.0e-28_wp))
   
      uR=u(1)
      a0 = third*(uR + sqrt(uR**2+1.0e-28_wp))
      a0_d=third*(1+uR/sqrt(uR**2+1.0e-28_wp))

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
      call amudia(vecl,1,D1,jD1,iD1,u,wrk1,jwrk1,iwrk1)                   !1
      wrk1(:)=-twothird*wrk1(:)                                              !1
      
      call diamua(vecl,1,D1,jD1,iD1,u,wrk2,jwrk2,iwrk2)                   !2
      wrk2(:)=-third*wrk2(:)                                                 !2
      
      call aplb(vecl,vecl,1,wrk1,jwrk1,iwrk1,wrk2,jwrk2,iwrk2,wrk3,&   !3
     &          jwrk3,iwrk3,nnz_max,iw,ierr(1))                               !3

      eps_d2=eps*D2(:)                                                       !4
      call aplb(vecl,vecl,1,eps_d2,jD2,iD2,wrk3,jwrk3,iwrk3,Jac,&     !5a
     &          jJac,iJac,nnz_max,iw,ierr(2))                                !5a

      call amux(vecl,-third*u,diag,D1,jD1,iD1) !First-derivative of u vec!5b
      call apldia(vecl,0,Jac,jJac,iJac,diag,Jac,jJac,iJac,iw)            !5b

     ! R - 5c
      Jac(1)=Jac(1)+sig0*Pinv(1)*(a0+a0_d*(uR- &
     &       NL_Burg_exactsolution(x(1),time,eps))) 
      Jac(1:4)=Jac(1:4)+sig0*Pinv(1)*(-eps)*d1vec0(:)/dx
      
      ! L - 5c
      Jac(nnz_max)=Jac(nnz_max)+sig1*Pinv(vecl)*(a1+a1_d*(uL- &
     &            NL_Burg_exactsolution(x(vecl),time,eps))) 
      Jac(nnz_max-3:nnz_max)=Jac(nnz_max-3:nnz_max)+ &
     &                     sig1*Pinv(vecl)*(-eps)*d1vec1(:)/dx
     
      wrk4=-akk*dt*Jac                                                       !6
      call aplsca(vecl,wrk4,jJac,iJac,1.0_wp,iw)                          !7

      !8
      iaJac=iJac
      jaJac=jJac
      aJac=wrk4
!-----------------------------------------------------------------------------  
      if (sum(ierr)/=0) then ! Catch errors
        print*,'Error building Jacobian'
        stop
      endif
   
      end subroutine Build_Jac
      
!==============================================================================
! 
      subroutine error(u,time,eps,dt)

      use SBP_Coef_Module, only: Pmat

      real(wp),                     intent(in   ) :: time, eps, dt
      real(wp), dimension(vecl), intent(inout) :: u

      integer                                    :: i
      real(wp)                                   :: errL2,errLinf,wrk,psum

!     calculate the rms residual over the domain                        

      psum    = 0.0_wp
      errL2   = 0.0_wp
      errLinf = 0.0_wp
      do i=1,vecl
              wrk = abs(u(i) - NL_Burg_exactsolution(x(i),time,eps))
            errL2 =  errL2 +  pmat(i)*wrk * wrk
          errLinf =  max(errLinf,wrk)
             psum = psum + pmat(i)
      enddo
      errL2  = ( sqrt(errL2/vecl/psum) )
      errLinf= ( errLinf)
!!      write( *,89)ixd-1,errL2,errLinf
   89 format('P=',I4,' ,,  L2',e18.10,' ,,  Linf',e18.10,' ,,')

      end subroutine error
!==============================================================================
!
      subroutine plot(punit,u,time,eps)

      integer,                      intent(in   ) :: punit
      real(wp),                     intent(in   ) :: time, eps
      real(wp), dimension(vecl), intent(inout) :: u

      integer                                    :: i
      real(wp)                                   :: wrk,err

!     write to plotter file                                             
      do i=1,vecl
        wrk = NL_Burg_exactsolution(x(i),time,eps)
        err = ( abs( wrk - u(i) ) + 1.0e-15 )
        write(punit,2)x(i),wrk,u(i),err
      enddo
    2 format(4(1x,e15.7))

      end subroutine plot
!==============================================================================
! DEFINES EXACT SOLUTION FOR PARTICULAR x AND t VALUE
      function NL_Burg_exactsolution(xin,tin,eps)

      real(wp), intent(in) :: xin, tin, eps
      real(wp), parameter  ::  a = -0.40_wp, b = 1.0_wp , d = +0.50_wp
      real(wp)             ::  c = (a+b)/2

      real(wp) :: NL_Burg_exactsolution

      integer  :: i
      real(wp) :: t1,t2, x0

      select case (exact_solution)
        case('Olver')         !          pp. 1190   Peter J. Olver

          t1  = (xin - c*tin - d )
          t2  = exp((b-a)*t1/(2.0_wp*eps))

          NL_Burg_exactsolution = (a * t2  + b ) / (1.0_wp * t2 + 1.0_wp)

        case('tanh')

          x0 =  0.1_wp
          NL_Burg_exactsolution = 1.1_wp - tanh(0.5_wp*(xin-1.1_wp*tin-x0)/eps)

      end select

      end function NL_Burg_exactsolution
!==============================================================================
! DEFINES EXACT DERIVATIVE FOR PARTICULAR x AND t VALUE
      function NL_Burg_exact_derivative(xin,tin,eps)

       real(wp), intent(in) :: xin, tin, eps
       real(wp), parameter  ::  a = -0.40_wp, b = 1.0_wp , d = +0.50_wp
       real(wp)             ::  c = (a+b)/2

       real(wp) :: NL_Burg_exact_derivative

       integer  :: i
       real(wp) :: t1,t2,t3, x0

       select case (exact_solution)
         case('Olver')         !          pp. 1190   Peter J. Olver

           t1  = (xin - c*tin - d )
           t2  = ((a - b)*t1)/(4.0_wp*eps)
           t3  = 2.0_wp / (exp(t2)+exp(-t2))

           NL_Burg_exact_derivative =  -((a - b)**2* t3**2) /(8.0_wp*eps)

         case('tanh')

           x0 =  0.1_wp
           NL_Burg_exact_derivative =  -1.0_wp/cosh((-1.1_wp*tin+xin-x0)/ &
     &                                 (2.0_wp*eps))**2 / (2.0_wp * eps)

       end select

      end function NL_Burg_exact_derivative
!==============================================================================
      end module Burgers_Module
