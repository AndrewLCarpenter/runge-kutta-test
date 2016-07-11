      module Burgers_Module 

      use precision_vars
      use SBP_Coef_Module
      use unary_mod, only: aplb,aplsca,amudia,diamua 

      implicit none
      
      private
      public    ::  grid,exact_Burg,error,plot,dx,Burgers_dudt,Build_Jac

      
      real(wp)  :: dx
      character(120), parameter :: exact_solution = 'tanh'
!      character(120), parameter :: exact_solution = 'Oliver'
      
      contains

!==============================================================================

      subroutine grid(x,xL,xR,nveclen)

      implicit none

      integer,                      intent(in   ) :: nveclen
      real(wp), dimension(nveclen), intent(  out) :: x
      real(wp),                     intent(in   ) :: xL,xR

      integer                                     :: i

      do i=1,nveclen
        x(i)= xL + (xR-xL)*(i-1.0_wp)/(nveclen-1.0_wp)
      enddo

      return
      end subroutine grid

!==============================================================================

      subroutine exact_Burg(nveclen,x,u,eps,time)

      integer,                     intent(in)    :: nveclen
      real(wp),                    intent(in)    :: eps,time
      real(wp), dimension(nveclen),   intent(in)    :: x
      real(wp), dimension(nveclen),   intent(  out) :: u
      
      integer                                    :: i

      do i = 1,nveclen
        u(i) =  NL_Burg_exactsolution(x(i),time,eps)
      enddo

      return
      end subroutine exact_Burg

!==============================================================================

      subroutine error(nveclen,x,u,time,eps,dt)

      integer,                      intent(in)    :: nveclen
      real(wp),                     intent(in)    :: time, eps, dt
      real(wp), dimension(nveclen), intent(inout) :: x,u

      integer                                    :: i
      real(wp)                                   :: errL2,errLinf,wrk,psum

!     calculate the rms residual over the domain                        

      psum    = 0.0_wp
      errL2   = 0.0_wp
      errLinf = 0.0_wp
      do i=1,nveclen
              wrk = abs(u(i) - NL_Burg_exactsolution(x(i),time,eps))
            errL2 =  errL2 +  pmat(i)*wrk * wrk
          errLinf =  max(errLinf,wrk)
             psum = psum + pmat(i)
      enddo
      errL2  = ( sqrt(errL2/nveclen/psum) )
      errLinf= ( errLinf)
!!      write( *,89)ixd-1,errL2,errLinf
   89 format('P=',I4,' ,,  L2',e18.10,' ,,  Linf',e18.10,' ,,')

      return
      end subroutine error

!==============================================================================

      subroutine plot(punit,nveclen,x,u,time,eps)

      integer,                      intent(in)    :: punit,nveclen
      real(wp),                     intent(in)    :: time, eps
      real(wp), dimension(nveclen), intent(inout) :: x,u

      integer                                    :: i
      real(wp)                                   :: wrk,err

!     write to plotter file                                             

      do i=1,nveclen
        wrk = NL_Burg_exactsolution(x(i),time,eps)
        err = ( abs( wrk - u(i) ) + 1.0e-15 )
        write(punit,2)x(i),wrk,u(i),err
      enddo
    2 format(4(1x,e15.7))
      return
      end subroutine plot

!==============================================================================

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
           NL_Burg_exact_derivative =  -1.0_wp/cosh((-1.1_wp*tin+xin-x0)/(2.0_wp*eps))**2 / (2.0_wp * eps)

       end select

      end function NL_Burg_exact_derivative


!==============================================================================

      subroutine Burgers_dUdt(nveclen,x,u,dudt,time,eps,dt)

      integer,                      intent(in   ) :: nveclen
      real(wp), dimension(nveclen), intent(in   ) :: x,u
      real(wp), dimension(nveclen), intent(  out) :: dudt
      real(wp),                     intent(in   ) :: time, eps, dt

      real(wp), dimension(nveclen)                :: f, df, dfv, gsat,wrk
      
      integer                                     :: i
      real(wp)                                    ::  uL,  uR, dx
      real(wp)                                    :: du0, du1
      real(wp)                                    ::  a0,  a1
      real(wp)                                    ::  g0,  g1
      real(wp),  parameter                        :: sig0 = -1.0_wp
      real(wp),  parameter                        :: sig1 = +1.0_wp

      real(wp), dimension(4), parameter           :: d1vec0= + reshape((/-24.0_wp/17.0_wp,  &
                                                                      &  +59.0_wp/34.0_wp,  &
                                                                      &   -4.0_wp/17.0_wp,  &
                                                                      &   -3.0_wp/34.0_wp/),&
                                                                      & (/4/))
      real(wp), dimension(4), parameter           :: d1vec1= + reshape((/+ 3.0_wp/34.0_wp,  &
                                                                      &  + 4.0_wp/17.0_wp,  &
                                                                      &  -59.0_wp/34.0_wp,  &
                                                                      &  +24.0_wp/17.0_wp/),&
                                                                      & (/4/))

      dx=x(2)-x(1)
!      print*,'calling in dudt'
      call Define_CSR_Operators(nveclen,dx)

      gsat = 0.0_wp ; df  = 0.0_wp ; dfv = 0.0_wp ;

      !     Canonical splitting for Inviscid flux in Burgers eqn
      !     f(u)_x = 2/3 (u u/2)_x + 1/3 u u_x

      f(:) = half*u(:)*u(:)

      call amux_local(nveclen,f,wrk,D1,jD1,iD1)

      df(:) = df(:) + twothirds*wrk(:)

      call amux_local(nveclen,u,wrk,D1,jD1,iD1)
      df(:) = df(:) + onethird*u(:)*wrk(:)

      !  Viscous Flux in Burgers eqn

      call amux_local(nveclen,u,dfv,D2,jD2,iD2)

      !  Inflow Boundary Condition

      uL = NL_Burg_exactsolution(x(1),time,eps)
      uR = u(1)

!      a0 = third*(uR + abs(uR))
      a0 = third*(uR + sqrt(uR**2+1.0e-28_wp))
      g0 = a0 * NL_Burg_exactsolution(x(1),time,eps)                &
         - eps* NL_Burg_exact_derivative(x(1),time,eps)

      du0 = dot_product(d1vec0(1:4),u(1:4)) / dx
      gsat(1) = gsat(1) + sig0 * Pinv(1) * (a0*uR - eps*du0 - g0)

      !  Outflow Boundary Condition

      uL = u(nveclen)
      uR = NL_Burg_exactsolution(x(nveclen),time,eps)

!      a1 = third*(uL - abs(uL))
      a1 = third*(uL - sqrt(uL**2+1.0e-28_wp))
      g1 = a1 * NL_Burg_exactsolution(x(nveclen),time,eps)                &
         - eps* NL_Burg_exact_derivative(x(nveclen),time,eps)
      du1 = dot_product(d1vec1(1:4),u(nveclen-3:nveclen)) / dx

      gsat(nveclen) = gsat(nveclen) + sig1 * Pinv(nveclen) * (a1*uL - eps*du1 - g1)
  !    gsat(:)=0.0_wp
      !  Sum all terms

      dudt(:) = dt*(eps*dfv(:) - df(:) + gsat(:))

      !dudt(:)=-df(:)
!      print*,'de in dudt'
      deallocate(iD1,iD2,Pmat,Pinv,jD1,jD2,D1,D2)     

      end subroutine Burgers_dUdt
!==============================================================================
      subroutine Build_Jac(nveclen,u,x,eps,dt,akk,iaxJac,jaxJac,axJac)
      
      integer,                         intent(in   ) :: nveclen
      real(wp), dimension(nveclen),    intent(in   ) :: u,x
      real(wp),                        intent(in   ) :: eps,dt,akk
      integer,  dimension(nveclen+1),  intent(  out) :: iaxJac
      integer,  dimension(nveclen**2), intent(  out) :: jaxJac
      real(wp), dimension(nveclen**2), intent(  out) :: axJac
      
      integer,  dimension(nveclen+1)  :: iwrk1,iwrk2,iwrk3,iJac
      integer,  dimension(nveclen**2) :: jwrk1,jwrk2,jwrk3,jJac
      real(wp), dimension(nveclen**2) :: wrk1,wrk2,wrk3a,wrk3b,Jac,wrk4,eps_d2      
      
      integer,  dimension(nveclen)                :: iw
      integer,  dimension(2)                      :: ierr
      real(wp)                                    :: dx
      integer :: i

      dx=x(2)-x(1)

      call Define_CSR_Operators(nveclen,dx)

!      jac=eps*d2-(onethird)*(d1*u+u*d1)  ||   +dgsat/du => not included yet
!      wrk1=d1*u
!      wrk2=u*d1
!      wrk3a=wrk1+wrk2
!      wrk3b=-1/3*wrk3a
!      eps_d2=eps*d2
!      jac=eps_d2+wrk3b

!      wrk4=-akk*dt*jac
!      xjac=I-akk*dt*jac=-akk*dt*jac+(1)*I
!      xjac=wrk4+(1)*I
! 
!      xjac-> wrk4, iJac, jJac

       call amudia(nveclen,1,D1,jD1,iD1,u,wrk1,jwrk1,iwrk1)
       call diamua(nveclen,1,D1,jD1,iD1,u,wrk2,jwrk2,iwrk2)
       call aplb(nveclen,nveclen,1,wrk1,jwrk1,iwrk1,wrk2,jwrk2,iwrk2,wrk3a,&
     &           jwrk3,iwrk3,nveclen**2,iw,ierr(1))
       wrk3b=-onethird*wrk3a(:)
       eps_d2=eps*D2(:)
       call aplb(nveclen,nveclen,1,eps_d2,jD2,iD2,wrk3b,jwrk3,iwrk3,Jac,&
     &           jJac,iJac,nveclen**2,iw,ierr(2))
     
       wrk4=-akk*dt*Jac
       call aplsca(nveclen,wrk4,jJac,iJac,1.0_wp,iw)

       iaxJac=iJac
 !      iaxJac(:)=0
   !    do i=1,17

    !    iaxJac(i)=0

       jaxJac=jJac
       axJac=wrk4

       if (sum(ierr)/=0) then
         print*,'Error building Jacobian'
         stop
       endif
      
!      print*,'de in jac'
      deallocate(iD1,iD2,Pmat,Pinv,jD1,jD2,D1,D2)  
      
      end subroutine Build_Jac
!==============================================================================
      end module Burgers_Module
