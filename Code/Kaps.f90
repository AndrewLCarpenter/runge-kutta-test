!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Kaps problem (Dekker 7.5.2 pp. 215 (Kaps problem   : Index 1)) 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
      module Kaps_mod
      private
      public :: Kaps
      contains
      subroutine Kaps(programStep,nveclen,ep,dt, &
     &                tfinal,iDT,resE_vec,resi_vec,akk)

      use precision_vars,    only: wp
      use control_variables, only: temporal_splitting,probname,xjac, &
     &                             tol,dt_error_tol,uvec,uexact

      implicit none; save
!-----------------------VARIABLES----------------------------------------------
      integer, parameter     :: vecl=2
      integer, intent(in   ) :: programStep

      !INIT vars
      real(wp), intent(in   ) :: ep
      real(wp), intent(inout) :: dt
      integer,  intent(inout) :: nveclen
      real(wp), intent(  out) :: tfinal
      integer,  intent(in   ) :: iDT

      real(wp)                :: tmp, epI

      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: resE_vec,resi_vec
      
      !Jacob vars
      real(wp), intent(in   ) :: akk
!------------------------------------------------------------------------------

      !**Pre-initialization. Get problem name and vector length**
      if (programStep==-1) then
        nvecLen = vecl
        probname='Kaps     '  
        tol=1.0e-12_wp
        dt_error_tol=1.0e-11_wp
              
      !**Initialization of problem information**
      elseif (programStep==0) then
       
        dt = 0.5_wp/10**((iDT-1)/20.0_wp) ! timestep
        tfinal = 1.0_wp                   ! final time

        !*Exact Solution**        
        tmp = exp(-tfinal)
        uexact(1) = tmp*tmp
        uexact(2) = tmp

        !**Equilibrium IC**
        uvec(1) = 1.0_wp
        uvec(2) = 1.0_wp

      !**RHS and Jacobian**
      elseif (programStep>=1) then
        epI = 1.0_wp / ep  !**Initialize 1/epsilon

        select case (Temporal_Splitting)
      
          case('IMEX') ! For IMEX schemes
            !**RHS**
            if (programStep==1 .or.programStep==2) then
              resE_vec(1) = dt*(-2.0_wp*uvec(1))
              resE_vec(2) = dt*(uvec(1) - uvec(2) - uvec(2)*uvec(2) )
              resi_vec(1) = dt*(-epI*uvec(1) + epI*uvec(2)*uvec(2))
              resi_vec(2) = 0.0_wp
            !**Jacobian**
            elseif (programStep==3) then
              xjac(1,1) = 1.0_wp-akk*dt*(-epI)
              xjac(1,2) = 0.0_wp-akk*dt*( epI)*2.0_wp*uvec(2)
              xjac(2,1) = 0.0_wp-akk*dt*(0.0_wp)
              xjac(2,2) = 1.0_wp-akk*dt*(0.0_wp)
            endif
            
          case('IMPLICIT') ! For fully implicit schemes
            !**RHS**
            if (programStep==1 .or.programStep==2) then
              resE_vec(:) = 0.0_wp
              resi_vec(1) = dt*(-(epI+2.0_wp)*uvec(1) + epI*uvec(2)*uvec(2))
              resi_vec(2) = dt*(uvec(1) - uvec(2) - uvec(2)*uvec(2) )
            !**Jacobian**
            elseif (programStep==3) then
              xjac(1,1) = 1.0_wp-akk*dt*(-(epI+2.0_wp))
              xjac(1,2) = 0.0_wp-akk*dt*(+epI*2.0_wp*uvec(2))
              xjac(2,1) = 0.0_wp-akk*dt*(1.0_wp)
              xjac(2,2) = 1.0_wp-akk*dt*(-(1.0_wp+2.0_wp*uvec(2)))
            endif
      
          case default ! To catch invald inputs
            write(*,*)'Invaild case entered. Enter "IMEX" or "IMPLICIT"'
            write(*,*)'Exiting'
            stop
            
        end select
          
      endif

      end subroutine Kaps
      end module Kaps_mod
