!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Burgers problem 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
! BURGERS_MOD.F90           *CONTAINS ROUTINES TO BUILD BURGER
!******************************************************************************

      subroutine Burgers(programStep,nveclen,ep,dt, &
     &                   tfinal,iDT,time,resE_vec,resI_vec,akk)

      use precision_vars
      use control_variables
      use Burgers_Module
      use Jacobian_CSR_Mod

      implicit none
!-----------------------VARIABLES----------------------------------------------
      integer,  parameter    :: vecl=16
      integer, intent(in   ) :: programStep

      !INIT vars
      real(wp), parameter :: xL=0.0_wp,xR=1.0_wp
      real(wp), dimension(vecl) :: x

      
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
        nvecLen = vecl
        probname='Burgers  '     
        tol=1.0e-12_wp  
        dt_error_tol=5.0e-14_wp
        
      !**Initialization of problem information**        
      elseif (programStep==0) then

        call grid(x,xL,xR,vecl)
        
        tinitial=0.0_wp
        call exact_Burg(vecl,x,uvec,ep,tinitial)
           
        dt = 0.25_wp*0.0001_wp/10**((iDT-1)/20.0_wp) ! timestep   
        dx = x(2)-x(1)
        tfinal = 0.5_wp                   ! final time
        
        call exact_Burg(vecl,x,uexact,ep,tfinal)
              
      !**RHS and Jacobian**
      elseif (programStep>=1) then

        select case (Temporal_Splitting)
        
          case('EXPLICIT')
            !**RHS**
            if (programStep==1 .or.programStep==2) then
              call Burgers_dUdt(vecl,x,uvec,resE_vec,time,ep,dt)
              
              resI_vec(:)=0.0_wp
            !**Jacobian**              
            elseif (programStep==3) then
              
            endif

          case('IMEX') ! For IMEX schemes
            !**RHS**
            if (programStep==1 .or.programStep==2) then

            !**Jacobian**
            elseif (programStep==3) then

            endif
            
          case('IMPLICIT') ! For fully implicit schemes
            !**RHS**
            if (programStep==1 .or.programStep==2) then
              call Burgers_dUdt(vecl,x,uvec,resI_vec,time,ep,dt)
              resE_vec(:)=0.0_wp

            !**Jacobian**
            elseif (programStep==3) then
              jac_case='SPARSE'
              call Allocate_CSR_Storage(vecl)
              call Build_Jac(vecl,uvec,x,ep,dt,akk,iaJac,jaJac,aJac)            

            endif
            
          case default ! To catch invald inputs
            write(*,*)'Invaild case entered. Enter "EXPLICIT", "IMEX", or "IMPLICIT"'
            write(*,*)'Exiting'
            stop
            
        end select
        
      endif
      
      end subroutine Burgers
