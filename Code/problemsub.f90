!******************************************************************************
! Subroutine to take in the problem number and all the appropriate problem
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90         *ONTAINS VARIABLES USED IN THE PROGRAM
! RUNGE_KUTTA.F90               *CONTAINS RK CONSTANTS
!******************************************************************************
      module problemsub_mod
      
      use precision_vars,    only: wp      
      
      implicit none;save
      
      private
      public :: problemsub
      
      contains
      
      subroutine problemsub(iprob,programStep,nveclen,ep,&
     &                      dt,tfinal,iDT,time,akk,L)

      use control_variables, only: resE,resI,uvec
      use Burgers_Module,    only: Burgers
   
      !PROBLEM PARAMETERS
      integer,  intent(in   ) :: iprob, programStep
      integer,  intent(inout) :: nveclen
      real(wp), intent(in   ) :: ep
      real(wp), intent(inout) :: dt
      real(wp), intent(  out) :: tfinal
      integer,  intent(in   ) :: iDT
      real(wp), intent(in   ) :: time
      real(wp), intent(in   ) :: akk
      integer,  intent(in   ) :: L
      real(wp), dimension(size(uvec)) :: resE_vec,resI_vec

      if     (iprob==1) then
       call vanderPol(   programStep,nveclen,ep,&
     &                dt,tfinal,iDT,resE_vec,resI_vec,akk)
      elseif (iprob==2) then
       call Pureschi(    programStep,nveclen,ep,&
     &                dt,tfinal,iDT,resE_vec,resI_vec,akk)
      elseif (iprob==3) then
       call Kaps(        programStep,nveclen,ep,&
     &                dt,tfinal,iDT,resE_vec,resI_vec,akk)
      elseif (iprob==4) then
       call Kreiss(      programStep,nveclen,ep,&
     &                dt,tfinal,iDT,time,resE_vec,resI_vec,akk)
      elseif (iprob==5) then 
       call Lorenz(       programStep,nveclen,ep,&
     &                dt,tfinal,iDT,resE_vec,resI_vec,akk)
      elseif (iprob==6) then 
       call Rossler_Chaos(programStep,nveclen,ep,&
     &                dt,tfinal,iDT,resE_vec,resI_vec,akk)
      elseif (iprob==7) then 
       call Oregonator(   programStep,nveclen,ep,&
     &                dt,tfinal,iDT,resE_vec,resI_vec,akk)
      elseif (iprob==8) then 
       call Brusselator(   programStep,nveclen,ep,&
     &                dt,tfinal,iDT,resE_vec,resI_vec,akk)
      elseif (iprob==9) then
       call Burgers(       programStep,nveclen,ep,&
     &                dt,tfinal,iDT,time,resE_vec,resI_vec,akk)
      endif
      
      if(programStep >= 0) then
  
        resE(:,L)=resE_vec(:)
        resI(:,L)=resI_vec(:)
      endif

      end subroutine problemsub
      end module problemsub_mod
