!******************************************************************************
! Subroutine to take in the problem number and call the appropriate problem
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90     *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90  *CONTAINS VARIABLES USED IN THE PROGRAM
! VANDERPOL.F90          *PROBLEM CONSTANTS FOR VANDERPOL
! PURESCHI.F90           *PROBLEM CONSTANTS FOR PURESCHI & RUSSO
! KAPS.F90               *PROBLEM CONSTANTS FOR KAPS
! KREISS.F90             *PROBLEM CONSTANTS FOR KREISS'
! ROSSLER_CHAOS.F90      *PROBLEM CONSTANTS FOR Rossler_Chaos
! OREGONATOR.F90         *PROBLEM CONSTANTS FOR OREGONATOR
! BRUSSELATOR.F90        *PROBLEM CONSTANTS FOR BRUSSELATOR 
! BURGERS_MOD.F90        *PROBLEM CONSTANTS AND ROUTINES FOR BURGERS
! BOSCARINO31_MOD.F90    *PROBLEM CONSTANTS AND ROUTINES FOR BOSCARINO-31
!******************************************************************************
      module problemsub_mod
      
      use precision_vars,    only: wp      
      use control_variables, only: programstep
      
      implicit none;save
      
      private
      public :: problemsub
      
      contains
      
      subroutine problemsub(iprob,nveclen,ep,dt,tfinal,iDT,time,akk,L)

      use control_variables, only: resE,resI,uvec
      use vanderPol_mod,     only: vanderPol
      use Pureschi_mod,      only: Pureschi
      use Kaps_mod,          only: Kaps
      use Kreiss_mod,        only: Kreiss
      use Lorenz_mod,        only: Lorenz
      use Rossler_mod,       only: Rossler_Chaos
      use Oregonator_mod,    only: Oregonator
      use Brusselator_mod,   only: Brusselator
      use Burgers_Module,    only: Burgers
      use Boscarino31_Mod,   only: Boscarino31
      use Broadwell_Mod,     only: Broadwell
   
      !PROBLEM PARAMETERS
      integer,  intent(in   ) :: iprob
      integer,  intent(  out) :: nveclen
      real(wp), intent(in   ) :: ep
      real(wp), intent(inout) :: dt
      real(wp), intent(  out) :: tfinal
      integer,  intent(in   ) :: iDT
      real(wp), intent(in   ) :: time
      real(wp), intent(in   ) :: akk
      integer,  intent(in   ) :: L
      real(wp), dimension(size(uvec)) :: resE_vec,resI_vec

      if     (iprob==1)  then
       call vanderPol(    nveclen,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==2)  then
       call Pureschi(     nveclen,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==3)  then
       call Kaps(         nveclen,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==4)  then
       call Kreiss(       nveclen,ep,dt,tfinal,iDT,time,resE_vec,resI_vec,akk)
      elseif (iprob==5)  then 
       call Lorenz(       nveclen,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==6)  then 
       call Rossler_Chaos(nveclen,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==7)  then 
       call Oregonator(   nveclen,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==8)  then 
       call Brusselator(  nveclen,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==9)  then
       call Burgers(      nveclen,ep,dt,tfinal,iDT,time,resE_vec,resI_vec,akk)
      elseif (iprob==10) then
       call Boscarino31(  nveclen,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==11) then
       call Broadwell(    nveclen,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      else
       print*,'Invalid problem number!'
       stop
      endif
      
      select case(programstep)
        case('BUILD_RHS')  
          resE(:,L)=resE_vec(:)
          resI(:,L)=resI_vec(:)
      end select

      end subroutine problemsub

      end module problemsub_mod
