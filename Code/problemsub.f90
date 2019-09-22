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
! BRUSSELATOR_2.F90      *PROBLEM CONSTANTS FOR BRUSSELATOR 
! BURGERS_MOD.F90        *PROBLEM CONSTANTS AND ROUTINES FOR BURGERS
! BOSCARINO31_MOD.F90    *PROBLEM CONSTANTS AND ROUTINES FOR BOSCARINO-31
! BROADWELL_MOD.F90      *PROBLEM CONSTANTS AND ROUTINES FOR BROADWELL
! CHARNEY_DeVORE6.F90    *PROBLEM CONSTANTS AND ROUTINES FOR CHARNEY_DeVore 6 eqn model
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From control_variables:
!   programstep -> string used in program_step_select, character(len=80), not modified 
! From Problem modules:
!   Use each problem's main subroutine to get problem constants, RHS and Jacobian
!******************************************************************************
! INPUTS:
! iprob   -> defines problem number,                                    integer
! ep      -> Stiffness epsilon value,                                   real(wp)
! iDT     -> Timestep counter from timestep loop to define dt,          integer
! time    -> Current solution time,                                     real(wp)
! akk     -> Diagonal term from RK scheme,                              real(wp)
! L       -> stage number,                                              integer
! INOUTS:
! dt      -> timestep: created and then used later in the problem case, real(wp)
! OUTPUTS:
! nveclen -> total length of u-vector used for problem                  integer
! neq     -> number of equations in problem                             integer
! tfinal  -> final time for iterative solver                            real(wp)
!******************************************************************************
      module problemsub_mod
      
      use precision_vars,    only: wp      
      use control_variables, only: programstep
      
      implicit none;save
      
      private
      public :: problemsub
      
      contains
      
      subroutine problemsub(iprob,nveclen,neq,ep,dt,tfinal,iDT,time,akk,L)

      use control_variables,       only: resE,resI,uvec,               &
                                       & resE_Tens_OTD,resI_Tens_OTD,  &
                                       & resE_OTD,resI_OTD
      use vanderPol_mod,           only: vanderPol
      use Pureschi_mod,            only: Pureschi
      use Kaps_mod,                only: Kaps
      use Kreiss_mod,              only: Kreiss
      use Lorenz_mod,              only: Lorenz
      use Rossler_mod,             only: Rossler_Chaos
      use Oregonator_mod,          only: Oregonator
      use Brusselator_mod,         only: Brusselator
      use Burgers_Module,          only: Burgers
      use Boscarino31_Mod,         only: Boscarino31
      use Broadwell_Mod,           only: Broadwell
      use Charney_DeVore6_mod,     only: Charney_DeVore6
      use Charney_DeVore10_mod,    only: Charney_DeVore10
      use Kuramoto_Sivashinsky_mod,only: Kuramoto_Sivashinsky
      use Boldrighini_Franceschini_mod,only: Boldrighini_Franceschini
   
      !PROBLEM PARAMETERS
      integer,  intent(in   ) :: iprob
      integer,  intent(  out) :: nveclen
      integer,  intent(  out) :: neq
      real(wp), intent(in   ) :: ep
      real(wp), intent(inout) :: dt
      real(wp), intent(  out) :: tfinal
      integer,  intent(in   ) :: iDT
      real(wp), intent(in   ) :: time
      real(wp), intent(in   ) :: akk
      integer,  intent(in   ) :: L
      real(wp), dimension(size(uvec)) :: resE_vec,resI_vec

      if     (iprob==1)  then
       call vanderPol(      nveclen,neq,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==2)  then
       call Pureschi(       nveclen,neq,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==3)  then
       call Kaps(           nveclen,neq,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==4)  then
       call Kreiss(         nveclen,neq,ep,dt,tfinal,iDT,time,resE_vec,resI_vec,akk)
      elseif (iprob==5)  then 
       call Lorenz(         nveclen,neq,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==6)  then 
       call Rossler_Chaos(  nveclen,neq,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==7)  then 
       call Oregonator(     nveclen,neq,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==8)  then 
       call Brusselator(    nveclen,neq,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==9)  then
       call Burgers(        nveclen,neq,ep,dt,tfinal,iDT,time,resE_vec,resI_vec,akk)
      elseif (iprob==10) then
       call Boscarino31(    nveclen,neq,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==11) then
       call Broadwell(      nveclen,neq,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==12) then
       call Charney_DeVore6(nveclen,neq,ep,dt,tfinal,iDT,     resE_vec,resI_vec,akk)
      elseif (iprob==13) then
       call Charney_DeVore10(nveclen,neq,ep,dt,tfinal,iDT,    resE_vec,resI_vec,akk)
      elseif (iprob==14) then
       call Kuramoto_Sivashinsky(nveclen,neq,ep,dt,tfinal,iDT,    resE_vec,resI_vec,akk)
      elseif (iprob==15) then
       call Boldrighini_Franceschini(nveclen,neq,ep,dt,tfinal,iDT,    resE_vec,resI_vec,akk)
      else
       print*,'Invalid problem number!'
       stop
      endif
      
      select case(programstep)
        case('BUILD_RHS')  
          resE(:,L)=resE_vec(:)
          resI(:,L)=resI_vec(:)
        case('BUILD_RHS_OTD')  
          resE_Tens_OTD(:,:,L) = resE_OTD(:,:)
          resI_Tens_OTD(:,:,L) = resI_OTD(:,:)
      end select

      end subroutine problemsub

      end module problemsub_mod
