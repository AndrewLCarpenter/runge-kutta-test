      subroutine problemsub(iprob,programStep,nveclen, &
     &      uvec,ep,uexact,dt,tfinal,iDT,time,resE,resI,akk,xjac)

      use precision_vars
      use control_variables

      implicit none
    
      !integer,  parameter                      :: vecl=2

      !PROBLEM PARAMETERS
      integer,                   intent(in   ) :: iprob, programStep
      integer,                   intent(inout) :: nveclen
      real(wp), dimension(nveclen), intent(inout) :: uvec
      real(wp),                  intent(in   ) :: ep
      real(wp), dimension(nveclen), intent(  out) :: uexact
      real(wp),                  intent(inout) :: dt
      real(wp),                  intent(  out) :: tfinal
      integer,                   intent(in   ) :: iDT
      real(wp),                  intent(inout) :: time
      real(wp), dimension(nveclen), intent(  out) :: resE,resI
      real(wp),                       intent(in   ) :: akk
      real(wp), dimension(nveclen,nveclen), intent(  out) :: xjac
      



      

      if     (iprob==1) then
       call vanderPol(programStep,nveclen,uvec,ep,&
     &                uexact,dt,tfinal,iDT,resE,resI,akk,xjac)
      elseif (iprob==2) then
       call Pureschi( programStep,nveclen,uvec,ep,&
     &                uexact,dt,tfinal,iDT,resE,resI,akk,xjac)
      elseif (iprob==3) then
       call Kaps(     programStep,nveclen,uvec,ep,&
     &                uexact,dt,tfinal,iDT,resE,resI,akk,xjac)
      elseif (iprob==4) then
       call Kreiss(   programStep,nveclen,uvec,ep,&
     &                uexact,dt,tfinal,iDT,time,resE,resI,akk,xjac)
      elseif (iprob==5) then 
       call Lorenz(   programStep,nveclen,uvec,ep,&
     &                uexact,dt,tfinal,iDT,resE,resI,akk,xjac)
      endif

      return
      end subroutine
