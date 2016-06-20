      subroutine problemsub(iprob,stage,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE,resI,akk,xjac)

      use precision_vars

      implicit none
    
      integer,  parameter                      :: vecl=4

      !PROBLEM PARAMETERS
      integer,                   intent(in   ) :: iprob, stage
      real(wp), dimension(vecl), intent(inout) :: uvec
      real(wp),                  intent(in   ) :: ep
      real(wp), dimension(vecl), intent(  out) :: uexact
      real(wp),                  intent(inout) :: dt
      integer,                   intent(inout) :: nveclen
      real(wp),                  intent(  out) :: tfinal
      integer,                   intent(in   ) :: iDT
      real(wp), dimension(vecl), intent(  out) :: resE,resI
      real(wp),                       intent(in   ) :: akk
      real(wp), dimension(vecl,vecl), intent(  out) :: xjac


      if     (iprob==1) then
        call vanderPol(stage,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE,resI,akk,xjac)
      elseif (iprob==2) then
        call Pureschi(stage,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE,resI,akk,xjac)
      elseif (iprob==3) then
        call Kaps(stage,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE,resI,akk,xjac)
      elseif (iprob==4) then
        call Kreiss(stage,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE,resI,akk,xjac)
      endif
      return
      end subroutine
