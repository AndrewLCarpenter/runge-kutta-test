      subroutine problemsub(iprob,problemStep,probname,nveclen,uvec,ep,uexact,dt,tfinal,iDT,resE,resI,akk,xjac)

      use precision_vars

      implicit none
    
      integer,  parameter                      :: vecl=2

      !PROBLEM PARAMETERS
      integer,                   intent(in   ) :: iprob, problemStep
      character(len=9), intent(  out) :: probname
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
        call vanderPol(problemStep,probname,nveclen,uvec,ep,uexact,dt,tfinal,iDT,resE,resI,akk,xjac)
      elseif (iprob==2) then
        call Pureschi(problemStep,probname,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE,resI,akk,xjac)
      elseif (iprob==3) then
        call Kaps(problemStep,probname,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE,resI,akk,xjac)
      elseif (iprob==4) then
        call Kreiss(problemStep,probname,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE,resI,akk,xjac)
      endif      
      return
      end subroutine
