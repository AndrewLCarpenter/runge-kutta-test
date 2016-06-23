      subroutine problemsub(iprob,programStep,probname,nveclen, &
     &      temporal_splitting,uvec,ep,uexact,dt,tfinal,iDT,resE,resI,akk,xjac)

      use precision_vars

      implicit none
    
      !integer,  parameter                      :: vecl=2

      !PROBLEM PARAMETERS
      integer,                   intent(in   ) :: iprob, programStep
      character(len=9),          intent(  out) :: probname
      integer,                   intent(inout) :: nveclen
      character(80),             intent(in   ) :: temporal_splitting
      real(wp), dimension(nveclen), intent(inout) :: uvec
      real(wp),                  intent(in   ) :: ep
      real(wp), dimension(nveclen), intent(  out) :: uexact
      real(wp),                  intent(inout) :: dt
      real(wp),                  intent(  out) :: tfinal
      integer,                   intent(in   ) :: iDT
      real(wp), dimension(nveclen), intent(  out) :: resE,resI
      real(wp),                       intent(in   ) :: akk
      real(wp), dimension(nveclen,nveclen), intent(  out) :: xjac

      if     (iprob==1) then
       call vanderPol(programStep,probname,nveclen,temporal_splitting,uvec,ep,&
     &                uexact,dt,tfinal,iDT,resE,resI,akk,xjac)
      elseif (iprob==2) then
       call Pureschi( programStep,probname,nveclen,temporal_splitting,uvec,ep,&
     &                uexact,dt,tfinal,iDT,resE,resI,akk,xjac)
      elseif (iprob==3) then
       call Kaps(     programStep,probname,nveclen,temporal_splitting,uvec,ep,&
     &                uexact,dt,tfinal,iDT,resE,resI,akk,xjac)
      elseif (iprob==4) then
       call Kreiss(   programStep,probname,nveclen,temporal_splitting,uvec,ep,&
     &                uexact,dt,tfinal,iDT,resE,resI,akk,xjac)
      elseif (iprob==5) then !some sort of problem, fully implicit doesnt converge and the imex doesn't work
       call Lorenz(   programStep,probname,nveclen,temporal_splitting,uvec,ep,&
     &                uexact,dt,tfinal,iDT,resE,resI,akk,xjac)
      endif

      return
      end subroutine
