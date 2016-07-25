!******************************************************************************
! Module to take data from test_cases and perform linear regression
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
!******************************************************************************

      module  poly_fit_Mod

      use precision_vars     , only : wp

      implicit none

      private
      public  :: fit

      contains

!==============================================================================
!  PERFORM LINEAR REGRESSION USING ROUTINES AND FUNCTIONS CONTAINED BELOW
      subroutine fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
      
      integer,                      intent(in   )  ::  mwt,ndata
      real(wp),                     intent(  out)  ::  a,b,siga,sigb,chi2,q
      real(wp),   dimension(ndata),  intent(in   ) ::  sig,x,y

      real(wp)                    ::  sigdat,ss,st2,sx,sxoss,sy,t,wt
      integer                     :: i

      sx=0.0_wp
      sy=0.0_wp
      st2=0.0_wp
      b=0.0_wp
      if(mwt /= 0) then 
        ss=0.0_wp
        do i=1,ndata
          wt=1.0_wp/(sig(i)**2)
          ss=ss+wt
          sx=sx+x(i)*wt
          sy=sy+y(i)*wt
        enddo
      else
        do i=1,ndata
          sx=sx+x(i)
          sy=sy+y(i)
        enddo 
        ss=float(ndata)
      endif
      sxoss=sx/ss
      if(mwt /= 0) then
        do i=1,ndata
          t=(x(i)-sxoss)/sig(i)
          st2=st2+t*t
          b=b+t*y(i)/sig(i)
        enddo
      else
        do i=1,ndata
          t=x(i)-sxoss
          st2=st2+t*t
          b=b+t*y(i)
        enddo 
      endif
      b=b/st2
      a=(sy-sx*b)/ss
      siga=sqrt((1.0_wp+sx*sx/(ss*st2))/ss)
      sigb=sqrt(1.0_wp/st2)
      chi2=0.0_wp
      if(mwt == 0) then
        do i=1,ndata
          chi2=chi2+(y(i)-a-b*x(i))**2
        enddo 
        q=1.0_wp
        sigdat=sqrt(chi2/(ndata-2) )
        siga=siga*sigdat
        sigb=sigb*sigdat
      else
        do  i=1,ndata
          chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
        enddo 
        q=gammq(0.5_wp*(ndata-2),0.5_wp*chi2)
      endif
      return
      end subroutine fit

! ==============================================================================

      function gammln(xx)

      real(wp),   intent(in   )  ::  xx
      real(wp)                   ::  gammln

      integer                    ::  j
      real(wp)                   ::  ser,stp,tmp,x,y
      real(wp),   dimension(6)   ::  cof

      SAVE cof,stp
      DATA cof,stp/76.18009172947146_wp,-86.50532032941677_wp,&
     &24.01409824083091_wp,-1.231739572450155_wp,.1208650973866179d-2,&
     &-.5395239384953d-5,2.5066282746310005_wp/
      x=xx
      y=x
      tmp=x+5.5_wp
      tmp=(x+0.5_wp)*log(tmp)-tmp
      ser=1.000000000190015_wp
      do j=1,6
      y=y+1._wp
      ser=ser+cof(j)/y
      enddo 
      gammln=tmp+log(stp*ser/x)
      return
      end function gammln

! ==============================================================================

      function gammq(a,x)

      real(wp),   intent(in   )  ::  a,x
      real(wp)                   ::  gammq

      real(wp)                   :: gamser,gln,gammcf
 
!     USES gcf,gser
!     Returns the incomplete gamma function Q(a, x) = 1 -P(a, x).
      gammq = 0.0_wp

      if(x <  0.0_wp .or. a <= 0.0_wp)then
        write(*,*)'bad arguments in gammq'
      endif
      if (x <  a+1.0_wp) then 
      call gser(gamser,a,x,gln)
      gammq=1.0_wp -gamser 
      call gcf(gammcf,a,x,gln)
      endif
      return
      end function gammq
 
! ==============================================================================

      subroutine gcf(gammcf,a,x,gln)

      integer,    parameter      ::  itkax=100
      real(wp),   parameter      ::  eps=3.e-7
      real(wp),   parameter      ::  FPMIN=1e-30

      real(wp),   intent(  out)  ::  gammcf
      real(wp),   intent(in   )  ::  a,x
      real(wp),   intent(  out)  ::  gln

      integer                    ::  i
      real(wp)                   ::  an,b,c,d,del,h

      gln=gammln(a)
      b=x+1.0_wp -a 
      c=1.0_wp /FPMIN 
      d=1.0_wp /b
      h=d
      do i=1, itkax 
      an=-i*(i-a)
      b=b+2.0_wp
      d=an*d+b
      if(abs(d) <  FPMIN)d=FPMIN
      c=b+an/c
      if(abs(c) <  FPMIN)c=FPMIN
      d=1.0_wp/d
      del=d*c
      h=h*del
      if(abs(del-1.0_wp) <  EPS)goto 1
      enddo
      write(*,*)'a too large, ITMAX too small in gcf'
    1 gammcf=exp(-x+a*log(x)-gln)*h
      return
      end subroutine gcf

! ==============================================================================

      subroutine gser(gamser,a,x,gln)

      integer,    parameter      ::  itkax=100
      real(wp),   parameter      ::  eps=3.e-7

      real(wp),   intent(  out)  ::  gamser
      real(wp),   intent(in   )  ::  a,x
      real(wp),   intent(  out)  ::  gln

      integer                    ::  n
      real(wp)                   ::  ap,del,summ

!     Returns the incomplete gamma function P(a, x) evaluated by its series 
!     representation as gamser. Also returns In r( a) as gln.

      gln=gammln (A)
      if(x <= 0.0_wp)then
      if(x <  0.0_wp)write(*,*)'x < 0 in gser'
      gamser=0.0_wp
      return
      endif
      ap=a
      summ=1.0_wp/a
      del=summ
      do n=1,itkax
      ap=ap+1.0_wp
      del=del*x/ap
      summ=summ+del
      if(abs(del) <  abs(summ)*eps)goto 1
      enddo 
      write(*,*)'a too large, itkax too small in gser'
    1 gamser=summ*exp(-x+a*log(x)-gln)
      return
      end subroutine gser

! ==============================================================================


      end module  poly_fit_Mod
