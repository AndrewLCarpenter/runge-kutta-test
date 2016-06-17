      SUBROUTINE gser(gamser,a,x,gln)

      implicit none
      
      integer,  parameter        ::  wp=8

      integer,    parameter      ::  itkax=100
      real(wp),   parameter      ::  eps=3.e-7

      real(wp),   intent(  out)  ::  gamser
      real(wp),   intent(in   )  ::  a,x
      real(wp),   intent(  out)  ::  gln

      integer                    ::  n
      real(wp)                   ::  ap,del,summ,gammln

!     USES gammln
!     Returns the incomplete gamma function P(a, x) evaluated by its series 
!     representation as gamser. Also returns In r( a) as gln.

      gln=gammln (A)
      if(x.le.0.)then
      if(x.lt.0.)pause 'x < 0 in gser'
      gamser=0.
      return
      endif
      ap=a
      summ=1./a
      del=summ
      do n=1,itkax
      ap=ap+1.
      del=del*x/ap
      summ=summ+del
      if(abs(del).lt.abs(summ)*eps)goto 1
      enddo 
      pause 'a too large, itkax too small in gser'
    1 gamser=summ*exp(-x+a*log(x)-gln)
      return
      end
