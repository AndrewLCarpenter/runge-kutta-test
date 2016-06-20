      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER itkax
      REAL a,gamser,gln,x,eps
      PARAMETER (itkax=100,eps=3.e-7)
c     USES gammln
c     Returns the incomplete gamma function P(a, x) evaluated by its series 
c     representation as gamser. Also returns In r( a) as gln.
      INTEGER n
      REAL ap,del,sum,gammln
      gln=gammln (A)
      if(x.le.0.)then
      if(x.lt.0.)pause 'x < 0 in gser'
      gamser=0.
      return
      endif
      ap=a
      sum=1./a
      del=sum
      do n=1,itkax
      ap=ap+1.
      del=del*x/ap
      sum=sum+del
      if(abs(del).lt.abs(sum)*eps)goto 1
      enddo 
      pause 'a too large, itkax too small in gser'
    1 gamser=sum*exp(-x+a*log(x)-gln)
      return
      end
