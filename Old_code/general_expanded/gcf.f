      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER itkax
      REAL a,gammcf,gln,x,eps,FPMIN
      PARAMETER (itkax=100,eps=3.e-7,FPMIN=1.e-30)
c     USES gammln
      INTEGER i
      REAL an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1. -a 
      c=1. /FPMIN 
      d=1./b
      h=d
      do i=1, itkax 
      an=-i*(i-a)
      b=b+2.
      d=an*d+b
      if(abs(d).lt.FPMIN)d=FPMIN
      c=b+an/c
      if(abs(c).lt.FPMIN)c=FPMIN
      d=1./d
      del=d*c
      h=h*del
      if(abs(del-1.).lt.EPS)goto 1
      enddo
      pause 'a too large, ITMAX too small in gcf'
    1 gammcf=exp(-x+a*log(x)-gln)*h
      return
      end
