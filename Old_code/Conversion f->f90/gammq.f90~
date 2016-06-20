      FUNCTION gammq(a,x)

      implicit none
      
      integer,  parameter        ::  wp=8

      real(wp),   intent(in   )  ::  a,x
      real(wp)     ::  gammq

      real(wp)                  :: gamser,gln,gammcf

!     USES gcf,gser
!     Returns the incomplete gamma function Q(a, x) = 1 -P(a, x).

      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammq'
      if (x.lt.a+1.) then 
      call gser(gamser,a,x,gln)
      gammq=1. -gamser 
      call gcf(gammcf,a,x,gln)
      endif
      return
      end
