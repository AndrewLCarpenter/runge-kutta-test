      FUNCTION gammq(a,x)
      REAL a,gammq,x
c     USES gcf,gser
c     Returns the incomplete gamma function Q(a, x) = 1 -P(a, x).
      REAL gammcf ,gamser ,gln
      if(x.lt.0. .or.a.le.0.)pause 'bad arguments in gammq'
      if (x .lt .a+1. ) then 
      call gser(gamser,a,x,gln)
      gammq=1. -gamser 
      call gcf(gammcf,a,x,gln)
      endif
      return
      end
