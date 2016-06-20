c--------------------------------------------------------------
      subroutine JACOB(uvec,xinv,dt,ep,akk,ttI,iprob,nvecLen)
      parameter(is=9,ivarlen=4)
      dimension uvec(ivarlen),xjac(ivarlen,ivarlen)
      dimension xinv(ivarlen,ivarlen)
      parameter(sigma=10.,rho=28.,beta=8./3.) !Lorenz constants for Prob. 5

   
      if(iprob.eq.1)then
        xjac(1,1) = 1.-akk*dt*(0.)
        xjac(1,2) = 0.-akk*dt*(1.)
        xjac(2,1) = 0.-akk*dt*(-2*uvec(1)*uvec(2)-1)/ep
        xjac(2,2) = 1.-akk*dt*(1-uvec(1)*uvec(1))/ep
      elseif(iprob.eq.2)then
        xjac(1,1) = 1.-akk*dt*(0.)
        xjac(1,2) = 0.-akk*dt*(0.)
        xjac(2,1) = 0.-akk*dt*(cos(uvec(1)) )/ep
        xjac(2,2) = 1.-akk*dt*(-1.)/ep
      elseif(iprob.eq.3)then
        xjac(1,1) = 1.-akk*dt*(-(1./ep))
        xjac(1,2) = 0.-akk*dt*(+1./ep)*2*uvec(2)
        xjac(2,1) = 0.-akk*dt*(0)
        xjac(2,2) = 1.-akk*dt*(0)
      elseif(iprob.eq.4)then
        xjac(1,1) = 1.-akk*dt
        xjac(1,2) = 0.-akk*dt/ep
        xjac(2,1) = 0.-akk*dt*(0)
        xjac(2,2) = 1.+akk*dt/ep !possible error here
       elseif(iprob.eq.5)then
        xjac(1,1) = 1.-akk*dt*(-sigma)/ep
        xjac(1,2) = 0.-akk*dt*(sigma)/ep
        xjac(1,3) = 0.-akk*dt*(0)/ep

        xjac(2,1) = 0.-akk*dt*(rho-uvec(3))
        xjac(2,2) = 1.-akk*dt*(-1.)
        xjac(2,3) = 0.-akk*dt*(-uvec(1))

        xjac(3,1) = 0.-akk*dt*(uvec(2))
        xjac(3,2) = 0.-akk*dt*(uvec(1))
        xjac(3,3) = 1.-akk*dt*(-beta)
      endif

      if(nvecLen.eq.2)then
        det = (xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1))
        xinv(1,1) =  xjac(2,2)/det
        xinv(1,2) = -xjac(1,2)/det
        xinv(2,1) = -xjac(2,1)/det
        xinv(2,2) =  xjac(1,1)/det

      elseif(nvecLen.eq.3)then

        x11 = xjac(1,1)
        x12 = xjac(1,2)
        x13 = xjac(1,3)
        x21 = xjac(2,1)
        x22 = xjac(2,2)
        x23 = xjac(2,3)
        x31 = xjac(3,1)
        x32 = xjac(3,2)
        x33 = xjac(3,3)

        det = - x13*x22*x31 + x12*x23*x31 +  x13*x21*x32 
     &        - x11*x23*x32 - x12*x21*x33 +  x11*x22*x33

        detI = 1./det

        xinv(1,1) =  (- x23*x32 + x22*x33) * detI
        xinv(1,2) =  (+ x13*x32 - x12*x33) * detI
        xinv(1,3) =  (- x13*x22 + x12*x23) * detI
        xinv(2,1) =  (+ x23*x31 - x21*x33) * detI
        xinv(2,2) =  (- x13*x31 + x11*x33) * detI
        xinv(2,3) =  (+ x13*x21 - x11*x23) * detI
        xinv(3,1) =  (- x22*x31 + x21*x32) * detI
        xinv(3,2) =  (+ x12*x31 - x11*x32) * detI
        xinv(3,3) =  (- x12*x21 + x11*x22) * detI

      elseif(nvecLen.eq.4)then

        x11 = xjac(1,1)
        x12 = xjac(1,2)
        x13 = xjac(1,3)
        x14 = xjac(1,4)
        x21 = xjac(2,1)
        x22 = xjac(2,2)
        x23 = xjac(2,3)
        x24 = xjac(2,4)
        x31 = xjac(3,1)
        x32 = xjac(3,2)
        x33 = xjac(3,3)
        x34 = xjac(3,4)
        x41 = xjac(4,1)
        x42 = xjac(4,2)
        x43 = xjac(4,3)
        x44 = xjac(4,4)

        det = 
     -  (x14*x23*x32*x41 - x13*x24*x32*x41 - 
     -   x14*x22*x33*x41 + x12*x24*x33*x41 + 
     -   x13*x22*x34*x41 - x12*x23*x34*x41 - 
     -   x14*x23*x31*x42 + x13*x24*x31*x42 + 
     -   x14*x21*x33*x42 - x11*x24*x33*x42 - 
     -   x13*x21*x34*x42 + x11*x23*x34*x42 + 
     -   x14*x22*x31*x43 - x12*x24*x31*x43 - 
     -   x14*x21*x32*x43 + x11*x24*x32*x43 + 
     -   x12*x21*x34*x43 - x11*x22*x34*x43 - 
     -   x13*x22*x31*x44 + x12*x23*x31*x44 + 
     -   x13*x21*x32*x44 - x11*x23*x32*x44 - 
     -   x12*x21*x33*x44 + x11*x22*x33*x44)

        detI = 1./det

        xinv(1,1) = (
     - -(x24*x33*x42) + x23*x34*x42 + x24*x32*x43 - 
     -   x22*x34*x43  - x23*x32*x44 + x22*x33*x44  ) * detI
        xinv(1,2) = (
     -   x14*x33*x42  - x13*x34*x42 - x14*x32*x43 + 
     -   x12*x34*x43  + x13*x32*x44 - x12*x33*x44  ) * detI
        xinv(1,3) = (
     - -(x14*x23*x42) + x13*x24*x42 + x14*x22*x43 - 
     -   x12*x24*x43  - x13*x22*x44 + x12*x23*x44  ) * detI
        xinv(1,4) = (
     -   x14*x23*x32  - x13*x24*x32 - x14*x22*x33 + 
     -   x12*x24*x33  + x13*x22*x34 - x12*x23*x34  ) * detI
        xinv(2,1) = (
     -   x24*x33*x41  - x23*x34*x41 - x24*x31*x43 + 
     -   x21*x34*x43  + x23*x31*x44 - x21*x33*x44  ) * detI
        xinv(2,2) = (
     - -(x14*x33*x41) + x13*x34*x41 + x14*x31*x43 - 
     -   x11*x34*x43  - x13*x31*x44 + x11*x33*x44  ) * detI
        xinv(2,3) = (
     -   x14*x23*x41  - x13*x24*x41 - x14*x21*x43 + 
     -   x11*x24*x43  + x13*x21*x44 - x11*x23*x44  ) * detI
        xinv(2,4) = (
     - -(x14*x23*x31) + x13*x24*x31 + x14*x21*x33 - 
     -   x11*x24*x33  - x13*x21*x34 + x11*x23*x34  ) * detI
        xinv(3,1) = (
     - -(x24*x32*x41) + x22*x34*x41 + x24*x31*x42 - 
     -   x21*x34*x42  - x22*x31*x44 + x21*x32*x44  ) * detI
        xinv(3,2) = (
     -   x14*x32*x41  - x12*x34*x41 - x14*x31*x42 + 
     -   x11*x34*x42  + x12*x31*x44 - x11*x32*x44  ) * detI
        xinv(3,3) = (
     - -(x14*x22*x41) + x12*x24*x41 + x14*x21*x42 - 
     -   x11*x24*x42  - x12*x21*x44 + x11*x22*x44  ) * detI
        xinv(3,4) = (
     -   x14*x22*x31  - x12*x24*x31 - x14*x21*x32 + 
     -   x11*x24*x32  + x12*x21*x34 - x11*x22*x34  ) * detI
        xinv(4,1) = (
     -   x23*x32*x41  - x22*x33*x41 - x23*x31*x42 + 
     -   x21*x33*x42  + x22*x31*x43 - x21*x32*x43  ) * detI
        xinv(4,2) = (
     - -(x13*x32*x41) + x12*x33*x41 + x13*x31*x42 - 
     -   x11*x33*x42  - x12*x31*x43 + x11*x32*x43  ) * detI
        xinv(4,3) = (
     -   x13*x22*x41  - x12*x23*x41 - x13*x21*x42 + 
     -   x11*x23*x42  + x12*x21*x43 - x11*x22*x43  ) * detI
        xinv(4,4) = (
     - -(x13*x22*x31) + x12*x23*x31 + x13*x21*x32 - 
     -   x11*x23*x32  - x12*x21*x33 + x11*x22*x33  ) * detI

      endif

  
      return
      end
