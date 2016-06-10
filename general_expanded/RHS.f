c-------------------------------------------------------------
      subroutine RHS(uvec,res,dt,ep,ttI,iprob)
      parameter(is=9,ivarlen=4)
      parameter(sigma=10.,rho=28.,beta=8./3.) !Lorenz constants for Prob. 5
      dimension uvec(ivarlen),res(ivarlen)
         
      if    (iprob.eq.1)then
        res(1) = dt*uvec(2)
        res(2) = dt*((1-uvec(1)*uvec(1))*uvec(2) - uvec(1))/ep
      elseif(iprob.eq.2)then
        res(1) = dt*(-uvec(2))
        res(2) = dt*( uvec(1) + (sin(uvec(1)) - uvec(2))/ep)
      elseif(iprob.eq.3)then
        res(1) = dt*(-(1./ep+2.)*uvec(1) + (1./ep)*uvec(2)*uvec(2))
        res(2) = dt*(uvec(1) - uvec(2) - uvec(2)*uvec(2) )
      elseif(iprob.eq.4)then
        res(1) = dt*(uvec(1) + uvec(2)/ep)
        res(2) = dt*(        - uvec(2))
      elseif(iprob.eq.5)then
        res(1)=dt*sigma*(uvec(2)-uvec(1))/ep
        res(2)=dt*(-uvec(1)*uvec(3)+rho*uvec(1)-uvec(2))
        res(3)=dt*(uvec(1)*uvec(2)-beta*uvec(3))
      endif

      return
      end
