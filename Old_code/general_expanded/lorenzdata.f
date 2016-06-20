      program test_cases

c     program 1) van der Pol (Hairer II, pp 403)
c     program 2) Pureshi and Russo 
c     program 3) Dekker 7.5.2 pp. 215 (Kaps problem   : Index 1)
c     program 4) Dekker 7.5.1 pp. 214 (Kreiss' problem: Index 2)
c     WARNING: Kreiss has problem with algebraic variable
c     program 5) Lorenz

      parameter(is=9,ivarlen=4)
      parameter(isamp=71,jmax=81,jactual=81)

      integer cases 				!input range of runge kutta cases
      integer problem				!input problem number
      integer counter

      dimension aE(is,is),aI(is,is),bE(is),bI(is),cE(is),cI(is)
      dimension bEH(is),bIH(is),bD(is,4)
      dimension stageE(is),stageI(is),maxiter(is)
      dimension svpB(is,is,0:is),alpha(is,is)
      dimension al3N(is,is),al3D(is,is),bint(is)
      dimension al4N(is,is),al4D(is,is)

      dimension ustage(ivarlen,is),res(ivarlen,is)
      dimension predvec(ivarlen,is)

      dimension uvec(ivarlen),uveco(ivarlen),usum(ivarlen)
      dimension uveciter(ivarlen),uorig(ivarlen) 
      dimension uexact(ivarlen),Rnewton(ivarlen)
      dimension errvec(ivarlen),errvecT(ivarlen),tmpvec(ivarlen)
      dimension xjacinv(ivarlen,ivarlen)

      dimension cost(isamp),sig(isamp)

      dimension error1(isamp,ivarlen),error1P(isamp,ivarlen)

      dimension a(ivarlen*2),b(ivarlen*2),siga1(ivarlen*2)
      dimension sigb1(ivarlen*2),chi2(ivarlen*2)

      dimension b1save(jmax,ivarlen),b1Psave(jmax,ivarlen)
      dimension epsave(jmax)
      dimension ustor(100000000,ivarlen)
      
      counter=0
      ipred=2
      do icase = 1,1                                 !   begin algorithms loop

        icount = 0        !  cost counters
        jcount = 0        !  cost counters

        do i = 1,nrk
          stageE(i) = 0.0
          stageI(i) = 0.0
          maxiter(i)= 0
        enddo
      
        call rungeadd(aE,aI,bE,bI,cE,cI,nrk,bEH,bIH,icase,bD, !icase is the only input
     &   svpB(1,1,0),alpha,al3N,al3D,al4N,al4D,ipred)
 
        do iprob = 5,5                                !   begin problems loop
      
c         loop over different values of stiffness epsilon 

          do jepsil = 1,1,1                              !  begin stiffness epsilon loop

            ep = 1.                          !  used for 81 values of ep

            do iDT = 1,1,1                                 !  timestep loop for vdP, Kaps, etc

            call INIT(uvec,uexact,dt,iDT,tfinal,ep,nvecLen,iprob)     !  initialize problem information

            dto = 0.00001
            t = 0.

            do ival = 1,nvecLen
              errvecT(ival) = 0.0
            enddo

            do i = 1,nrk                                       !  initialize stage value preditor
              do ival = 1,nvecLen                              !  with trivial guess for first stage
                predvec(ival,i) = uvec(ival)
              enddo
            enddo


            do ktime = 1,100000000                      ! advance solution (time advancement loop )
              counter=counter+1
              if(t+dt.gt.tfinal)dt = tfinal-t + 1.0e-11

              do ival = 1,nvecLen
                uveco(ival) = uvec(ival)
              enddo

              jcount = jcount + (nrk-1)                         !  keep track of total RK stages 

              do L = 1,nrk                                      !  begin RK loop

                ttI = t + cI(L)*dt
                
                call RHS(uvec,res(1,L),dt,ep,ttI,iprob)

                do ival = 1,nvecLen                             !  write the solution into a storage register
                  ustage(ival,L) = uvec(ival)
                enddo

              if(L.ne.nrk)then

                ttI = t + cI(L+1)*dt

                do ival = 1,nvecLen
                  usum(ival) = uveco(ival)
                enddo
                do LL = 1,L 
                  do ival = 1,nvecLen
                    usum(ival) = usum(ival) + aI(L+1,LL)*res(ival,LL)
                  enddo
                enddo

                do ival = 1,nvecLen

                  if(ipred.eq.2)predvec(ival,L+1)=uvec(ival)    ! previous guess as starter
                  uvec(ival)  = predvec(ival,L+1)               ! put predict into start guess
                  uorig(ival) = uvec(ival)                      ! put predict into storage for testing
                enddo

                if(L.gt.1.and.L.lt.nrk)then
                  do ival = 1,nvecLen
                    uvec(ival)  = uveco(ival)               ! put predict into start guess
                    do jpre = 2,L
                       Z = ustage(ival,jpre)-uveco(ival)
                       uvec(ival) = uvec(ival) + alpha(L+1,jpre)*Z
                       uorig(ival) = uvec(ival)               ! put predict into storage for testing
                    enddo
                  enddo
                endif

c U^{(n+1,4)} = al4_{1}*U^{(n  ,4)} + al4_{2}*U^{(n  ,5)} +
c             + al4_{3}*U^{(n  ,6)} + al4_{4}*U^{(n+1,2)} +
c             + al4_{5}*U^{(n+1,3)}


c al4_{i} = \sum_{j=1}^{2*(order)} al4N(i,j)*r^{j} /
c           \sum_{j=1}^{2*(order)} al4D(i,j)*r^{j}

                if(L.eq.2.and.ktime.ne.1)then
                  the1 = 1.
                  the2 = the1*the
                  the3 = the2*the
                  the4 = the3*the
                  the5 = the4*the
                  the6 = the5*the
                  the7 = the6*the
                  the8 = the7*the
                  do LL = 1,5
                    xnum = al3N(LL,1)*the1
     &                   + al3N(LL,2)*the2
     &                   + al3N(LL,3)*the3
     &                   + al3N(LL,4)*the4
     &                   + al3N(LL,5)*the5
     &                   + al3N(LL,6)*the6
     &                   + al3N(LL,7)*the7
     &                   + al3N(LL,8)*the8
                    xden = al3D(LL,1)*the1
     &                   + al3D(LL,2)*the2
     &                   + al3D(LL,3)*the3
     &                   + al3D(LL,4)*the4
     &                   + al3D(LL,5)*the5
     &                   + al3D(LL,6)*the6
     &                   + al3D(LL,7)*the7
     &                   + al3D(LL,8)*the8
                    bint(LL) = xnum/xden
                  enddo
                    bint(1) =  9.9518675264213746 
                    bint(2) =  4.8366852488953721 
                    bint(3) =-24.163405114569394 
                    bint(4) = 14.152132944153401
                    bint(5) =  0.94399768676237158
                  do ival = 1,nvecLen
                    Z1 = ustage(ival,3)-ustage(ival,3)
                    Z2 = ustage(ival,4)-ustage(ival,3)
                    Z3 = ustage(ival,5)-ustage(ival,3)
                    Z4 = ustage(ival,6)-ustage(ival,3)
                    Z5 = ustage(ival,2)-ustage(ival,3)
                    uvec(ival) = ustage(ival,3) + bint(1)*z1
     &                                          + bint(2)*z2
     &                                          + bint(3)*z3
     &                                          + bint(4)*z4
     &                                          + bint(5)*z5
                    uorig(ival) = uvec(ival)               ! put predict into storage for testing
                  enddo

c                 do ival = 1,nvecLen
c                   uvec(ival) = bint(1)*ustage(ival,3)
c    &                         + bint(2)*ustage(ival,4)
c    &                         + bint(3)*ustage(ival,5)
c    &                         + bint(4)*ustage(ival,6)
c    &                         + bint(5)*ustage(ival,2)
c                   uvec(ival) = bint(1)*ustage(ival,4)
c    &                         + bint(2)*ustage(ival,5)
c    &                         + bint(3)*ustage(ival,6)
c    &                         + bint(4)*ustage(ival,2)
c    &                         + bint(5)*ustage(ival,3)
c                   uorig(ival) = uvec(ival)               ! put predict into storage for testing
c                 enddo
                endif

                do k = 1,20

                  icount = icount + 1

                  do ival = 1,nvecLen
                    uveciter(ival) = uvec(ival)
                  enddo

                  call RHS(uvec,res(1,L+1),dt,ep,ttI,iprob)

                  do ival = 1,nvecLen
                    Rnewton(ival) =  
     &              uvec(ival)  - aI(L+1,L+1)*res(ival,L+1) - usum(ival)
                  enddo

            call JACOB(uvec,xjacinv,dt,ep,aI(L+1,L+1),ttI,iprob,nvecLen)

                  do i = 1,nvecLen
                    do j = 1,nvecLen
                      uvec(i) = uvec(i) - xjacinv(i,j)*Rnewton(j)
                    enddo
                  enddo

                  tmp = 0.0
                  do ival = 1,nvecLen
                     tmp = tmp + abs(uvec(ival)-uveciter(ival))
                  enddo
                  if(tmp.lt.1.e-12) go to 160                 !  kick out of newton iteration

                enddo
  160           continue
 
             if((k.gt.maxiter(L+1)).and.(ktime.ne.1))maxiter(L+1)=k

              do ival = 1,nvecLen                    !  write the solution into a storage register
                ustage(ival,L+1) = uvec(ival)
              enddo

              xnorm = 0.0                                     !  assess how good initial guess was
              snorm = 0.0                                     !  assess how good initial guess was
              do ival = 1,nvecLen
                tmpD = (uvec(ival) - uorig(ival))
                xnorm = xnorm + tmpD*tmpD
                snorm = snorm + uvec(ival)*uvec(ival)
              enddo                                         
              xnorm = sqrt(xnorm/snorm)
              stageE(L+1) = stageE(L+1) + xnorm
              stageI(L+1) = stageI(L+1) + 1.*k

              elseif(L.eq.nrk) then
            
                do ival = 1,nvecLen
                  uvec(ival) = uveco(ival)
                enddo

                do ival = 1,nvecLen
                  do LL = 1,nrk 
                    uvec(ival) = uvec(ival) + bI(LL)*res(ival,LL)
                  enddo
                enddo

                if(time.lt.tfinal-1.0e-11)then               !  begin predicted error

                do ival = 1,nvecLen
                  errvec(ival) = 0.0
                  do LL = 1,nrk 
                    errvec(ival) = errvec(ival) 
     &                       + dt*( (bI(LL)-bIH(LL))*res(ival,LL) )
                  enddo
                  errvec(ival) = abs(errvec(ival))
                enddo
                rat = 1.0

                endif                                        !  end predicted error
              endif
	      do ii=1,nveclen
                ustor(ktime,ii)=uvec(ii)
              enddo
              print *,L  
            enddo                                            ! end RK loop
            do ival = 1,nvecLen
              errvecT(ival) = errvecT(ival) + errvec(ival)
            enddo                                         
            t = t + dt
            print *,ktime
            if(t.ge.tfinal) go to 100
            
          enddo                                              ! end time advancement loop
  100  continue  
         enddo                                             !  end  loop over different dt
         enddo                                              !  end stiffness epsilon loop
        enddo                                            !  end problems loop
       enddo                                             !  end algorithms loop
      !print *, ustor(:,:)
      stop
      end


    
