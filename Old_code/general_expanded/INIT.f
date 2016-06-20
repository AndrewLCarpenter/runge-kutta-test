      subroutine INIT(uvec,uexact,dt,iDT,tfinal,ep,nvecLen,iprob)
      parameter(is=9,ivarlen=4)
        dimension uvec(ivarlen),uexact(ivarlen)
        dimension ExactTot(81,ivarlen)

        if(iprob.eq.1)then                 !  van der Pol (Hairer II, pp 403)
c    vanderPol has 8 decades with 10 levels of epsilon per decade plus 1st pt of 9th decade
 
         open(unit=39,file='exact.vanderpol.data')

          rewind(39)
          do i=1,81
            read(39,*)ExactTot(i,1),ExactTot(i,2)
            ExactTot(i,3) = 1./10**((i-1)/(10.))                  !  used for 81 values of ep
          enddo

          do i=1,81
            diff = abs(ExactTot(i,3) - ep)
            if(diff.le.1.0e-10)then
              uexact(1) = ExactTot(i,1)
              uexact(2) = ExactTot(i,2)
              go to 1 
            endif
          enddo
   1      continue

          dt = 0.5/10**((iDT-1.)/20.)
          nvecLen = 2
          tfinal = 0.5

          uvec(1) = 2.
          uvec(2) = -0.6666654321121172

        elseif(iprob.eq.2)then                     !  Pureshi and Russo 
 
          open(unit=39,file='exact.pureschi.1.data')
c         open(unit=39,file='exact.pureschi.2.data')


          rewind(39)
          do i=1,81
            read(39,*)ExactTot(i,1),ExactTot(i,2)
            ExactTot(i,3) = 1./10**((i-1)/(10.))                  !  used for 81 values of ep
          enddo

          do i=1,81
            diff = abs(ExactTot(i,3) - ep)
            if(diff.le.1.0e-10)then
              uexact(1) = ExactTot(i,1)
              uexact(2) = ExactTot(i,2)
              go to 2
            endif
          enddo
   2      continue

          dt = 0.5/10**((iDT-1.)/20.)
          nvecLen = 2
          tfinal = 5.0

          pi = acos(-1.)

c  IC: problem 1   :  equilibrium IC
          uvec(1) = pi/2.
          uvec(2) = sin(pi/2.)
c  IC: problem 2   :  non-equilibrium IC
c         uvec(1) = pi/2.
c         uvec(2) = 1./2.

        elseif(iprob.eq.3)then                   ! Dekker 7.5.2 pp. 215 (Kaps problem   : Index 1)
 
          dt = 0.5/10**((iDT-1.)/20.)
          nvecLen = 2

          tfinal = 1.0
          tmp = exp(-tfinal)
          uexact(1) = tmp*tmp
          uexact(2) = tmp

c  IC: problem 1   :  equilibrium IC
          uvec(1) = 1.
          uvec(2) = 1.
c  IC: problem 2   :  non-equilibrium IC
c         uvec(1) = ?
c         uvec(2) = ?

        elseif(iprob.eq.4)then          ! Dekker 7.5.1 pp. 214 (Kreiss' problem: Index 2)
 
          dt = 0.25/10**((iDT-1.)/20.)
          nvecLen = 2

          tfinal = 1.0
          tmpP = exp(+tfinal   )
          tmpM = exp(-tfinal/ep)
          tmp = exp(tfinal)
          uo = 1.
          vo = 1.
          uexact(1) = uo*tmpP + vo*(tmpP - tmpM)/(1. + ep)
          uexact(2) = vo*tmpM

c  IC: problem 1   :  equilibrium IC
          uvec(1) = uo
          uvec(2) = vo
c  IC: problem 2   :  non-equilibrium IC
c         uvec(1) = ?
c         uvec(2) = ?

	elseif(iprob.eq.5)then		!Lorenz
	
	  dt = 0.5/10**((iDT-1.)/20.)
	  nveclen = 3 

	  tfinal = 1.0 !!arbitrary

	  open(unit=39,file='exact.lorenz.data')		! need file still

          rewind(39)
          do i=1,81
            read(39,*)ExactTot(i,1),ExactTot(i,2),ExactTot(i,3)
            ExactTot(i,4) = 1./10**((i-1)/(10.))                  !  used for 81 values of ep
          enddo

          do i=1,81
            diff = abs(ExactTot(i,4) - ep)
            if(diff.le.1.0e-10)then
              uexact(1) = ExactTot(i,1)
              uexact(2) = ExactTot(i,2)
              uexact(3) = ExactTot(i,3)
              go to 5 
            endif
          enddo
   5      continue

c  IC: problem 1 : Lorenz(1963)
	  uvec(1) = 0
          uvec(2) = 1
          uvec(3) = 0

        endif

      return
      end
