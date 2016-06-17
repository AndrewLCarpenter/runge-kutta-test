!####################################################################
!### Program to test Runge-Kutta IMEX schemes on various problems ###
!####################################################################
!
!----------------------------PROBLEM LIST----------------------------
!     problem 1) van der Pol (Hairer II, pp 403)
!     problem 2) Pureshi and Russo 
!     problem 3) Dekker 7.5.2 pp. 215 (Kaps problem   : Index 1)
!     problem 4) Dekker 7.5.1 pp. 214 (Kreiss' problem: Index 2)
!     WARNING: Kreiss has problem with algebraic variable
!     problem 5) Lorenz
!     problem 6) Burger's
!     problem 7) Black-Scholes
!--------------------------------------------------------------------
!
!-----------------------------REQUIRED FILES-------------------------
! RUNGEADD.F90  *CONTAINS RK CONSTANTS
! INIT.F90      *CONTAINS INITIALIZATIONS FOR PROBLEMS
! RHS.F90       *CONTAINS THE RIGHT-HAND SIDE OF THE PROBLEMS
! JACOB.F90     *COMPUTES THE JACOBIAN FOR EACH PROBLEM
! FIT.F90       *STATISTICS?
! GSER.F90      *STATISTICS?
! GCF.F90       *STATISTICS?
! GAMMQ.F90     *STATISTICS?
! GAMMLN.F90    *STATISTICS?
!--------------------------------------------------------------------
!
!--------------------------HOW TO USE--------------------------------
!**Adding test cases:
!      RHS,JACOB,& INIT need to be updated with information about each
!      test case. Also update probname array with problem name (9 char)
!**Compiling
!      Compile all files listed above together in the same directory
!**Running test cases
!      Program will prompt user for which predictor, case, and problem
!      to run.
!**Output results
!      The program outputs results to files with the form:
! /"probname"/"casename"/"probname"_"casename"_"variable number".dat
! /"probname"/"casename"/"probname"_"casename"_"variable number"P.dat
! /"probname"/"casename"/"probname"_"casename"_conv.dat
!
!***************************BEGIN PROGRAM****************************
      program test_cases

      use precision_vars
      use CSR_Variables

      implicit none     
!------------------------------PARAMETERS----------------------------
!     integer, parameter :: wp=8               !working precision
      integer, parameter :: is=9               !max constant length
      integer, parameter :: ivarlen=4          !max variable length
      integer, parameter :: isamp=71           !?
      integer, parameter :: jmax=81            !?
      integer, parameter :: jactual=81         !?
!-----------------------------VARIABLES------------------------------      
      integer             :: ipred,cases,problem         !user inputs
      integer             :: icase,icount,jcount,i,ii    !do loops
      integer             :: iprob,nrk,jepsil,iDT,ival   !do loops
      integer             :: ktime,L,LL,nveclen          !do loops
      integer             :: jpre,inew,j,k,jsamp,istage  !do loops

      real(wp)            :: dt,itmp,ep,t,sigma,rho,beta,tfinal,ttE
      real(wp)            :: ttI,Z,time,rat,bb,xnum,xden,z1,z2,z3,z4,z5
      real(wp)            :: the,the1,the2,the3,the4,the5,the6,the7,the8,q
      real(wp)            :: tmp,xnorm,snorm,tmpD,dto,totalerror,totalerrorp

      real(wp), dimension(is,is) :: aE,aI,alpha,al3N,al3D,al4N,al4D
      real(wp), dimension(is) :: bE,bI,cE,cI,bEH,bIH,stageE,stageI
      real(wp), dimension(is) :: maxiter,bint
      real(wp), dimension(is,ivarlen) :: bD
      real(wp), dimension(is,is,0:is) :: svpB
      real(wp), dimension(ivarlen,is) :: ustage,resE,resI,predvec
      real(wp), dimension(ivarlen) :: uvec,uveco,usum,uveciter,uorig
      real(wp), dimension(ivarlen) :: uexact,Rnewton,errvec,errvecT,tmpvec
      real(wp), dimension(ivarlen,ivarlen) :: xjac,xjacinv
      real(wp), dimension(isamp) :: cost,sig
      real(wp), dimension(isamp,ivarlen) :: error1,error1P
      real(wp), dimension(ivarlen*2) :: a,b,siga1,sigb1,chi2
      real(wp), dimension(jmax,ivarlen) :: b1save,b1Psave
      real(wp), dimension(jmax) :: epsave
      
      character(len=80) :: filename,fileloc   !file name and location
      character(len=25) :: casename           !name of RK case
      character(len=9), dimension(7) :: probname!array of problem names
      character(len=4)  :: ext='.dat'         !file extension
      character         :: istr               !loop index placeholder

      logical           :: dirExists          !check if directory exists
     
      !**ENTER PROBLEM NAMES**
      probname(1)="vanderPol"
      probname(2)="Pureshi  "
      probname(3)="Kaps     "
      probname(4)="Kreiss   "
      probname(5)="Lorenz   "
      probname(6)="Burgers  "
      probname(7)="BlackScho"      
!-----------------------------USER INPUT-----------------------------
      write(*,*)'what is ipred?' !input predictor number
      read(*,*)ipred
      write(*,*)'what is case?'  !input range of runge kutta cases
      read(*,*)cases
      write(*,*)'which problem?' !input problem number
      read(*,*)problem
!-------------------------ALGORITHMS LOOP----------------------------
      do icase = 7,18 

        icount = 0                                  !cost counters
        jcount = 0                                  !cost counters
        
        !**GET RK COEFFICIENTS**
        call rungeadd(aE,aI,bE,bI,cE,cI,nrk,bEH,bIH,icase,bD, &   
     &   svpB(1,1,0),alpha,al3N,al3D,al4N,al4D,casename)

        !**initilizations?**
        do i = 1,nrk
          stageE(i) = 0.0
          stageI(i) = 0.0
          maxiter(i)= 0
        enddo
!--------------------------PROBLEMS LOOP------------------------------
        do iprob = problem,problem
          !**INIT. FILE PATH FOR OUTPUTS**                           
          fileloc=trim(probname(iprob))//'/'//trim(casename)//'/'
          inquire( file=trim(fileloc)//'/.', exist=dirExists )  
          if (.not.dirExists)call system('mkdir '//trim(fileloc)) ! create directory

          call Allocate_CSR_Storage(problem, nvecLen)

!--------------------------STIFFNESS LOOP-----------------------------
          do jepsil = 1,jactual,1                              

            itmp = 11 - jmax/jactual         !used for 81 values of ep
            ep = 1.0_wp/10**((jepsil-1)/(itmp*1.0_wp))           
            
            !**INIT. OUTPUT FILES**
            do ii = 1,ivarlen
              write(istr,"(I1.1)")ii
              filename=trim(probname(iprob))//'_'//trim(casename)//'_'//istr//ext
              open(49+ii,file=trim(fileloc)//filename)
              write(49+ii,*)'zone T = "ep = ',ep,'",'
              filename=trim(probname(iprob))//'_'//trim(casename)//'_'//istr//'P'//ext
              open(59+ii,file=trim(fileloc)//filename)
              write(59+ii,*)'zone T = "ep = ',ep,'",'
            enddo
!--------------------------TIMESTEP LOOP------------------------------
            do iDT = 1,isamp,1                     

            !**INITIALIZE PROBLEM INFORMATION**
            call INIT(uvec,uexact,dt,iDT,tfinal,ep,nvecLen,iprob,sigma,rho,beta)     
            dto = dt        !store time step
            t = 0.0_wp      !init. start time

            !**INIT. ERROR VECTOR**
            do ival = 1,nvecLen
              errvecT(ival) = 0.0_wp
            enddo

            do i = 1,nrk            !initialize stage value preditor
              do ival = 1,nvecLen   !with trivial guess for first stage
                predvec(ival,i) = uvec(ival)
              enddo
            enddo
!--------------------------TIME ADVANCEMENT LOOP----------------------
            do ktime = 1,100000000                      
              if(t+dt.gt.tfinal)dt = tfinal-t + 1.0d-11 !check if dt>tfinal
              
              !**STORE VALUES OF UVEC**
              do ival = 1,nvecLen
                uveco(ival) = uvec(ival)
              enddo
              jcount = jcount + (nrk-1)    !keep track of total RK stages 
!--------------------------------RK LOOP------------------------------
              do L = 1,nrk    
                ttI = t + cI(L)*dt
                ttE = t + cE(L)*dt

                !**GET INFORMATION ABOUT RHS OF PROBLEM**
                call RHS(uvec,resI(1,L),resE(1,L),dt,ep,iprob,sigma,rho,beta)

                do ival = 1,nvecLen !write the solution into a storage register
                  ustage(ival,L) = uvec(ival)
                enddo

              if(L.ne.nrk)then !True for all but last loop
                ttI = t + cI(L+1)*dt   
                ttE = t + cE(L+1)*dt

                do ival = 1,nvecLen !start the summation vector
                  usum(ival) = uveco(ival)
                enddo

                do LL = 1,L 
                  do ival = 1,nvecLen
                    usum(ival) = usum(ival) + aI(L+1,LL)*resI(ival,LL)+ &
     &                           aE(L+1,LL)*resE(ival,LL) !imex summation
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

! U^{(n+1,4)} = al4_{1}*U^{(n  ,4)} + al4_{2}*U^{(n  ,5)} +
!             + al4_{3}*U^{(n  ,6)} + al4_{4}*U^{(n+1,2)} +
!             + al4_{5}*U^{(n+1,3)}


! al4_{i} = \sum_{j=1}^{2*(order)} al4N(i,j)*r^{j} /
!           \sum_{j=1}^{2*(order)} al4D(i,j)*r^{j}

                if(L.eq.2.and.ktime.ne.1)then
!                  the1 = 1.0_wp             !"the" is undefined
!                  the2 = the1*the           !bint gets set and then overwritten  
!***********************************************************************
!                  the3 = the2*the           !unless nrk> 5 then bint gets set for L[1,nrk], LL=L
!                  the4 = the3*the
!                  the5 = the4*the
!                  the6 = the5*the
!                  the7 = the6*the
!                  the8 = the7*the
!                  do LL = 1,5
!                    xnum = al3N(LL,1)*the1&
!     &                   + al3N(LL,2)*the2&
!     &                   + al3N(LL,3)*the3&
!     &                   + al3N(LL,4)*the4&
!     &                   + al3N(LL,5)*the5&
!     &                   + al3N(LL,6)*the6&
!     &                   + al3N(LL,7)*the7&
!     &                   + al3N(LL,8)*the8
!                    xden = al3D(LL,1)*the1&
!     &                   + al3D(LL,2)*the2&
!     &                   + al3D(LL,3)*the3&
!     &                   + al3D(LL,4)*the4&
!     &                   + al3D(LL,5)*the5&
!     &                   + al3D(LL,6)*the6&
!     &                   + al3D(LL,7)*the7&
!     &                   + al3D(LL,8)*the8
!                    bint(LL) = xnum/xden
!                  enddo
                    bint(1) =  9.9518675264213746_wp 
                    bint(2) =  4.8366852488953721_wp 
                    bint(3) =-24.163405114569394_wp 
                    bint(4) = 14.152132944153401_wp
                    bint(5) =  0.94399768676237158_wp
                  do ival = 1,nvecLen
                    Z1 = ustage(ival,3)-ustage(ival,3)
                    Z2 = ustage(ival,4)-ustage(ival,3)
                    Z3 = ustage(ival,5)-ustage(ival,3)
                    Z4 = ustage(ival,6)-ustage(ival,3)
                    Z5 = ustage(ival,2)-ustage(ival,3)
                    uvec(ival) = ustage(ival,3) + bint(1)*z1&
     &                                          + bint(2)*z2&
     &                                          + bint(3)*z3&
     &                                          + bint(4)*z4&
     &                                          + bint(5)*z5
                    uorig(ival) = uvec(ival)               ! put predict into storage for testing
                  enddo

!                 do ival = 1,nvecLen
!                   uvec(ival) = bint(1)*ustage(ival,3)
!    &                         + bint(2)*ustage(ival,4)
!    &                         + bint(3)*ustage(ival,5)
!    &                         + bint(4)*ustage(ival,6)
!    &                         + bint(5)*ustage(ival,2)
!                   uvec(ival) = bint(1)*ustage(ival,4)
!    &                         + bint(2)*ustage(ival,5)
!    &                         + bint(3)*ustage(ival,6)
!    &                         + bint(4)*ustage(ival,2)
!    &                         + bint(5)*ustage(ival,3)
!                   uorig(ival) = uvec(ival)               ! put predict into storage for testing
!                 enddo
                endif
!-------------------------------NEWTON ITERATION----------------------
                do k = 1,20

                  icount = icount + 1

                  do ival = 1,nvecLen
                    uveciter(ival) = uvec(ival) !store old uvec
                  enddo

                  call RHS(uvec,resI(1,L+1),resE(1,L+1),dt,ep,iprob,sigma,rho,beta) !get res() for newton iteration

                  do ival = 1,nvecLen
                    Rnewton(ival) =uvec(ival)-aI(L+1,L+1)*resI(ival,L+1)-usum(ival)
                  enddo
           
                  call Jacobian(uvec,xjac,dt,ep,aI(L+1,L+1),iprob,nvecLen,sigma,rho,beta)
                  call Invert_Jacobian(nvecLen,xjac,xjacinv)

!                 call Jacobian_CSR(uvec,dt,ep,aI(L+1,L+1),iprob,nvecLen,sigma,rho,beta)
!                 call ilu0(nvecLen, aJac, jaJac, iaJac, aluJac, jluJac, juJac, iw, ierr) 

                  !  Form LU-decomposition
                  !  Back-Solve to get solution
                  
                  do i = 1,nvecLen
                    do j = 1,nvecLen
                      uvec(i) = uvec(i) - xjacinv(i,j)*Rnewton(j)     !u^n+1=u^n-J^-1*F
                    enddo
                  enddo

                  tmp = 0.0_wp
                  do ival = 1,nvecLen
                     tmp = tmp + abs(uvec(ival)-uveciter(ival)) !check accuracy of zeros
                  enddo
                  if(tmp.lt.1.e-12) go to 160                 !  kick out of newton iteration

                enddo
!-----------------------------------END NEWTON ITERATION--------------
  160           continue
 
              if((k.gt.maxiter(L+1)).and.(ktime.ne.1))maxiter(L+1)=k

              do ival = 1,nvecLen                    !  write the solution into a storage register
                ustage(ival,L+1) = uvec(ival)
              enddo

              xnorm = 0.0_wp                                     !  assess how good initial guess was
              snorm = 0.0_wp                                     !  assess how good initial guess was
              do ival = 1,nvecLen
                tmpD = (uvec(ival) - uorig(ival))
                xnorm = xnorm + tmpD*tmpD
                snorm = snorm + uvec(ival)*uvec(ival)
              enddo                                         
              xnorm = sqrt(xnorm/snorm)
              stageE(L+1) = stageE(L+1) + xnorm
              stageI(L+1) = stageI(L+1) + 1.0_wp*k

              elseif(L.eq.nrk) then
            
                do ival = 1,nvecLen
                  uvec(ival) = uveco(ival)
                enddo

                do ival = 1,nvecLen
                  do LL = 1,nrk 
                    uvec(ival) = uvec(ival) + bI(LL)*resI(ival,LL)+bE(LL)*resE(ival,LL)
                  enddo
                enddo
!---------------------------PREDICTED ERROR---------------------------
                if(time.lt.tfinal-1.0d-11)then               
                do ival = 1,nvecLen
                  errvec(ival) = 0.0_wp
                  do LL = 1,nrk 
                    errvec(ival) = errvec(ival)& 
     &                       + dt*( (bE(LL)-bEH(LL))*resE(ival,LL)+(bI(LL)-bIH(LL))*resI(ival,LL) )!?dt?
                  enddo
                  errvec(ival) = abs(errvec(ival))
                enddo
                rat = 1.0_wp
                endif                            
!-----------------------END PREDICTED ERROR---------------------------
!-----------------------PREDICT NEXT STAGE VALUES---------------------
!**About 100 different kinds
!**Note that ipred=2 is accomplished elsewhere         
                if(ipred.eq.1)then !begin with dense output
                do K=2,nrk
                  do ival = 1,nvecLen
                    predvec(ival,K) = uvec(ival)
                  enddo
                enddo

                elseif(ipred.eq.3)then
                do K=2,nrk
                  do ival = 1,nvecLen
                    predvec(ival,K) = uveco(ival)
                    the = (1.+cI(K)*rat)
                    do LL = 1,nrk
                      bb = bD(LL,1)*the&
     &                   + bD(LL,2)*the*the&
     &                   + bD(LL,3)*the*the*the&
     &                   + bD(LL,4)*the*the*the*the
                    predvec(ival,K) = predvec(ival,K) + bb*resI(ival,LL)
                    enddo
                  enddo
                enddo

                elseif(ipred.eq.4.or.ipred.eq.5)then
!              stage value predictors
!  U^(n+1,i+1) =                  \sum_{k=0}^{order} (etah_{i k}*r^(k)) * U^{n-1} +
!                \sum_{j=1}^{s-1} \sum_{k=0}^{order} ( BBh_{ijk}*r^(k)) * U^{n,j+1}
                do ival = 1,nvecLen
                  do inew=2,nrk
                    predvec(ival,inew) = 0.0_wp
                    do j = 1,nrk
                      the = 1.0_wp
                      bb = svpB(inew,j,0)&
                        + svpB(inew,j,1)*the&
                        + svpB(inew,j,2)*the*the&
                        + svpB(inew,j,3)*the*the*the&
                        + svpB(inew,j,4)*the*the*the*the
                      predvec(ival,inew) = predvec(ival,inew)& 
     &                                   + bb*ustage(ival,j)
                    enddo
                  enddo
                enddo
                endif
!-------------------END PREDICT NEXT STAGE VALUES---------------------
              endif
            enddo                                          
!-----------------------------END RK LOOP-----------------------------
            do ival = 1,nvecLen
              errvecT(ival) = errvecT(ival) + errvec(ival)
            enddo                                         
            t = t + dt                  !increment time
            if(t.ge.tfinal) go to 100   !if end time then exit loop

          enddo                                          
!-----------------------END TIME ADVANCEMENT LOOP---------------------
  100  continue
       
          cost(iDT) = log10((nrk-1)/dto)    !  nrk - 1 implicit stages

          totalerror  = 0.0_wp
          totalerrorP = 0.0_wp
          do ival = 1,nvecLen
            tmpvec(ival) = abs(uvec(ival)-uexact(ival))
            if(tmpvec(ival).eq.0.0)tmpvec(ival)=1.0d-15
            totalerror  = totalerror  + tmpvec(ival)**2 
            totalerrorP = totalerrorP + errvecT(ival)**2 
          enddo
          totalerror = sqrt(totalerror/nvecLen)

          do ii = 1,nveclen
            error1(iDT,ii)  = log10(tmpvec(ii))
            error1P(iDT,ii) = log10(errvecT(ii))
            write(49+ii,*)cost(iDT),error1(iDT,ii)
            write(59+ii,*)cost(iDT),error1P(iDT,ii)
          enddo
  
         enddo
!----------------------------END TIMESTEP LOOP------------------------
           jsamp = 41 

!          if(icase.eq.1)then                        !  make sure that data has not hit machine precision
!            jsamp = 51      
!          elseif(icase.eq.2)then
!            jsamp = 51      
!          elseif(icase.eq.3)then
!            jsamp = 51      
!          elseif(icase.eq.4)then
!            jsamp = 51      
!          endif
  
           do i=1,isamp
             sig(i) = 0.0_wp
           enddo

!--------------------------OUTPUTS------------------------------------
           !**GATHER OUTPUT VALUES
               do ii = 1,nveclen
                 call fit(cost,error1(:,ii),jsamp,sig,0,a(ii),&
     &       b(ii),siga1(ii),sigb1(ii),chi2(ii),q)
                 call fit(cost,error1P(:,ii),jsamp,sig,0,a(ii+nveclen),&
     &       b(ii+nveclen),siga1(ii+nveclen),sigb1(ii+nveclen),&
     &       chi2(ii+nveclen),q)
               enddo
           
           !**OUTPUT TO TERMINAL**
           write(*,60,advance="no")ep
           do ii = 1,nveclen*2-1
             write(*,60,advance="no")a(ii),b(ii)
           enddo
           write(*,60)a(nveclen*2),b(nveclen*2)

           epsave(jepsil) = log10(ep)
           do ii = 1,nveclen
            b1save(jepsil,ii) = -b(ii)
            b1Psave(jepsil,ii) = -b(ii+nveclen)
           enddo
         enddo           
!---------------------END STIFFNESS LOOP------------------------------
         !**OUTPUT CONVERGENCE VS STIFFNESS**
         filename=trim(probname(iprob))//'_'//trim(casename)//'_conv'//ext
         open(35,file=trim(fileloc)//filename)
         do ii=1,nveclen
           write(35,*)'zone T = "Var ',ii,': Implicit",'
           do j=1,jactual
             write(35,50)epsave(j),b1save(j,ii)
           enddo
           write(35,*)'zone T = "Var ',ii,': Predicted",'
           do j=1,jactual
             write(35,50)epsave(j),b1Psave(j,ii)
           enddo
         enddo

         !**OUTPUT TO TERMINAL**
         tmp = 1.0_wp*icount/jcount
         write(*,*)'average iterations per step',tmp
           
         do istage = 2,nrk
           stageE(istage) = (nrk-1)*stageE(istage)/jcount
           stageI(istage) = (nrk-1)*stageI(istage)/jcount
         enddo
         write(*,*)'error of initial guess is '
         do istage = 2,nrk
           write(*,*) istage,stageE(istage),maxiter(istage),stageI(istage)
         enddo
!----------------------END OUTPUTS------------------------------------
        enddo                                       
!----------------------END PROBLEMS LOOP------------------------------
       enddo                                    
!----------------------END ALGORITHMS LOOP----------------------------
   60  format( e12.5,1x,12(f8.3,1x))
   50  format( 10(e12.5,1x))

      END PROGRAM test_cases


   
