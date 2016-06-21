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
! NEWTON_ITERATION.F90
! PROBLEMSSUB.F90
! ...
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
      use poly_fit_Mod

      implicit none     

!------------------------------PARAMETERS----------------------------
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
      integer             :: programStep                 !counts step of program

      real(wp)            :: dt,itmp,ep,t,sigma,rho,beta,tfinal,ttE
      real(wp)            :: ttI,Z,time,rat,bb,xnum,xden,z1,z2,z3,z4,z5
      real(wp)            :: the,the1,the2,the3,the4,the5,the6,the7,the8,q
      real(wp)            :: tmp,xnorm,snorm,tmpD,dto,totalerror,totalerrorp

      real(wp), dimension(is,is) :: aE,aI,alpha,al3N,al3D,al4N,al4D
      real(wp), dimension(is) :: bE,bI,cE,cI,bEH,bIH,stageE,stageI
      real(wp), dimension(is) :: maxiter,bint
      real(wp), dimension(is,4) :: bD
      real(wp), dimension(is,is,0:is) :: svpB
      real(wp), dimension(ivarlen,is) :: ustage,resE,resI,predvec
      real(wp), dimension(ivarlen) :: uvec,uveco,usum,uveciter,uorig
      real(wp), dimension(ivarlen) :: uexact,Rnewton,errvec,errvecT,tmpvec
      real(wp), dimension(ivarlen,ivarlen) :: xjacinv,xjac
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

      logical           :: dirExists,filex,test    !check if directory exists
     
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
      do icase = cases,cases

        icount = 0                                  !cost counters
        jcount = 0                                  !cost counters
        
        !**GET RK COEFFICIENTS**
        call rungeadd(aE,aI,bE,bI,cE,cI,nrk,bEH,bIH,icase,bD, &   
     &   svpB(1,1,0),alpha,al3N,al3D,al4N,al4D,casename)

        !**initilizations?**
        do i = 1,nrk
          stageE(i) = 0.0_wp
          stageI(i) = 0.0_wp
          maxiter(i)= 0
        enddo
!--------------------------PROBLEMS LOOP------------------------------
        do iprob = problem,problem
          !**INIT. FILE PATH FOR OUTPUTS**                           
          fileloc=trim(probname(iprob))//'/'//trim(casename)//'/'
          inquire( file=trim(fileloc)//'/.', exist=dirExists ) 
          if (.not.dirExists)call system('mkdir -p '//trim(fileloc)) ! create directory
          
!!          call Allocata_CSR_Storage(problem,nveclen)
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
            programStep=0
            call problemsub(iprob,programStep,uvec,ep,uexact,dt,nveclen,tfinal,iDT)   !,resE(1,1),resI(1,1),aI(1,1),xjac)
           ! call INIT(uvec,uexact,dt,iDT,tfinal,ep,nvecLen,iprob,sigma,rho,beta)     
            dto = dt        !store time step
            t = 0.0_wp      !init. start time

            !**INIT. ERROR VECTOR**
            errvecT(:) = 0.0_wp

            do i = 1,nrk            !initialize stage value preditor
              predvec(:,i) = uvec(:)
            enddo
!--------------------------TIME ADVANCEMENT LOOP----------------------
            do ktime = 1,100000000                      

              if(t+dt.gt.tfinal)dt = tfinal-t + 1.0e-11_wp !check if dt>tfinal
              
              !**STORE VALUES OF UVEC**
              uveco(:) = uvec(:)

              jcount = jcount + (nrk-1)    !keep track of total RK stages 

!--------------------------------RK LOOP------------------------------
              ttI = t + cI(1)*dt
              ttE = t + cE(1)*dt

              programStep=1
              call problemsub(iprob,programStep,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE(1,1),resI(1,1))
!             call RHS(uvec,resI(1,1),resE(1,1),dt,ep,iprob,sigma,rho,beta)
              ustage(:,1) = uvec(:)

              do L = 2,nrk

                ttI = t + cI(L)*dt   
                ttE = t + cE(L)*dt

                usum(:) = uveco(:)
                do LL = 1,L-1 
                  usum(:) = usum(:) + aI(L,LL)*resI(:,LL) &
                          &         + aE(L,LL)*resE(:,LL)
                enddo

!---------------BEG NEWTON ITERATION------------------------------------

                uvec(:) = predvec(:,L)

                call Newton_Iteration(uvec,iprob,L,ep,dt,nveclen,iDT,resE,resI,aI(L,L),usum,icount,k)

!---------------END NEWTON ITERATION------------------------------------
             
                ustage(:,L) = uvec(:)                !  Save the solution at each stage
                ! Fill in resE and resI with the converged data
                programStep=2
                call problemsub(iprob,programStep,uvec,ep,uexact,dt,nveclen,tfinal,iDT,resE(1,L),resI(1,L))

                call Stage_Value_Predictor(L,nrk,ustage,predvec)

                if((k.gt.maxiter(L)).and.(ktime.ne.1)) maxiter(L) = k
 
                xnorm = sqrt(dot_product(uvec(:)-predvec(:,L),uvec(:)-predvec(:,L)) /  &
                      &      dot_product(uvec(:)             ,uvec(:)             ) )

                stageE(L) = stageE(L) + xnorm
                stageI(L) = stageI(L) + 1.0_wp*k

            enddo   !-----------------------------END of A_{k,j} portion of RK LOOP---------------

            uvec(:) = uveco(:)
            do LL = 1,nrk 
              uvec(:) = uvec(:) + bI(LL)*resI(:,LL)+bE(LL)*resE(:,LL)
            enddo
!----------------------------- Final Sum of RK loop using the b_{j} 

            ! ERROR ESTIMATE
            if(time.lt.tfinal-1.0e-11_wp)then               
              errvec(:) = 0.0_wp
              do LL = 1,nrk 
                errvec(:) = errvec(:) + (bE(LL)-bEH(LL))*resE(:,LL) &
                                    & + (bI(LL)-bIH(LL))*resI(:,LL)
              enddo
              errvec(:) = abs(errvec(:))
              rat = 1.0_wp
            endif                            

            errvecT(:) = errvecT(:) + errvec(:)

            t = t + dt                  !increment time
            if(t.ge.tfinal) exit

          enddo                                          
!-----------------------END TIME ADVANCEMENT LOOP---------------------
       
          cost(iDT) = log10((nrk-1)/dto)    !  nrk - 1 implicit stages

          totalerror  = 0.0_wp
          totalerrorP = 0.0_wp
          do ival = 1,nvecLen
            tmpvec(ival) = abs(uvec(ival)-uexact(ival))
            if(tmpvec(ival) == 0.0_wp)tmpvec(ival)=1.0e-15_wp
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


   
