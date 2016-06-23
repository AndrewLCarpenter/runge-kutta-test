!###############################################################################
!######## Program to test Runge-Kutta IMEX schemes on various problems #########
!###############################################################################
!
!----------------------------PROBLEM LIST---------------------------------------
!     problem 1) van der Pol (Hairer II, pp 403)
!     problem 2) Pureshi and Russo 
!     problem 3) Dekker 7.5.2 pp. 215 (Kaps problem   : Index 1)
!     problem 4) Dekker 7.5.1 pp. 214 (Kreiss' problem: Index 2)
!     WARNING: Kreiss has problem with algebraic variable
!     problem 5) Lorenz
!     problem 6) Burger's
!     problem 7) Black-Scholes
!-------------------------------------------------------------------------------
!
!-----------------------------REQUIRED FILES------------------------------------
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! POLY_FIT_MOD.F90          *CONTAINS SUBROUTINES TO OUTPUT DATA
! CSR_VARIABLES.F90         *?
! INVERT_JACOBIAN.F90       *INVERTS JACOBIAN (NOT USED?)
! ALLOCATE_CSR_STORAGE.F90  *ALLOCATES STOREAGE OF VARIABLES IN CSR
! JACOB2.F90                *INVERTS JACOBIAN
! JACOBIAN.F90              *DEFINES JACOBIANS (NOT USED?)
! NEWTON_ITERATION.F90      *PERFORMS NEWTON ITERATION
! STAGE_VALUE_PREDICTOR.F90 *PREDICTS NEXT STAGE VALUES FOR NEWTON ITERATIONS
! RUNGEADD.F90              *CONTAINS RK CONSTANTS
! VANDERPOL.F90             *PROBLEM CONSTANTS FOR VANDERPOL
! PURESCHI.F90              *PROBLEM CONSTANTS FOR PURESCHI & RUSSO
! KAPS.F90                  *PROBLEM CONSTANTS FOR KAPS
! KREISS.F90                *PROBLEM CONSTANTS FOR KREISS'
! PROBLEMSUB.F90            *DEFINES WHICH PROBLEM IS RELATED TO USER INPUT
!-------------------------------------------------------------------------------
!
!--------------------------HOW TO USE-------------------------------------------
!**Adding test cases:
!      Subroutine file can be added to the directory with form of vanderpol.f90,
!      etc. problemsub.f90 must be updated for the problem to be included
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
!********************************BEGIN PROGRAM**********************************
      program test_cases

      use precision_vars
      use poly_fit_Mod
!------------------------------VARIABLES----------------------------------------
      implicit none     
!------------------------------PARAMETERS---------------------------------------
      integer, parameter :: is=9               !max constant length
      integer, parameter :: ivarlen=4          !max variable length
      integer, parameter :: isamp=71           !?
      integer, parameter :: jmax=81            !?
      integer, parameter :: jactual=81         !?
!-----------------------------VARIABLES-----------------------------------------   
      !user inputs
      integer          :: ipred,cases,problem      
             
      !do loops
      integer          :: icase,i,iprob,jepsil,iDT,ival  
      integer          :: j,k,istage,ktime,L,LL
      
      !counters      
      integer          :: icount,jcount                    
      
      !problemsub variables
      integer                         :: programStep
      character(len=9)                :: probname
      real(wp), dimension(:), ALLOCATABLE     :: uvec !ivarlen
      real(wp)                        :: ep
      real(wp), dimension(:), ALLOCATABLE     :: uexact !ivarlen
      real(wp)                        :: dt
      integer                         :: nveclen
      real(wp)                        :: tfinal     
      real(wp), dimension(:,:), ALLOCATABLE  :: resE,resI !ivarlen,is
      
      !rungeadd variables
      real(wp),   dimension(is,is)      :: aE,aI
      real(wp),   dimension(is)         :: bE,bI,cE,cI
      integer                           :: nrk
      real(wp),   dimension(is)         :: bEH,bIH
      real(wp),   dimension(is,ivarlen) :: bD 
      real(wp),   dimension(is,is,0:is) :: svpB 
      real(wp),   dimension(is,is)      :: alpha,al3N,al3D,al4N,al4D
      character(len=25)                 :: casename     

      !fit variables
      real(wp), dimension(isamp)         :: cost      
      real(wp), dimension(:,:), ALLOCATABLE  :: error,errorP    !isamp,ivarlen  
      integer                            :: jsamp
      real(wp), dimension(isamp)         :: sig
      real(wp), dimension(:), ALLOCATABLE      :: a,b,siga1,sigb1,chi2 !ivarlen*2
      real(wp)                           :: q
      
      !newton_iteration variables
      !add variables
      
      !internal variables**includes newton
      real(wp)                          :: itmp,t,time,rat,tmp,xnorm
      real(wp)                          :: dto,totalerror,totalerrorp
      real(wp), dimension(is)           :: stageE,stageI,maxiter,bint
      real(wp), dimension(:,:), ALLOCATABLE    :: ustage,predvec !ivarlen,is
      real(wp), dimension(:), ALLOCATABLE       :: uveco,usum,uveciter,uorig !ivarlen
      real(wp), dimension(:), ALLOCATABLE       :: errvec,errvecT,tmpvec !ivarlen
      real(wp), dimension(:,:), ALLOCATABLE  :: b1save,b1Psave !jmax,ivarlen
      real(wp), dimension(jmax)         :: epsave
      
      character(len=80) :: filename,fileloc   !file name and location
      character(len=4)  :: ext='.dat'         !file extension
      character         :: istr               !loop index placeholder

      logical           :: dirExists   !check if directory exists
         
!-----------------------------USER INPUT----------------------------------------
      write(*,*)'what is ipred?' !input predictor number
      read(*,*)ipred
      write(*,*)'what is case?'  !input range of runge kutta cases
      read(*,*)cases
      write(*,*)'which problem?' !input problem number
      read(*,*)problem
!-------------------------ALGORITHMS LOOP---------------------------------------
      do icase = cases,cases
        
        !**initilizations?**
        stageE(:) = 0.0_wp
        stageI(:) = 0.0_wp
        maxiter(:)= 0
        icount = 0                                  !cost counters
        jcount = 0                                  !cost counters
        
        !**GET RK COEFFICIENTS**
        call rungeadd(aE,aI,bE,bI,cE,cI,nrk,bEH,bIH,icase,bD, &   
     &       svpB(1,1,0),alpha,al3N,al3D,al4N,al4D,casename)

!--------------------------PROBLEMS LOOP----------------------------------------
        do iprob = problem,problem
          !**CALL FOR PROBNAME**
          programStep=-1
          call problemsub(iprob,programStep,probname,nveclen)!,uvec,ep,uexact,dt,&
    ! &    tfinal,iDT)
          !**ALLOCATE VARIABLES**
          !problemsub
          AllOCATE(uvec(nveclen),uexact(nveclen))
          ALLOCATE(resE(nveclen,is),resI(nveclen,is))
          
          !fit
          ALLOCATE(error(isamp,nveclen),errorP(isamp,nveclen))
          ALLOCATE(a(nveclen*2),b(nveclen*2),siga1(nveclen*2))
          ALLOCATE(sigb1(nveclen*2),chi2(nveclen*2))
          
          !Newton_iteration
          
          !internal
          ALLOCATE(ustage(nveclen,is),predvec(nveclen,is))
          ALLOCATE(uveco(nveclen),usum(nveclen),uveciter(nveclen),uorig(nveclen))
          ALLOCATE(errvec(nveclen),errvecT(nveclen),tmpvec(nveclen))
          ALLOCATE(b1save(jmax,nveclen),b1Psave(jmax,nveclen))

          !**INIT. FILE PATH FOR OUTPUTS**                           
          fileloc=trim(probname)//'/'//trim(casename)//'/'
          inquire( file=trim(fileloc)//'/.', exist=dirExists ) 
          if (.not.dirExists)call system('mkdir -p '//trim(fileloc)) ! create directory
          
!         call Allocata_CSR_Storage(problem,nveclen)
!--------------------------STIFFNESS LOOP---------------------------------------
          do jepsil = 1,jactual,1     
                         
            itmp = 11 - jmax/jactual         !used for 81 values of ep
            ep = 1.0_wp/10**((jepsil-1)/(itmp*1.0_wp))           
            
            !**INIT. OUTPUT FILES**
            do i = 1,ivarlen
              write(istr,"(I1.1)")i
              filename= &
     &           trim(probname)//'_'//trim(casename)//'_'//istr//ext
              open(49+i,file=trim(fileloc)//filename)
              write(49+i,*)'zone T = "ep = ',ep,'",'
              filename= &
     &           trim(probname)//'_'//trim(casename)//'_'//istr//'P'//ext
              open(59+i,file=trim(fileloc)//filename)
              write(59+i,*)'zone T = "ep = ',ep,'",'
            enddo
!--------------------------TIMESTEP LOOP----------------------------------------
            do iDT = 1,isamp,1                     

              !**INITIALIZE PROBLEM INFORMATION**
              programStep=0
              call problemsub(iprob,programStep,probname,nveclen,uvec,ep,uexact,dt,&
     &                        tfinal,iDT)

              dto = dt        !store time step
              t = 0.0_wp      !init. start time
              
              !**INIT. ERROR VECTOR**
              errvecT(:) = 0.0_wp

              do j = 1,nveclen
                do i = 1,nrk            !initialize stage value preditor
                  predvec(j,i) = uvec(j)
                enddo
              enddo
!--------------------------TIME ADVANCEMENT LOOP--------------------------------
              do ktime = 1,1000000                      
                !print*,uvec(:),ktime
                if(t+dt.gt.tfinal)dt = tfinal-t + 1.0e-11_wp !check if dt>tfinal
               
                !**STORE VALUES OF UVEC**
                uveco(:) = uvec(:)

                jcount = jcount + (nrk-1)    !keep track of total RK stages 

!--------------------------------RK LOOP----------------------------------------
                programStep=1
                call problemsub(iprob,programStep,probname,nveclen,uvec,ep,uexact,dt,&
     &                          tfinal,iDT,resE(1,1),resI(1,1))
                ustage(:,1) = uvec(:)
  
                do L = 2,nrk
                  usum(:) = uveco(:)
                  do LL = 1,L-1 
                    usum(:) = usum(:) + aI(L,LL)*resI(:,LL) &
     &                                + aE(L,LL)*resE(:,LL)
                  enddo
 
                  if(ipred.eq.2) predvec(:,L)=uvec(:)  ! previous guess as starter
                  if(ipred.ne.2) uvec(:) = predvec(:,L)  
!---------------BEG NEWTON ITERATION -------------------------------------------
                 ! print*,'newtonstart',uvec
                  call Newton_Iteration(uvec,iprob,L,ep,dt,nveclen,iDT,resE, &
     &                                  resI,aI(L,L),usum,icount,k)
                 ! print*,'newtonend  ',uvec,k

!---------------END NEWTON ITERATION--------------------------------------------
             
                  ustage(:,L) = uvec(:)       !  Save the solution at each stage
                  ! Fill in resE and resI with the converged data
                  programStep=2
                  call problemsub(iprob,programStep,probname,nveclen,uvec,ep,uexact,dt,&
     &                            tfinal,iDT,resE(1,L),resI(1,L))

                  if(ipred.ne.2)call Stage_Value_Predictor(L,nrk,ustage,predvec)

                  if((k.gt.maxiter(L)).and.(ktime.ne.1)) maxiter(L) = k

                  xnorm = sqrt(dot_product(uvec(:)-predvec(:,L),  &
     &                                    uvec(:)-predvec(:,L))/  &
     &                        dot_product(uvec(:),uvec(:)))

                  stageE(L) = stageE(L) + xnorm
                  stageI(L) = stageI(L) + 1.0_wp*k

                  if (ipred.eq.2) then
                    uvec(:)  = uveco(:)          ! put predict into start guess
                    do j = 2,L-1
                      uvec(:) = uvec(:) + alpha(L,j)*(ustage(:,j)-uvec(:))
                    enddo
                
                    if(L.eq.2.and.ktime.ne.1)then
                      bint(1) =  9.9518675264213746_wp 
                      bint(2) =  4.8366852488953721_wp 
                      bint(3) =-24.163405114569394_wp 
                      bint(4) = 14.152132944153401_wp
                      bint(5) =  0.94399768676237158_wp
                      uvec(:) = ustage(:,3) + bint(1)*(ustage(:,3)-ustage(:,3))&
     &                                      + bint(2)*(ustage(:,4)-ustage(:,3))&
     &                                      + bint(3)*(ustage(:,5)-ustage(:,3))&
     &                                      + bint(4)*(ustage(:,6)-ustage(:,3))&
     &                                      + bint(5)*(ustage(:,2)-ustage(:,3))
                    endif
                  endif
                enddo
!-----------------------------END of A_{k,j} portion of RK LOOP-----------------

                uvec(:) = uveco(:)
                do LL = 1,nrk 
                  uvec(:) = uvec(:) + bI(LL)*resI(:,LL)+bE(LL)*resE(:,LL)
                enddo

                
!-----------------------Final Sum of RK loop using the b_{j}--------------------

                ! ERROR ESTIMATE
                if(time.lt.tfinal-1.0e-11_wp)then               
                  errvec(:) = 0.0_wp
                  do LL = 1,nrk 
                    errvec(:) = errvec(:) + (bE(LL)-bEH(LL))*resE(:,LL) &
     &                                    + (bI(LL)-bIH(LL))*resI(:,LL)
                  enddo
                  errvec(:) = abs(errvec(:))
                  rat = 1.0_wp
                endif                            
  
                errvecT(:) = errvecT(:) + errvec(:)
  
                t = t + dt                  !increment time
                if(t.ge.tfinal) exit
              enddo                                          
!-----------------------END TIME ADVANCEMENT LOOP-------------------------------
           
              cost(iDT) = log10((nrk-1)/dto)    !  nrk - 1 implicit stages

              tmpvec(:) = abs(uvec(:)-uexact(:))
              do ival = 1,nvecLen
                if(tmpvec(ival) == 0.0_wp)tmpvec(ival)=1.0e-15_wp
              enddo
  
              totalerror  = sum( tmpvec(:)**2 )
              totalerrorP = sum( errvecT(:)**2 )
              totalerror = sqrt(totalerror/nvecLen)
              error(iDT,:)  = log10(tmpvec(:))
              errorP(iDT,:) = log10(errvecT(:))
              
              do i = 1,nveclen
                write(49+i,*)cost(iDT),error(iDT,i)
                write(59+i,*)cost(iDT),errorP(iDT,i)
              enddo
            enddo
!----------------------------END TIMESTEP LOOP----------------------------------
!----------------------------OUTPUTS--------------------------------------------
            jsamp = 41 
            sig(:) = 0.0_wp
            !**GATHER OUTPUT VALUES
            do i = 1,nveclen
         !     print*,'i',i,'nveclen',nveclen
         !     print*,cost
         !     print*,error(:,i)
         !     print*,jsamp
         !     print*,sig
         !     print*,0
         !     print*,a(i)
         !     print*,b(i)
         !     print*,siga1(i)
         !     print*,sigb1(i)
         !     print*,chi2(i)
         !     print*,q
              call fit(cost,error(:,i),jsamp,sig,0,a(i),&
     &         b(i),siga1(i),sigb1(i),chi2(i),q)
              call fit(cost,errorP(:,i),jsamp,sig,0,a(i+nveclen),&
     &         b(i+nveclen),siga1(i+nveclen),sigb1(i+nveclen),&
     &         chi2(i+nveclen),q)
            enddo

            !**OUTPUT TO TERMINAL**
            write(*,60,advance="no")ep
            do i = 1,nveclen*2-1
              write(*,60,advance="no")a(i),b(i)
            enddo
            write(*,60)a(nveclen*2),b(nveclen*2)

            epsave(jepsil) = log10(ep)
            do i = 1,nveclen
              b1save(jepsil,i) = -b(i)
              b1Psave(jepsil,i) = -b(i+nveclen)
            enddo
          enddo           
!---------------------END STIFFNESS LOOP----------------------------------------
          !**OUTPUT CONVERGENCE VS STIFFNESS**
          filename=trim(probname)//'_'//trim(casename)//'_conv'//ext
          open(35,file=trim(fileloc)//filename)
          do i=1,nveclen
            write(35,*)'zone T = "Var ',i,': Implicit",'
            do j=1,jactual
              write(35,50)epsave(j),b1save(j,i)
            enddo
            write(35,*)'zone T = "Var ',i,': Predicted",'
            do j=1,jactual
              write(35,50)epsave(j),b1Psave(j,i)
            enddo
          enddo

          !**OUTPUT TO TERMINAL**
          tmp = 1.0_wp*icount/jcount
          write(*,*)'average iterations per step',tmp
             
          stageE(2:) = (nrk-1)*stageE(2:)/jcount
          stageI(2:) = (nrk-1)*stageI(2:)/jcount

          write(*,*)'error of initial guess is '
          do istage = 2,nrk
            write(*,*) istage,stageE(istage),maxiter(istage),stageI(istage)
          enddo
!----------------------END OUTPUTS----------------------------------------------
           !**DEALLOCATE VARIABLES**
          !problemsub
          DEAllOCATE(uvec,uexact)
          DEALLOCATE(resE,resI)
          
          !fit
          DEALLOCATE(error,errorP)
          DEALLOCATE(a,b,siga1)
          DEALLOCATE(sigb1,chi2)
          
          !Newton_iteration
          
          !internal
          DEALLOCATE(ustage,predvec)
          DEALLOCATE(uveco,usum,uveciter,uorig)
          DEALLOCATE(errvec,errvecT,tmpvec)
          DEALLOCATE(b1save,b1Psave)
        enddo                                       
!----------------------END PROBLEMS LOOP----------------------------------------
      enddo                                    
!----------------------END ALGORITHMS LOOP--------------------------------------
   60 format( e12.5,1x,12(f8.3,1x))
   50 format( 10(e12.5,1x))
!----------------------END PROGRAM----------------------------------------------
      END PROGRAM test_cases
