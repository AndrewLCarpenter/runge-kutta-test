!##############################################################################
!######## Program to test Runge-Kutta IMEX schemes on various problems ########
!##############################################################################
!
!----------------------------PROBLEM LIST--------------------------------------
!     problem 1) van der Pol (Hairer II, pp 403)
!     problem 2) Pureschi and Russo 
!     problem 3) Dekker 7.5.2 pp. 215 (Kaps problem   : Index 1)
!     problem 4) Dekker 7.5.1 pp. 214 (Kreiss' problem: Index 2)
!     WARNING: Kreiss has problem with algebraic variable
!     problem 5) Lorenz
!     WARNING: Lorenz has problem as well with something and wont converge
!     problem 6) Burger's
!     problem 7) Black-Scholes
!------------------------------------------------------------------------------
!
!-----------------------------REQUIRED FILES-----------------------------------
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! POLY_FIT_MOD.F90          *CONTAINS SUBROUTINES TO OUTPUT DATA
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES USED IN THE PROGRAM
! CSR_VARIABLES.F90         *?
! INVERT_JACOBIAN.F90       *INVERTS JACOBIAN 
! ALLOCATE_CSR_STORAGE.F90  *ALLOCATES STOREAGE OF VARIABLES IN CSR
! NEWTON_ITERATION.F90      *PERFORMS NEWTON ITERATION
! STAGE_VALUE_PREDICTOR.F90 *PREDICTS NEXT STAGE VALUES FOR NEWTON ITERATIONS
! RUNGE_KUTTA.F90           *CONTAINS RK CONSTANTS
! VANDERPOL.F90             *PROBLEM CONSTANTS FOR VANDERPOL
! PURESCHI.F90              *PROBLEM CONSTANTS FOR PURESCHI & RUSSO
! KAPS.F90                  *PROBLEM CONSTANTS FOR KAPS
! KREISS.F90                *PROBLEM CONSTANTS FOR KREISS'
! PROBLEMSUB.F90            *DEFINES WHICH PROBLEM IS RELATED TO USER INPUT
!------------------------------------------------------------------------------
!
!--------------------------HOW TO USE------------------------------------------
!**Adding test cases:
!     Subroutine file can be added to the directory with form of vanderpol.f90,
!     etc. problemsub.f90 must be updated for the problem to be included
!**Compiling
!     Compile all files listed above together in the same directory
!**Running test cases
!     Program will prompt user for which predictor, case, and problem
!     to run.
!**Output results
!     The program outputs results to files with the form:
! /"probname"/"casename"/"probname"_"casename"_"variable number".dat
! /"probname"/"casename"/"probname"_"casename"_"variable number"P.dat
! /"probname"/"casename"/"probname"_"casename"_conv.dat
!
!********************************BEGIN PROGRAM*********************************
      program test_cases

      use precision_vars
      use output_module
      use Stage_value_module
      use control_variables
      use runge_kutta
      use Newton
!------------------------------VARIABLES---------------------------------------
      implicit none     
!------------------------------PARAMETERS--------------------------------------
      integer, parameter :: isamp=71           !?
      integer, parameter :: jmax=81            !?
      integer, parameter :: jactual=81         !?
!-----------------------------VARIABLES----------------------------------------  
      !internal variables
      real(wp)                              :: itmp,t,time,dto,tt
      real(wp), dimension(:,:), ALLOCATABLE :: ustage,predvec 
      real(wp), dimension(:),   ALLOCATABLE :: uveco,uveciter,uorig 
      real(wp), dimension(:),   ALLOCATABLE :: errvec,errvecT,tmpvec 
      real(wp), dimension(:,:), ALLOCATABLE :: b1save,b1Psave
      real(wp), dimension(jmax)             :: epsave
      real(wp)                             :: cputime1,cputime2
     
      !user inputs
      integer            :: ipred,problem,cases      
             
      !do loops
      integer            :: icase,i,iprob,jepsil,iDT 
      integer            :: j,k,istage,ktime,L,LL
      
      !counters      
      integer            :: icount,jcount                    
      
      !problemsub variables
      integer                               :: programStep
      real(wp), dimension(:),   ALLOCATABLE :: uvec
      real(wp)                              :: ep
      real(wp), dimension(:),   ALLOCATABLE :: uexact
      real(wp)                              :: dt
      integer                               :: nveclen
      real(wp)                              :: tfinal     
      real(wp), dimension(:,:), ALLOCATABLE :: resE,resI   

      !data out variables
      real(wp), dimension(isamp)            :: cost      
      real(wp), dimension(:,:), ALLOCATABLE :: error,errorP    
      integer                               :: jsamp
      real(wp), dimension(isamp)            :: sig
      real(wp), dimension(:),   ALLOCATABLE :: b
      real(wp), dimension(is)               :: stageE,stageI,maxiter
      
      !newton_iteration variables
      real(wp), dimension(:), ALLOCATABLE :: usum
      
      !dummy variables
      real(wp), dimension(:,:), ALLOCATABLE :: xjac
              
!-----------------------------USER INPUT---------------------------------------
      write(*,*)'what is ipred?' !input predictor number
      read(*,*)ipred
!     write(*,*)'what is case?'  !input range of runge kutta cases
!     read(*,*)cases
      write(*,*)'which problem?' !input problem number
      read(*,*)problem
!-------------------------ALGORITHMS LOOP--------------------------------------
!     do icase = cases,cases
      do icase = 15,15
        
        !**initilizations?**
        stageE(:) = 0.0_wp
        stageI(:) = 0.0_wp
        maxiter(:)= 0
        icount = 0                                  !cost counters
        jcount = 0                                  !cost counters
        
        !**GET RK COEFFICIENTS**
        call rungeadd(icase)

!--------------------------PROBLEMS LOOP---------------------------------------
        do iprob = problem,problem
        
          call cpu_time(cputime1)

          !**CALL FOR PROBNAME & NVECLEN**
          programStep=-1
          call problemsub(iprob,programStep,nveclen,uvec,ep, &
     &             uexact,dt,tfinal,iDT,tt,resE(:,1),resI(:,1),aI(1,1),xjac)
                              
          !**ALLOCATE VARIABLES**
          !problemsub
          AllOCATE(uvec(nveclen),uexact(nveclen))
          ALLOCATE(resE(nveclen,is),resI(nveclen,is))

          !data outs
          ALLOCATE(error(isamp,nveclen),errorP(isamp,nveclen))
          ALLOCATE(b(nveclen*2))
        
          !Newton_iteration
          ALLOCATE(usum(nveclen))
       
          !internal
          ALLOCATE(ustage(nveclen,is),predvec(nveclen,is))
          ALLOCATE(uveco(nveclen),uveciter(nveclen),uorig(nveclen))
          ALLOCATE(errvec(nveclen),errvecT(nveclen),tmpvec(nveclen))
          ALLOCATE(b1save(jmax,nveclen),b1Psave(jmax,nveclen))
          
          !dummy variables
          ALLOCATE(xjac(nveclen,nveclen))

          !**INIT. FILE PATH FOR OUTPUTS**
          call create_file_paths
          call output_names                           

!         call Allocata_CSR_Storage(problem,nveclen)

!--------------------------STIFFNESS LOOP--------------------------------------
          do jepsil = 1,jactual,1     
                         
            itmp = 11 - jmax/jactual         !used for 81 values of ep
            ep = 1.0_wp/10**((jepsil-1)/(itmp*1.0_wp))           
            
            !**INIT. OUTPUT FILES**
            call init_output_files(nveclen,ep)
!--------------------------TIMESTEP LOOP----------------------------------------
            do iDT = 1,isamp,1      
             

              !**INITIALIZE PROBLEM INFORMATION**
              programStep=0
              call problemsub(iprob,programStep,nveclen,uvec,ep, &
     &             uexact,dt,tfinal,iDT,tt,resE(:,1),resI(:,1),aI(1,1),xjac)  

              dto = dt        !store time step
              t = 0.0_wp      !init. start time

              !**INIT. ERROR VECTOR**
              errvecT(:) = 0.0_wp

                do i = 1,ns            !initialize stage value preditor
                  predvec(:,i) = uvec(:)
                enddo
!--------------------------TIME ADVANCEMENT LOOP-------------------------------
              do ktime = 1,1000000                      
                if(t+dt > tfinal)dt = tfinal-t+1.0e-11_wp !check if dt>tfinal
                tt=t+Ci(1)*dt
                !**STORE VALUES OF UVEC**
                uveco(:) = uvec(:)

                jcount = jcount + (ns-1)    !keep track of total RK stages  
!--------------------------------RK LOOP---------------------------------------
                programStep=1
                call problemsub(iprob,programStep,nveclen, &
     &     uvec,ep,uexact,dt,tfinal,iDT,tt,resE(1,1),resI(1,1),aI(1,1),xjac) 
     
                ustage(:,1) = uvec(:)

                do L = 2,ns
                tt=t+Ci(L)*dt

                  usum(:) = uveco(:)
                  do LL = 1,L-1 
                    usum(:) = usum(:) + aI(L,LL)*resI(:,LL) &
     &                                + aE(L,LL)*resE(:,LL)
                  enddo
               
                  if(ipred==2) predvec(:,L)=uvec(:)!previous guess as starter
                  if(ipred/=2) uvec(:) = predvec(:,L)  

!---------------BEG NEWTON ITERATION ------------------------------------------
                  call Newton_Iteration(uvec,iprob,L,ep,dt,&
     &                             nveclen,iDT,tt,resE,resI,aI(L,L),usum,icount,k)
!---------------END NEWTON ITERATION-------------------------------------------

                  ustage(:,L) = uvec(:)     !  Save the solution at each stage
                  ! Fill in resE and resI with the converged data
                  programStep=2
                  call problemsub(iprob,programStep,nveclen,&
     &     uvec,ep,uexact,dt,tfinal,iDT,tt,resE(1,L),resI(1,L),aI(1,1),xjac)

                  if(ipred/=2)call Stage_Value_Predictor(ipred,L,ustage,&
     &                                      predvec,uvec,uveco,alpha,ktime)

                  if((k > maxiter(L)).and.(ktime/=1)) maxiter(L) = k

                  stageE(L) = stageE(L) + xnorm(uvec,predvec(:,L))
                  stageI(L) = stageI(L) + 1.0_wp*k

                  if(ipred==2)call Stage_Value_Predictor(ipred,L,ustage,&
     &                                      predvec,uvec,uveco,alpha,ktime)
                enddo
!-----------------------------END of A_{k,j} portion of RK LOOP----------------
     
                uvec(:) = uveco(:)

                do LL = 1,ns 
                  uvec(:) = uvec(:) + bI(LL)*resI(:,LL)+bE(LL)*resE(:,LL)
                enddo

!-----------------------Final Sum of RK loop using the b_{j}-------------------

                ! ERROR ESTIMATE
                if(t <= tfinal-1.0e-11_wp)then               
                  errvec(:) = 0.0_wp
                  do LL = 1,ns 
                    errvec(:) = errvec(:) + (bE(LL)-bEH(LL))*resE(:,LL) &
     &                                    + (bI(LL)-bIH(LL))*resI(:,LL)

                  enddo
                  errvec(:) = abs(errvec(:))
                endif                            
  
                errvecT(:) = errvecT(:) + errvec(:)
  
                t = t + dt                  !increment time
                if(t >= tfinal) exit
              enddo                                          
!-----------------------END TIME ADVANCEMENT LOOP------------------------------
     
              cost(iDT) = log10((ns-1)/dto)    !  ns - 1 implicit stages
              
              call output_conv_error(cost(iDT),uvec,uexact,errvecT)
              
!
              tmpvec(:) = abs(uvec(:)-uexact(:))
              do i = 1,nvecLen
                if(tmpvec(i) == 0.0_wp)tmpvec(i)=1.0e-15_wp
              enddo
              error(iDT,:)  = log10(tmpvec(:))

              errorP(iDT,:) = log10(errvecT(:))

            enddo
!----------------------------END TIMESTEP LOOP---------------------------------
!----------------------------OUTPUTS-------------------------------------------
            jsamp = 41 
            sig(:) = 0.0_wp

            call output_terminal_iteration(cost,error,errorP,jsamp,sig,0, &
     &           ep,nveclen,b)
                    
                    
            epsave(jepsil) = log10(ep)
            do i = 1,nveclen
              b1save(jepsil,i) = -b(i)
              b1Psave(jepsil,i) = -b(i+nveclen)
            enddo
          enddo  
!---------------------END STIFFNESS LOOP---------------------------------------
          !**OUTPUT CONVERGENCE VS STIFFNESS**
          call output_conv_stiff(nveclen,jactual,epsave, &
    &                           b1save,b1Psave)
          
          !**OUTPUT TO TERMINAL**
          call output_terminal_final(icount,jcount,stageE,stageI,maxiter)
          
          call cpu_time(cputime2)
          write(*,*)'Total time elapsed for this case: ',cputime2-cputime1,'sec'  

!----------------------END OUTPUTS---------------------------------------------
           !**DEALLOCATE VARIABLES**
          !problemsub
          DEAllOCATE(uvec,uexact)
          DEALLOCATE(resE,resI)
          
          !data out
          DEALLOCATE(error,errorP)
          DEALLOCATE(b)
          
          !Newton_iteration
          DEALLOCATE(usum)
          
          !internal
          DEALLOCATE(ustage,predvec)
          DEALLOCATE(uveco,uveciter,uorig)
          DEALLOCATE(errvec,errvecT,tmpvec)
          DEALLOCATE(b1save,b1Psave)
          
          !dummy variables
          DEALLOCATE(xjac)
        enddo                                       
!----------------------END PROBLEMS LOOP---------------------------------------
      enddo                                 
!----------------------END ALGORITHMS LOOP-------------------------------------
!----------------------END PROGRAM---------------------------------------------
      END PROGRAM test_cases
