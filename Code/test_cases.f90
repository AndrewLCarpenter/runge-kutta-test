!##############################################################################
!######## Program to test Runge-Kutta IMEX schemes on various problems ########
!##############################################################################
!
!----------------------------PROBLEM LIST--------------------------------------
!     problem 1 ) van der Pol (Hairer II, pp 403)
!     problem 2 ) Pureschi and Russo 
!     problem 3 ) Dekker 7.5.2 pp. 215 (Kaps problem   : Index 1)
!     problem 4 ) Dekker 7.5.1 pp. 214 (Kreiss' problem: Index 2)
!     problem 5 ) Lorenz attractor
!     problem 6 ) Rossler_Chaos (Wolf.Swift.Swinney.Vastano. Physica 16D,(1985),285-317
!     problem 7 ) Oregonator
!     problem 8 ) Brusselator
!     problem 9 ) Burgers
!     problem 10) Boscarino31 
!     problem 11) Broadwell Model
!------------------------------------------------------------------------------
!
!-----------------------------REQUIRED FILES-----------------------------------
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
! RUNGE_KUTTA.F90           *CONTAINS RK CONSTANTS
! OUTPUT_MODULE.F90         *CONTAINS ROUTINES TO OUTPUT DATA TO FILES AND TERMINAL
! STAGE_VALUE_MODULE.F90    *PREDICTS NEXT STAGE VALUES FOR NEWTON ITERATIONS
! PROBLEMSUB.F90            *DEFINES WHICH PROBLEM IS RELATED TO USER INPUT
! NEWTON.F90                *PERFORMS NEWTON ITERATION
!------------------------------------------------------------------------------
!
!----------------------------GLOBAL VARIABLES/ROUTINES-------------------------
! From precision_variables:
!   wp  -> working precision
! From output_module:
!   create_file_paths         -> Creates file paths for output files if they do not exsist
!   output_names              -> Outputs RK case, IMEX/IMPLICIT/EXPLICIT, and problem names to the terminal
!   init_output_files         -> Create zones in output files for Tecplot/Tec360
!   output_conv_stiff         -> Output data to convergence vs. stiffness files
!   output_terminal_iteration -> Outputs data to terminal every epsilon iteration
!   output_terminal_final     -> Outputs data to terminal at end of RK case/problem
!   output_conv_error         -> Outputs data to files initilized earilier in the program
!   write_time_depen_sol      -> Outputs time dependant solution
!   output_iDT_sol            -> Outputs solution for each iDT value for use in finding exact solution
!   write_exact_sol           -> Outputs exact solution to a file in the appropriate directory
! From Stage_value_module:
!   Stage_Value_Predictor -> returns appropriate data for predicting stage value
!   xnorm                 -> Function that takes a norm for stage value storage,  real(wp)
! From control_variables:
!   allocate_vars   -> Subroutine to allocate global variables
!   deallocate_vars -> Subroutine to deallocate global variables at end of problem loop
!   isamp       -> number of dt's,                               integer,                                              not modified
!   jmax        -> number of epsilon values,                     integer,                                              not modified 
!   uveco       -> newton iteration / usum storage (?),          real(wp), dimension(u-vector length),                     modified
!   errvec      -> temporary error storage,                      real(wp), dimension(u-vector length),                     modified
!   errvecT     -> error storage for each time step,             real(wp), dimension(u-vector length),                     modified
!   tmpvec      -> (uvec-uexact),                                real(wp), dimension(u-vector length),                     modified
!   resE        -> Explicit RHS vector,                          real(wp), dimension(u-vector length, num. stages),    not modified
!   resI        -> Implicit RHS vector,                          real(wp), dimension(u-vector length, num. stages),    not modified
!   error       -> solution error,                               real(wp), dimension(num. dt's, u-vector length),          modified
!   errorP      -> solution error, predicted,                    real(wp), dimension(num. dt's, u-vector length),          modified
!   b1save      -> convergence rate storage,                     real(wp), dimension(num. ep's, u-vector length),          modified
!   b1Psave     -> convergence rate storage, predicted,          real(wp), dimension(num. ep's, u-vector length),          modified
!   ustage      -> stage value predictor,                        real(wp), dimension(u-vector length, num. stages),        modified
!   predvec     -> stage value predictor,                        real(wp), dimension(u-vector length, num. stages),        modified
!   jactual     -> actual number of epsilon values,              integer,                                              not modified
!   uvec        -> Array containing variables,                   real(wp), dimension(u-vector length),                     modified
!   uexact      -> Array containing exact solution to variables, real(wp), dimension(u-vector length),                 not modified
!   b           -> convergence rate,                             real(wp), dimension(2 * u-vector length + num. eq's), not modified          
!   usum        -> Array containing summation from test_cases,   real(wp), dimension(u-vector length),                     modified
!   programstep -> string used in program_step_select,           character(len=80),                                    not modified 
!   errorL2     -> solution error, L2 norm,                      real(wp), dimension(num. dt's, num. eq's),                modified
!   b1L2save    -> convergence rate storage, L2 norm ,           real(wp), dimension(num. ep's, num. eq's),                modified
! From runge_kutta:
!   aE  -> Explicit a  constants, real(wp), dimension(max. num. stages, max. num. stages), not modified
!   aI  -> Implicit a  constants, real(wp), dimension(max. num. stages, max. num. stages), not modified
!   bE  -> Explicit b  constants, real(wp), dimension(max. num. stages),                   not modified
!   bI  -> Implicit b  constants, real(wp), dimension(max. num. stages),                   not modified
!   bEH -> Explicit bH constants, real(wp), dimension(max. num. stages),                   not modified
!   bIH -> Implicit bH constants, real(wp), dimension(max. num. stages),                   not modified
!   cI  -> Implicit c  constants, real(wp), dimension(max. num. stages),                   not modified
!   is  -> Maximum number stages, integer,                                                 not modified
!   ns  -> Actualy number stages, integer,                                                 not modified 
!   rungeadd -> Subroutine to set constants for particular RK schemes
! From Newton:
!   Newton_Iteration -> Subroutine to perform Newton iteration 
! From problemsub_mod
!   problemsub -> Subroutine to take in the problem number and call the appropriate problem
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
!**Find exact solution
!     To find the exact solution:
!     1) set temporal_splitting (in control_variables) to "IMPLICIT"
!     2) turn iDT_sol_flag=.true. and run desired problem and RK case.
!     3) complile convergence_testing and move executable to IMPLICIT/probname/casename
!     4) Take note of terminal output "Location = ##"
!     5) set that number to exact_sol_iDT=##
!     6) set exact_sol_flag=.true.
!     6) compile this program and run problem & RK case again
!     7) exact solution will be output to IMPLICIT/probname/casename, move/rename file to parent directory as appropriate
!     NOTE: remember to turn all flags==.FALSE. as soon as they are no longer needed
! 
!********************************BEGIN PROGRAM*********************************
      program test_cases

      use precision_vars,     only: wp
      use output_module,      only: create_file_paths, output_names,          & 
     &                              init_output_files, output_conv_stiff,     &
     &                              output_terminal_iteration,                &
     &                              output_terminal_final,output_conv_error,  &
     &                              write_time_depen_sol,output_iDT_sol,      &
     &                              write_exact_sol
      use Stage_value_module, only: Stage_Value_Predictor,xnorm
      use control_variables,  only: allocate_vars,deallocate_vars,isamp,jmax, &
     &                              uveco,errvec,errvecT,tmpvec,resE,resI,    &
     &                              error,errorP,b1save,b1Psave,ustage,       &
     &                              predvec,jactual,uvec,uexact,b,usum,       &
     &                              programstep,errorL2,b1L2save
      use runge_kutta,        only: aE,aI,bE,bI,bEH,bIH,cI,is,ns,rungeadd
      use Newton,             only: Newton_Iteration
      use problemsub_mod,     only: problemsub

      implicit none     
!-----------------------------VARIABLES----------------------------------------  
      !internal variables
      real(wp)                              :: itmp,t,dto,tt,cputime1,cputime2
      real(wp), dimension(jmax)             :: epsave
     
      integer            :: ipred,problem,cases                   !user inputs     
      integer            :: icase,i,iprob,jepsil,iDT,k,ktime,L,LL !do loops
      integer            :: icount,jcount                         !counters                     
      
      !problemsub variables
      real(wp)                              :: ep,dt,tfinal 
      integer                               :: nveclen,neq

      !data out variables
      real(wp), dimension(isamp)            :: cost         
      real(wp), dimension(is)               :: stageE,stageI,maxiter
      logical,  parameter                   :: time_sol_flag=.false. !True to output time solution
      integer,  parameter                   :: time_sol_iDT=45      !which dt value to output (You only want one timestep in the file,
                                                                    !                          system could be improved)
      logical,  parameter                   :: iDT_sol_flag=.false.               
      logical,  parameter                   :: exact_sol_flag=.false.
!     integer,  parameter                   :: exact_sol_iDT=67
                                                     
      
!      real(wp),dimension(40) :: tmpv
                          
!-----------------------------USER INPUT---------------------------------------
      write(*,*)'what is ipred?' !input predictor number
      read(*,*)ipred
      write(*,*)'what is case?'  !input range of runge kutta cases
      read(*,*)cases
      write(*,*)'which problem?' !input problem number
      read(*,*)problem
!-------------------------ALGORITHMS LOOP--------------------------------------  
      do icase = cases,cases
      
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
          programStep='INITIALIZE_PROBLEM_INFORMATION'
          call problemsub(iprob,nveclen,neq,ep,dt,tfinal,iDT,tt,aI(1,1),1)
        
          call allocate_vars(nveclen,neq,is) 
                                      
          !**INIT. FILE PATH FOR OUTPUTS**
          call create_file_paths
          call output_names                           

!--------------------------STIFFNESS LOOP--------------------------------------
          do jepsil = 1,jactual,1    
           
            cost(:)=0.0_wp       
                                    
            itmp = 11 - jmax/jactual         !used for 81 values of ep
            ep = 1.0_wp/10**((jepsil-1)/(itmp*1.0_wp))           
            
            !**INIT. OUTPUT FILES**
            call init_output_files(neq,ep)

!--------------------------TIMESTEP LOOP----------------------------------------
            do iDT =1,isamp,1   

              !**INITIALIZE PROBLEM INFORMATION**
              programStep='SET_INITIAL_CONDITIONS'
              call problemsub(iprob,nveclen,neq,ep,dt,tfinal,iDT,tt,aI(1,1),1)  
              dto = dt        !store time step
              t = 0.0_wp      !init. start time

              !**INIT. ERROR VECTOR**
              errvecT(:) = 0.0_wp
        
              do i = 1,ns            !initialize stage value preditor
                predvec(:,i) = uvec(:)
              enddo
!--------------------------TIME ADVANCEMENT LOOP-------------------------------
              do ktime = 1,100000000
                if(t+dt > tfinal)dt = tfinal-t+1.0e-11_wp !check if dt>tfinal
                tt=t+Ci(1)*dt
                !**STORE VALUES OF UVEC**
                uveco(:) = uvec(:)
             
                jcount = jcount + (ns-1)    !keep track of total RK stages  
!--------------------------------RK LOOP---------------------------------------
                programStep='BUILD_RHS'
                call problemsub(iprob,nveclen,neq,ep,dt,tfinal,iDT,tt,aI(1,1),1) 

                ustage(:,1) = uvec(:)

                !  ESDIRK Loop
                do L = 2,ns                        
                  tt=t+Ci(L)*dt

                  usum(:) = uveco(:)
                  do LL = 1,L-1 
                    usum(:) = usum(:)+ aI(L,LL)*resI(:,LL)+ aE(L,LL)*resE(:,LL)                               
                  enddo
              
                  if(ipred==2) predvec(:,L)=uvec(:)!previous guess as starter
                  if(ipred/=2) uvec(:) = predvec(:,L)  
!-----------------BEG NEWTON ITERATION ------------------------------------------
                  call Newton_Iteration(iprob,L,ep,dt,nveclen,&
     &                                      tt,aI(L,L),icount,k)
!-----------------END NEWTON ITERATION-------------------------------------------
                  ustage(:,L) = uvec(:)     !  Save the solution at each stage
                 
                  ! Fill in resE and resI with the converged data
                  programStep='BUILD_RHS'
                  call problemsub(iprob,nveclen,neq,ep,dt,tfinal,iDT,tt,aI(1,1),L)

                  if(ipred/=2)call Stage_Value_Predictor(ipred,L,ktime)

                  if((k > maxiter(L)).and.(ktime/=1)) maxiter(L) = k

                  stageE(L) = stageE(L) + xnorm(uvec,predvec(:,L))
                  stageI(L) = stageI(L) + 1.0_wp*k

                  if(ipred==2)call Stage_Value_Predictor(ipred,L,ktime)
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

                if (time_sol_flag .and. iDT==time_sol_iDT) then
                  call write_time_depen_sol(t,ep,nveclen,neq)
                endif
                
                t = t + dt                  !increment time
                if(t >= tfinal) exit        
              enddo                                                     
!-----------------------END TIME ADVANCEMENT LOOP------------------------------
              if (iDT_sol_flag .and. jepsil==jmax) call output_iDT_sol(iDT)           
     
              cost(iDT) = log10((ns-1)/dto)    !  ns - 1 implicit stages
              
              call output_conv_error(cost(iDT),nveclen,neq,iprob)
              
              tmpvec(:) = abs(uvec(:)-uexact(:))

              do i = 1,nvecLen
                if(tmpvec(i) <= 0.0_wp)tmpvec(i)=1.0e-15_wp
              enddo
              error(iDT,:)  = log10(tmpvec(:))
              errorP(iDT,:) = log10(errvecT(:))
              do i = 1,neq
                errorL2(iDT,i)= log10(sqrt(dot_product(tmpvec(i:nveclen:neq), &
     &                                   tmpvec(i:nveclen:neq))/(nveclen/neq)))
              enddo
            enddo  
!----------------------------END TIMESTEP LOOP---------------------------------
            if(exact_sol_flag) call write_exact_sol()
!----------------------------OUTPUTS-------------------------------------------

            call output_terminal_iteration(cost,0,ep,nveclen,neq)

            epsave(jepsil) = log10(ep)
            do i = 1,nveclen
              b1save(jepsil,i) = -b(i)
              b1Psave(jepsil,i) = -b(i+nveclen)
            enddo
            do i=1,neq
              b1L2save(jepsil,i) = -b(2*nveclen+i)
            enddo
            
          enddo      
!---------------------END STIFFNESS LOOP---------------------------------------
          !**OUTPUT CONVERGENCE VS STIFFNESS**
          call output_conv_stiff(nveclen,neq,epsave)
          
          !**OUTPUT TO TERMINAL**
          call output_terminal_final(icount,jcount,stageE,stageI,maxiter)
          
          call cpu_time(cputime2)
          write(*,*)'Total time elapsed for this case: ',cputime2-cputime1,'sec'  

!----------------------END OUTPUTS---------------------------------------------
          call deallocate_vars
        enddo                                       
!----------------------END PROBLEMS LOOP---------------------------------------
      enddo                                 
!----------------------END ALGORITHMS LOOP-------------------------------------
!----------------------END PROGRAM---------------------------------------------
      END PROGRAM test_cases
