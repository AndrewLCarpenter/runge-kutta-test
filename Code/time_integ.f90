!*****************************************************************************
! Module to perform integration over one timestep 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *ONTAINS VARIABLES USED IN THE PROGRAM
! Stage_value_module 
! Runge_Kutta.F90    
! NEWTON.F90         
! PROBLEMSUB.F90            *DEFINES WHICH PROBLEM IS RELATED TO USER INPUT
!******************************************************************************
      module time_integ_mod

      use precision_vars, only: wp

      implicit none; save

      public :: DIRK, DIRK_OTD, FIRK
      private

!--------------------------------CONSTANTS-------------------------------------         
      real(wp), parameter :: sqrt2 = sqrt(2.0_wp)


      contains

!========================================================================================

      subroutine DIRK(iprob,nveclen,neq,iDT,ktime,icount,ep,t,dt,     &
                     & stageE,stageI,maxiter)

        use precision_vars,     only: wp
        use control_variables,  only: uvec,uveco,ustage,resE,resI,usum,    &
                                    & programstep,ipred,predvec,        &
                                    & errvec
        use Stage_value_module, only: Stage_Value_Predictor,xnorm
        use runge_kutta,        only: aE,aI,bE,bI,bEH,bIH,cI,ns
        use Newton,             only: Newton_Iteration
        use problemsub_mod,     only: problemsub

        integer,  intent(in)    :: iprob, iDT, ktime
        integer,  intent(inout) :: nveclen, neq, icount
        real(wp), intent(inout) :: ep, t, dt 

        integer,  dimension(:), intent(inout)  :: maxiter
        real(wp), dimension(:), intent(inout)  :: stageE, stageI

!-----------------------------VARIABLES----------------------------------------  

        !internal variables

        real(wp)                :: tfinal, tt
        integer                 :: k,L,LL

!------------------------------------------------------------------------------   
           continue
!------------------------------------------------------------------------------   

           tt = t + cI(1)*dt                        ! intermediate time (stage value)

           uveco(:)    = uvec(:)                    !**STORE VALUES OF UVEC**
           ustage(:,1) = uvec(:)                    ! Save the solution at first stage

           programStep='BUILD_RHS'
           call problemsub(iprob,nveclen,neq,ep,dt,tfinal,iDT,tt,aI(1,1),1) 

           do L = 2,ns                         !===== Begin: A_{k,j} portion of ESDIRK LOOP

             tt = t + cI(L)*dt                      ! intermediate time

             usum(:) = uveco(:)                     ! sum explicit and implicit rhs vector
             do LL = 1,L-1 
               usum(:) = usum(:)+ aI(L,LL)*resI(:,LL)+ aE(L,LL)*resE(:,LL)                               
             enddo

             if(ipred==2) predvec(:,L)=uvec(:)      !previous guess as starter
             if(ipred/=2) uvec(:) = predvec(:,L)  
                                                    ! BEG NEWTON ITERATION ----------------
             call Newton_Iteration(iprob,L,ep,dt,nveclen,tt,aI(L,L),icount,k)

                                                    ! END NEWTON ITERATION-----------------
             ustage(:,L) = uvec(:)                  !  Save the solution at each stage
                 
             programStep='BUILD_RHS'                ! Fill resE and resI with the converged data
             call problemsub(iprob,nveclen,neq,ep,dt,tfinal,iDT,tt,aI(L,L),L)

             if(ipred/=2)call Stage_Value_Predictor(ipred,L,ktime)

             if((k > maxiter(L)).and.(ktime/=1)) maxiter(L) = k

             stageE(L) = stageE(L) + xnorm(uvec,predvec(:,L))
             stageI(L) = stageI(L) + 1.0_wp*k

             if(ipred==2)call Stage_Value_Predictor(ipred,L,ktime)

           enddo                               !===== End: A_{k,j} portion of ESDIRK LOOP
     
             uvec(:) = uveco(:)                      !----Final b_{j} Sum of ESDIRK loop
           errvec(:) = 0.0_wp
           do LL = 1,ns 
               uvec(:) =   uvec(:) + (bE(LL)        )*resE(:,LL) &
                                   + (bI(LL)        )*resI(:,LL)

             errvec(:) = errvec(:) + (bE(LL)-bEH(LL))*resE(:,LL) &
                              &    + (bI(LL)-bIH(LL))*resI(:,LL)
           enddo

         end subroutine DIRK

!========================================================================================

      subroutine DIRK_OTD(iprob,nveclen,neq,iDT,ep,t,dt)

        use precision_vars,     only: wp
        use control_variables,  only: ustage,programstep,OTDN,Ovec,Oveco,Osum,err_OTD, &
                                    & resE_Tens_OTD, resI_Tens_OTD
        use runge_kutta,        only: aE,aI,bE,bI,bEH,bIH,cI,ns
        use Newton,             only: Newton_Iteration_OTD
        use problemsub_mod,     only: problemsub

        integer,  intent(in)    :: iprob, iDT
        integer,  intent(inout) :: nveclen, neq
        real(wp), intent(inout) :: ep, t, dt 

!-----------------------------VARIABLES----------------------------------------  

        !internal variables

        real(wp)                :: tfinal, tt
        integer                 :: k,L,LL

!------------------------------------------------------------------------------   
           continue
!------------------------------------------------------------------------------   

           tt = t + cI(1)*dt                        ! intermediate time (stage value)

           Oveco(:,:) = Ovec(:,:)                   !**STORE VALUES OF Ovec**

           resE_Tens_OTD(:,:,:) = 0.0_wp ; resI_Tens_OTD(:,:,:) = 0.0_wp ; 

           programStep='BUILD_RHS_OTD'
           call problemsub(iprob,nveclen,neq,ep,dt,tfinal,iDT,tt,aI(1,1),1) 

           do L = 2,ns                         !===== Begin: A_{k,j} portion of ESDIRK LOOP

             tt = t + cI(L)*dt                      ! intermediate time

             Osum(:,:) = Oveco(:,:)                    
             do LL = 1,L-1 
               
               Osum(:,:) = Osum(:,:) + aE(L,LL)*resE_Tens_OTD(:,:,LL) &
                                     + aI(L,LL)*resI_Tens_OTD(:,:,LL)
             enddo
                                                    ! BEG NEWTON ITERATION ----------------
             call Newton_Iteration_OTD(iprob,L,ep,dt,nveclen,tt,aI(L,L))
                                                    ! END NEWTON ITERATION-----------------

             programStep='BUILD_RHS_OTD'                ! Fill resE and resI with the converged data
             call problemsub(iprob,nveclen,neq,ep,dt,tfinal,iDT,tt,aI(L,L),L)

           enddo                               !===== End: A_{k,j} portion of ESDIRK LOOP
     
              Ovec(:,:) = Oveco(:,:)                    !----Final b_{j} Sum of ESDIRK loop
           err_OTD(:,:) = 0.0_wp
           do LL = 1,ns 
                Ovec(:,:) =    Ovec(:,:) + (bE(LL)        )*resE_Tens_OTD(:,:,LL) &
                        &                + (bI(LL)        )*resI_Tens_OTD(:,:,LL)

             err_OTD(:,:) = err_OTD(:,:) + (bE(LL)-bEH(LL))*resE_Tens_OTD(:,:,LL) &
                        &                + (bI(LL)-bIH(LL))*resI_Tens_OTD(:,:,LL)
           enddo

         end subroutine DIRK_OTD

!========================================================================================

      subroutine FIRK(iprob,nveclen,neq,iDT,ktime,icount,ep,t,dt,     &
                     & stageE,stageI,maxiter)

        use precision_vars,     only: wp
        use control_variables,  only: uvec,uveco,resI,programstep,errvec, &
                                      uvecS,uvecoS
        use Stage_value_module, only: Stage_Value_Predictor,xnorm
        use runge_kutta,        only: aE,aI,bE,bI,bEH,bIH,cI,ns
        use Newton,             only: Newton_Iteration_FIRK
        use problemsub_mod,     only: problemsub

        integer,  intent(in)    :: iprob, iDT, ktime
        integer,  intent(inout) :: nveclen, neq, icount
        real(wp), intent(inout) :: ep, t, dt 

        integer,  dimension(:), intent(inout)  :: maxiter
        real(wp), dimension(:), intent(inout)  :: stageE, stageI

!-----------------------------VARIABLES----------------------------------------  

        !internal variables

        real(wp)                :: tfinal, tt
        integer                 :: k,L,LL,iL,iH, nvecS

!------------------------------------------------------------------------------   
           continue
!------------------------------------------------------------------------------   

           nvecS = nveclen * ns
           uveco(:)    = uvec(:)                    !**STORE VALUES OF UVEC**
           
           do L = 1,ns                              ! Initialize stage values
             iL = (L-1)*nveclen + 1
             iH = (L-1)*nveclen + nveclen
             uvecS(iL:iH) = uvec(:)
           enddo
           uvecoS(:) = uvecS(:)

           call Newton_Iteration_FIRK(iprob,ep,dt,nveclen,nvecS,t,icount,k)

             uvec(:) = uveco(:)                     !----Final b_{j} Sum of ESDIRK loop
           errvec(:) = 0.0_wp                       !----Error estimation vector
           do LL = 1,ns 
               uvec(:) =   uvec(:) + (bI(LL)        )*resI(:,LL)
             errvec(:) = errvec(:) + (bI(LL)-bIH(LL))*resI(:,LL)
           enddo

         end subroutine FIRK

!========================================================================================

         end module time_integ_mod
