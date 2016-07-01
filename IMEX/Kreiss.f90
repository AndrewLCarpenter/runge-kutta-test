!******************************************************************************
! Subroutine to Initialize, calculate the RHS, and calculate the Jacobian
! of the Kreiss problem (Dekker 7.5.1 pp. 214 (Kreiss' problem: Index 2)) 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************

      subroutine Kreiss(programStep,nveclen,ep,dt, &
     &                  tfinal,iDT,time,rese_vec,resi_vec,akk)

      use precision_vars
      use control_variables

      implicit none
!-----------------------VARIABLES----------------------------------------------
      integer, parameter     :: vecl=2
      integer, intent(in   ) :: programStep

      !INIT vars
      real(wp),                  intent(in   ) :: ep
      real(wp),                  intent(inout) :: dt
      integer,                   intent(inout) :: nveclen
      real(wp),                  intent(  out) :: tfinal
      integer,                   intent(in   ) :: iDT
      real(wp),                  intent(in   ) :: time

      real(wp) :: epi,lambda_p,lambda_m,sin_t,cos_t  
      real(wp), dimension(nveclen,nveclen) :: E_mat,wrk_matrix
      real(wp), dimension(nveclen)         :: wrk_vec   
      
      !RHS vars
      real(wp), dimension(vecl), intent(  out) :: rese_vec,resi_vec
      
      !Jacob vars
      real(wp), intent(in   ) :: akk
!------------------------------------------------------------------------------

      !**Pre-initialization. Get problem name and vector length**
      if (programStep==-1) then
        nvecLen = vecl
        probname='Kreiss   '   
        tol=1.0e-10_wp
        
      else if (programStep==0) then
        !**Initialize constants**
        epi =1.0_wp/ep
        cos_t=cos(time)
        sin_t=sin(time)
        
        dt = 0.25_wp/10**((iDT-1)/20.0_wp) ! timestep
        tfinal = 1.0_wp                    ! final time

        !**Exact Solution**
        lambda_p=0.5_wp*(-1.0_wp-epi) + 0.5_wp*sqrt((1.0_wp-epi)**2-4.0_wp)
        lambda_m=0.5_wp*(-1.0_wp-epi) - 0.5_wp*sqrt((1.0_wp-epi)**2-4.0_wp)
       
        E_mat(1,1)=cos(tfinal)
        E_mat(1,2)=-sin(tfinal)
        E_mat(2,1)=sin(tfinal)
        E_mat(2,2)=cos(tfinal)
        
        wrk_matrix(1,1)=ep
        wrk_matrix(1,2)=ep
        wrk_matrix(2,1)=1+lambda_p*ep
        wrk_matrix(2,2)=1+lambda_m*ep
        
        wrk_vec(1)=exp(lambda_p*tfinal)
        wrk_vec(2)=exp(lambda_m*tfinal)
        
        uexact=matmul(matmul(E_mat,wrk_matrix),wrk_vec)
        
        !**IC**
        uvec(1) = 2*ep
        uvec(2) = -1.0_wp*ep+1.0_wp
        
      !**RHS and Jacobian**
      elseif (programStep>=1) then
 
      !**Initialize constants**
      cos_t=cos(time)
      sin_t=sin(time)
      
        select case (Temporal_Splitting)

          case('IMEX') ! For IMEX schemes
            !**RHS**
            if (programStep==1 .or.programStep==2) then
              rese_vec(1) = dt*(-sin_t*sin_t*uvec(1)+sin_t*cos_t*uvec(2))
              rese_vec(2) = dt*( sin_t*cos_t*uvec(1)-cos_t*cos_t*uvec(2))
              resi_vec(1) = dt*(-cos_t*cos_t*uvec(1)-sin_t*cos_t*uvec(2))/ep
              resi_vec(2) = dt*(-sin_t*cos_t*uvec(1)-sin_t*sin_t*uvec(2))/ep
            !**Jacobian**
            elseif (programStep==3) then
              xjac(1,1) = 1.0_wp-akk*dt*(-cos_t*cos_t)/ep
              xjac(1,2) = 0.0_wp-akk*dt*(-sin_t*cos_t)/ep
              xjac(2,1) = 0.0_wp-akk*dt*(-sin_t*cos_t)/ep
              xjac(2,2) = 1.0_wp-akk*dt*(-sin_t*sin_t)/ep
            endif
            
          case('IMPLICIT') ! For fully implicit schemes
            !**RHS**
            if (programStep==1 .or.programStep==2) then
              rese_vec(:) = 0.0_wp
              resi_vec(1) =  dt*(-sin_t*sin_t*uvec(1)+sin_t*cos_t*uvec(2))+ &
     &                       dt*(-cos_t*cos_t*uvec(1)-sin_t*cos_t*uvec(2))/ep
              resi_vec(2) =  dt*( sin_t*cos_t*uvec(1)-cos_t*cos_t*uvec(2))+ &
     &                       dt*(-sin_t*cos_t*uvec(1)-sin_t*sin_t*uvec(2))/ep
            !**Jacobian**
            elseif (programStep==3) then
              xjac(1,1) = 1.0_wp-akk*dt*(-cos_t*cos_t-ep*sin_t*sin_t)/ep
              xjac(1,2) = 0.0_wp-akk*dt*(-cos_t*sin_t+ep*sin_t*cos_t)/ep
              xjac(2,1) = 0.0_wp-akk*dt*(-cos_t*sin_t+ep*cos_t*sin_t)/ep
              xjac(2,2) = 1.0_wp-akk*dt*(-sin_t*sin_t-ep*cos_t*cos_t)/ep
            endif
      
          case default ! To catch invald inputs
            write(*,*)'Invaild case entered. Enter "IMEX" or "IMPLICIT"'
            write(*,*)'Exiting'
            stop
            
        end select
        
      endif
      
      end subroutine Kreiss

      

