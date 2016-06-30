      subroutine Kreiss(programStep,nveclen,ep,dt, &
                       &  tfinal,iDT,time,rese_vec,resi_vec,akk)

      use precision_vars
      use control_variables

      implicit none

      integer,  parameter                      :: vecl=2

      integer,                   intent(in   ) :: programStep

      !INIT vars
      !real(wp), dimension(vecl), intent(inout) :: uvec
      real(wp),                  intent(in   ) :: ep
     ! real(wp), dimension(vecl), intent(  out) :: uexact
      real(wp),                  intent(inout) :: dt
      integer,                   intent(inout) :: nveclen
      real(wp),                  intent(  out) :: tfinal
      integer,                   intent(in   ) :: iDT
      real(wp),                  intent(in   ) :: time

      real(wp) :: epi,lambda_p,lambda_m,sin_t,cos_t  
      real(wp), dimension(nveclen,nveclen) :: E_mat,wrk_matrix
      real(wp), dimension(nveclen)         :: wrk_vec   
      !RHS vars
      real(wp), dimension(vecl+1), intent(  out) :: rese_vec,resi_vec
      
      !Jacob vars
      real(wp),                       intent(in   ) :: akk
      !real(wp), dimension(vecl,vecl), intent(  out) :: xjac
      !print*,'in Kreiss'

      
      if (programStep==-1) then
       ! print*,'in kreiss step=-1'
        nvecLen = vecl
        probname='Kreiss   '   
        tol=1.0e-10_wp
      else if (programStep==0) then
            epi =1/ep
      cos_t=cos(time)
      sin_t=sin(time)
       ! print*,'in kreiss step=0'
        dt = 0.25_wp/10**((iDT-1)/20.0_wp)

        tfinal = 1.0_wp


        lambda_p=0.5_wp*(-1.0_wp-epi) + 0.5_wp*sqrt((1.0_wp-epi)**2-4.0_wp)
        lambda_m=0.5_wp*(-1.0_wp-epi) - 0.5_wp*sqrt((1.0_wp-epi)**2-4.0_wp)
        uvec(1) = 2*ep
        uvec(2) = -1.0_wp*ep+1.0_wp
        !uexact(1) = uvec(1)*tmpP + uvec(2)*(tmpP - tmpM)/(1.0_wp + ep)
        !uexact(2) = uvec(2)*tmpM
       ! a = ep*(exp(lambda_p*tfinal)+exp(lambda_m*tfinal))
      ! b = (1.0_wp+ep*lambda_p)*exp(lambda_p*tfinal)+(1.0_wp+ep*lambda_m)*exp(lambda_m*tfinal)
       ! uexact(1) = cos(tfinal)*a-sin(tfinal)*b
       !uexact(2) = sin(tfinal)*a+cos(tfinal)*b
        !uexact(1) = -exp(tfinal)*sin(tfinal)
        !uexact(2) = exp(tfinal)*cos(tfinal)
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

      elseif (programStep>=1 .and. programStep<=3) then
            epi =1/ep
      cos_t=cos(time)
      sin_t=sin(time)
        select case (Temporal_Splitting)

          case('IMEX')
                if (programStep==1 .or.programStep==2) then
              rese_vec(1) = dt*(-sin_t*sin_t*uvec(1)+sin_t*cos_t*uvec(2))
              rese_vec(2) = dt*( sin_t*cos_t*uvec(1)-cos_t*cos_t*uvec(2))
              resi_vec(1) = dt*(-cos_t*cos_t*uvec(1)-sin_t*cos_t*uvec(2))/ep
              resi_vec(2) = dt*(-sin_t*cos_t*uvec(1)-sin_t*sin_t*uvec(2))/ep
            elseif (programStep==3) then
              xjac(1,1) = 1.0_wp-akk*dt*(-cos_t*cos_t)/ep
              xjac(1,2) = 0.0_wp-akk*dt*(-sin_t*cos_t)/ep
              xjac(2,1) = 0.0_wp-akk*dt*(-sin_t*cos_t)/ep
              xjac(2,2) = 1.0_wp-akk*dt*(-sin_t*sin_t)/ep
            endif
          case('IMPLICIT')
            if (programStep==1 .or.programStep==2) then
              rese_vec(1) = 0.0_wp
              rese_vec(2) = 0.0_wp
              resi_vec(1) =  dt*(-sin_t*sin_t*uvec(1)+sin_t*cos_t*uvec(2))+ &
     &                   dt*(-cos_t*cos_t*uvec(1)-sin_t*cos_t*uvec(2))/ep
              resi_vec(2) =  dt*( sin_t*cos_t*uvec(1)-cos_t*cos_t*uvec(2))+ &
     &                   dt*(-sin_t*cos_t*uvec(1)-sin_t*sin_t*uvec(2))/ep

            elseif (programStep==3) then
              xjac(1,1) = 1.0_wp-akk*dt*(-cos_t*cos_t-ep*sin_t*sin_t)/ep
              xjac(1,2) = 0.0_wp-akk*dt*(-cos_t*sin_t+ep*sin_t*cos_t)/ep
              xjac(2,1) = 0.0_wp-akk*dt*(-cos_t*sin_t+ep*cos_t*sin_t)/ep
              xjac(2,2) = 1.0_wp-akk*dt*(-sin_t*sin_t-ep*cos_t*cos_t)/ep
            endif
        end select
      endif
      
      return
      end subroutine

      

