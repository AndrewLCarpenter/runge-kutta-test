      subroutine Kreiss(programStep,probname,nveclen,temporal_splitting,uvec,ep,uexact,dt, &
                       &  tfinal,iDT,time,resE,resI,akk,xjac)

      use precision_vars

      implicit none

      integer,  parameter                      :: vecl=2

      integer,                   intent(in   ) :: programStep
      character(80),             intent(in   ) :: temporal_splitting

      !INIT vars
      character(len=9),          intent(  out) :: probname
      real(wp), dimension(vecl), intent(inout) :: uvec
      real(wp),                  intent(in   ) :: ep
      real(wp), dimension(vecl), intent(  out) :: uexact
      real(wp),                  intent(inout) :: dt
      integer,                   intent(inout) :: nveclen
      real(wp),                  intent(  out) :: tfinal
      integer,                   intent(in   ) :: iDT
      real(wp),                  intent(in   ) :: time

      real(wp)                                 :: tmp,tmpM,tmpP
      real(wp) :: epi,lambda_p,lambda_m,a,b,sin_t,cos_t  
      real(wp), dimension(nveclen,nveclen) :: E_mat,wrk_matrix
      real(wp), dimension(nveclen)         :: wrk_vec   
      !RHS vars
      real(wp), dimension(vecl+1), intent(  out) :: resE,resI
      
      !Jacob vars
      real(wp),                       intent(in   ) :: akk
      real(wp), dimension(vecl,vecl), intent(  out) :: xjac

      epi =1/ep
      cos_t=cos(time)
      sin_t=sin(time)
      
      if (programStep==-1) then
        nvecLen = vecl
        probname='Kreiss   '   
      else if (programStep==0) then
        dt = 0.25_wp/10**((iDT-1)/20.0_wp)

        tfinal = 1.0_wp
        tmpP = exp(+tfinal   )
        tmpM = exp(-tfinal/ep)
        tmp = exp(tfinal)
        uvec(1) = 1.0_wp
        uvec(2) = 1.0_wp
        lambda_p=0.5_wp*(-1.0_wp-epi) + 0.5_wp*sqrt((1.0_wp-epi)**2-4.0_wp)
        lambda_m=0.5_wp*(-1.0_wp-epi) - 0.5_wp*sqrt((1.0_wp-epi)**2-4.0_wp)
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
        select case (Temporal_Splitting)

          case('IMEX')
                if (programStep==1 .or.programStep==2) then
              resE(1) = dt*(-sin_t*sin_t*uvec(1)-sin_t*cos_t*uvec(2))
              resE(2) = dt*(-sin_t*cos_t*uvec(1)-cos_t*cos_t*uvec(2))
              resI(1) = dt*(-cos_t*cos_t*uvec(1)-sin_t*cos_t*uvec(2))/ep
              resI(2) = dt*(-sin_t*cos_t*uvec(1)-sin_t*sin_t*uvec(2))/ep
            elseif (programStep==3) then
              xjac(1,1) = 1.0_wp-akk*dt*(-cos_t*cos_t)/ep
              xjac(1,2) = 0.0_wp-akk*dt*(-sin_t*cos_t)/ep
              xjac(2,1) = 0.0_wp-akk*dt*(-sin_t*cos_t)/ep
              xjac(2,2) = 1.0_wp-akk*dt*(-sin_t*sin_t)/ep
            endif
          case('IMPLICIT')
            if (programStep==1 .or.programStep==2) then
            !  resE(1) = 0.0_wp
            !  resE(2) = 0.0_wp
           !   resI(1) = dt
           !   resI(2) = 
            elseif (programStep==3) then
            !  xjac(1,1) = 1.0_wp-akk*dt*(0.0_wp)
            !  xjac(1,2) = 0.0_wp-akk*dt/ep
           !   xjac(2,1) = 0.0_wp-akk*dt*( 0.0_wp)
           !   xjac(2,2) = 1.0_wp-akk*dt*(0.0_wp) !possible error here
            endif
        end select
      endif
      
      return
      end subroutine

      

