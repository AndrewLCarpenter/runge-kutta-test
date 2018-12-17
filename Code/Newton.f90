!*****************************************************************************
! Module to perform Newton Iteration in one of three modes: QR decomp, line 
! search, or normal Newton Iteration. 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! QR_MODULE.F90             *CONTAINS QR ROUTINES
! CONTROL_VARIABLES.F90     *ONTAINS VARIABLES USED IN THE PROGRAM
! ILUT_MODULE.F90           *CONTAINS ROUTINES TO PERFORM LU DECOMPOSITION
! JACOBIAN_CSR_MOD.F90      *ALLOCATES JACOBIAN CSR VARIABLES
! PROBLEMSUB.F90            *DEFINES WHICH PROBLEM IS RELATED TO USER INPUT
!******************************************************************************
      module Newton
            
      use precision_vars, only:wp

      implicit none; save
      
      public :: Newton_Iteration
      private      
      
      contains    
!==============================================================================
!******************************************************************************
! Subroutine to perform Newton iteration 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From control_variables:
!   temporal_splitting -> string used in choose_RHS_type and choose_Jac_type,         character(len=80),                       not modified
!   Jac_case           -> string used in to determine CSR Jacobian or dense Jacobian, character(len=6),                        not modified
!   uvec               -> Array containing variables,                                 real(wp), dimension(u-vector length),        modified
!   usum               -> Array containing summation from test_cases,                 real(wp), dimension(u-vector length),    not modified
!   xjac               -> matrix containing dense Jacobian,                           real(wp), dimension(u-vector length**2), not modified
!
! MODULE ROUTINES:
! Build_Jac        -> Calls problemsub to build the problem's jacobian
! LU_Solver        -> Performs LU decomposition and solving
! newt_line_search -> Performs newton iteration with line search
! Convert_to_CSR   -> converts dense matrix to CSR
! QR_decomp        -> Performs QR decomposition and solving
! Nonlinear_Residual    -> Calls problemsub to build RHS and RNewton
! check_exit       -> checks whether newton iteration 
! Mat_invert       -> inverts dense matrix up to rank 4
! 
!******************************************************************************
! INPUTS:
! iprob   -> defines problem number,                                    integer
! L       -> stage number,                                              integer
! ep      -> Stiffness epsilon value,                                   real(wp)
! dt      -> timestep: created and then used later in the problem case, real(wp)
! nveclen -> total length of u-vector used for problem                  integer
! time    -> current solution time,                                     integer
! aI      -> Diagonal term from RK scheme,                              real(wp)
! INOUTS:
! icount  -> counter used to determine average iterations,              integer
! OUTPUTS:
! k       -> total number of newton iterations used for max iteration,  integer
!******************************************************************************
      subroutine Newton_Iteration(iprob,L,ep,dt,nveclen,time,aI,icount,k)
      
      use control_variables, only: temporal_splitting,jac_case,uvec,usum,xjac

!     integer, parameter :: iter_max=1500
      integer, parameter :: iter_max=100
      logical, parameter :: Line_search=.false. !set to TRUE for line search
      logical, parameter :: QR = .true. !set to TRUE for QR factorization

      integer,  intent(in   ) :: iprob,L
      real(wp), intent(in   ) :: ep,dt
      integer,  intent(in   ) :: nveclen
      real(wp), intent(in   ) :: time,aI
      integer,  intent(inout) :: icount
      integer,  intent(  out) :: k
      
      integer                              :: i,j  !Do loop variables
      real(wp), dimension(nveclen)         :: Rnewton,uveciter!Newton 
      real(wp), dimension(nveclen,nveclen) :: xjacinv !Jacobian

!------------------------------------------------------------------------------
      temporal_select: select case(temporal_splitting)
        case('EXPLICIT')
          do k = 1,iter_max
            icount = icount + 1
            uveciter(:) = uvec(:)
            uvec(:)=usum(:)
            if(check_exit(uveciter,iter_max,k)) exit            
          enddo   
        case default !IMEX or IMPLICIT
          line: if (Line_search) then
            call newt_line_search(iprob,L,ep,dt,nveclen,time,aI,icount,k)    
!           call Inexact_Newton_Dogleg(iprob,L,ep,dt,nveclen,time,aI,icount,k)
          else
            Jac_select: select case(Jac_case)
              case('SPARSE')
                do k = 1,iter_max
                  icount = icount + 1
                  uveciter(:) = uvec(:) !store old uvec

                  Rnewton=Nonlinear_Residual(ep,dt,time,aI,iprob,L)
                  call Build_Jac(ep,dt,time,aI,iprob,L)
                
                  uvec(:)=uvec(:)-LU_solver(Rnewton)
                  
                  if(check_exit(uveciter,iter_max,k)) exit
                              
                enddo  
              case('DENSE')
                do k = 1,iter_max
                  icount = icount + 1
                  uveciter(:) = uvec(:) !store old uvec   
                                
                  Rnewton=Nonlinear_Residual(ep,dt,time,aI,iprob,L)
!                 write(*,*)k,sqrt(dot_product(rnewton,rnewton))
                  call Build_Jac(ep,dt,time,aI,iprob,L)
                    
                  if (nveclen>4 .and. .not. QR) then !No explicit inverse
                    call Convert_to_CSR(xjac) !convert dense xjac to csr
                    uvec(:)=uvec(:)-LU_solver(Rnewton)
                      
                  elseif (nveclen>4 .and. QR) then !QR factorization
!                   write(*,*)'xjac'
!                   do i = 1,nvecLen
!                      write(*,*)i,(j,xjac(i,j),j=1,nvecLen)
!                   enddo
                    call QR_decomp(xjac,nveclen,rnewton) 
                       
                  elseif (nveclen<=4) then !Explicit inverse
                    xjacinv=Mat_invert(xjac)                  
                    do i = 1,nvecLen
                      do j = 1,nvecLen
                        uvec(i) = uvec(i) - xjacinv(i,j)*Rnewton(j)     !u^n+1=u^n-J^-1*F
                      enddo
                    enddo  
                  endif
          
                  if(check_exit(uveciter,iter_max,k)) exit
                    
                enddo
!               write(*,*)k,sqrt(dot_product(rnewton,rnewton))
            end select Jac_select 
          endif line
      end select temporal_select
      end subroutine Newton_Iteration
      
!==============================================================================
!******************************************************************************
! Function to check exit conditions for Newton iteration 
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From control_variables:
!   tol  -> Newton iteration exit tolerance, real(wp),                                 set
!   uvec -> Array containing variables,      real(wp), dimension(u-vector length), not modified
!
!******************************************************************************
! INPUTS:
! uveciter -> u vector iteration value,                           real(wp), dimension(u-vector length)
! k        -> newton iteration counter, used for error reporting, integer
! OUTPUTS:
! check_exit -> logical output function for if statement in newton iteration, logical
!******************************************************************************
      function check_exit(uveciter,iter_max,k)
      
      use control_variables, only: tol,uvec
      
      real(wp), dimension(:), intent(in) :: uveciter
      integer,                intent(in) :: iter_max,k
      real(wp)                           :: tmp
      logical                            :: check_exit
      
      check_exit=.false.
      tmp = sum(abs(uvec(:)-uveciter(:))) !check accuracy of zeros         
      if (k>=iter_max) write(*,*)'tmp',tmp,'k',k
      if (tmp/=tmp) then
        print*,'stopping NaN k=',k
        stop
      endif      
      if(tmp < tol) check_exit=.true. 
      
      end function check_exit
!==============================================================================
!******************************************************************************
! Function to create Rnewton from the RHS
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
! PROBLEMSUB.F90            *DEFINES WHICH PROBLEM IS RELATED TO USER INPUT
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From control_variables:
!   uvec         -> Array containing variables,                 real(wp), dimension(u-vector length),              not modified
!   usum         -> Array containing summation from test_cases, real(wp), dimension(u-vector length),              not modified
!   resI         -> Implicit RHS vector,                        real(wp), dimension(u-vector length, num. stages), not modified
!   programstep  -> string used in program_step_select,         character(len=80),                                 not modified 
! From problemsub_mod:
!   problemsub   -> Subroutine to set problem information, determined by program step
!******************************************************************************
! INPUTS:
! ep      -> Stiffness epsilon value,                                   real(wp)
! dt      -> timestep: created and then used later in the problem case, real(wp)
! time    -> current solution time,                                     integer
! aI      -> Diagonal term from RK scheme,                              real(wp)
! iprob   -> defines problem number,                                    integer
! L       -> stage number,                                              integer
! OUTPUTS: 
! Nonlinear_Residual -> array containing modified RHS for newton iteration,  real(wp), dimension(u-vector length)
!******************************************************************************    
      function Nonlinear_Residual(ep,dt,time,aI,iprob,L)
      
      use control_variables, only: uvec,usum,resI,programstep
      use problemsub_mod,    only: problemsub
            
      real(wp), dimension(size(uvec))       :: Nonlinear_Residual
      real(wp),               intent(in   ) :: ep,dt,time,aI
      integer,                intent(in   ) :: iprob,L

      integer  :: iDT,nveclen,neq
      real(wp) :: tfinal,dt_in
      
      dt_in=dt
      
      programStep='BUILD_RHS'
      call problemsub(iprob,nveclen,neq,ep,dt_in,tfinal,iDT,time,aI,L)
      Nonlinear_Residual(:) = uvec(:)-aI*resI(:,L)-usum(:)

      end function Nonlinear_Residual
!==============================================================================
!******************************************************************************
! Subroutine to get Jacobian for newton iteration
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
! PROBLEMSUB.F90            *DEFINES WHICH PROBLEM IS RELATED TO USER INPUT
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From control_variables:
!   programstep  -> string used in program_step_select, character(len=80), not modified 
! From problemsub_mod:
!   problemsub   -> Subroutine to set problem information, determined by program step
!******************************************************************************
! INPUTS:
! ep      -> Stiffness epsilon value,                                   real(wp)
! dt      -> timestep: created and then used later in the problem case, real(wp)
! time    -> current solution time,                                     integer
! aI      -> Diagonal term from RK scheme,                              real(wp)
! iprob   -> defines problem number,                                    integer
! L       -> stage number,                                              integer
!******************************************************************************        
      subroutine Build_Jac(ep,dt,time,aI,iprob,L)
      
      use control_variables, only: programstep
      use problemsub_mod,    only: problemsub     
       
      real(wp), intent(in) :: ep,dt,time,aI
      integer,  intent(in) :: iprob,L
    
      integer  :: nveclen,iDT,neq
      real(wp) :: tfinal,dt_in
      
      dt_in=dt
  
      programStep='BUILD_JACOBIAN'
      call problemsub(iprob,nveclen,neq,ep,dt_in,tfinal,iDT,time,aI,L)
      
      end subroutine Build_Jac
      
!==============================================================================
!******************************************************************************
! Subroutine to perform LU decomposition and solving
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! UNARY_MOD.F90             *PERFORMS SPARSE MATRIX OPERATIONS
! JACOBIAN_CSR_MOD.F90      *ALLOCATE AND STORE CSR JACOBIAN VARIABLES
! ILUT_MODULE.F90           *CONTAINS ROUTINES TO PERFORM LU DECOMPOSITION
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From Jacobian_CSR_Mod:
!   iaJac  -> ia matrix for global storage of Jacobian, integer,  dimension(u-vector length + 1),                                  not modified
!   jaJac  -> ja matrix for global storage of Jacobian, integer,  dimension(dependant on temporal_splitting, see problem routine), not modified
!    aJac  ->  a matrix for global storage of Jacobian, real(wp), dimension(dependant on temporal_splitting, see problem routine), not modified
!    jUJac -> matrix for storage of LU Jacobian,        integer,  dimension(u-vector length),                                          set & modified
!   jLUJac -> matrix for storage of LU Jacobian,        integer,  dimension(size of jaJac*4),                                          set & modified
!   aLUJac -> matrix for storage of LU Jacobian,        real(wp), dimension(size of aJac*4),                                           set & modified
! From ilut_module:
!   lusol -> Performs LU solving
!   ilutp -> Performs Preconditioning
!
!******************************************************************************
! INPUTS:
! Rnewton -> array containing modified RHS for newton iteration,  real(wp), dimension(u-vector length)
!******************************************************************************
      function LU_solver(Rnewton)
        
      use Jacobian_CSR_Mod,  only: iaJac,jaJac,aJac,aLUJac,jLUJac,jUJac
      use ilut_module,       only: lusol,ilutp
      use matvec_module,     only: amux
      
      real(wp), dimension(:), intent(in)          :: Rnewton
      integer                              :: nveclen,ierr=0
      real(wp), dimension(size(Rnewton))   :: LU_solver
      real(wp), dimension(size(Rnewton))   :: w  
      integer,  dimension(2*size(Rnewton)) :: jw,iperm

      nveclen=size(Rnewton)   

      call ilutp(nveclen,aJac,jaJac,iaJac,nveclen,1e-13_wp,0.1_wp,nveclen, &
     &           aLUJac,jLUJac,jUJac,size(alujac),w,jw,iperm,ierr)

      call lusol(nveclen,Rnewton,LU_solver,aLUJac,jLUJac,jUJac)     

      if (ierr/=0) then 
        print*,'Error in LU_solver!'
        print*,'ierr=',ierr
        stop
      endif
      end function LU_solver 

!==============================================================================
!******************************************************************************
! Subroutine to perform Newton iteration with line search
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From control_variables:
!   uvec -> Array containing variables,       real(wp), dimension(u-vector length),        modified
!   xjac -> matrix containing dense Jacobian, real(wp), dimension(u-vector length**2), not modified
!
! MODULE ROUTINES:
! Build_Jac     -> Calls problemsub to build the problem's jacobian
! Nonlinear_Residual -> Calls problemsub to build RHS and RNewton
! check_exit    -> checks whether newton iteration 
! Mat_invert    -> inverts dense matrix up to rank 4
! 
!******************************************************************************
! INPUTS:
! iprob   -> defines problem number,                                    integer
! L       -> stage number,                                              integer
! ep      -> Stiffness epsilon value,                                   real(wp)
! dt      -> timestep: created and then used later in the problem case, real(wp)
! nveclen -> total length of u-vector used for problem                  integer
! time    -> current solution time,                                     integer
! aI      -> Diagonal term from RK scheme,                              real(wp)
! INOUTS:
! icount  -> counter used to determine average iterations,              integer
! OUTPUTS:
! k       -> total number of newton iterations used for max iteration,  integer
!******************************************************************************
      subroutine newt_line_search(iprob,L,ep,dt,nveclen,time,aI,icount,k)
      
      use control_variables, only: uvec,xjac,Jac_case
      
      integer,  intent(in   ) :: iprob,L
      real(wp), intent(in   ) :: ep,dt
      integer,  intent(in   ) :: nveclen
      real(wp), intent(in   ) :: time,aI
      integer,  intent(inout) :: icount,k
      
      integer                              :: j  !Do loop variables
      integer                              :: ierr
      real(wp), dimension(nveclen)         :: Rnewton,ustor !Newton 
      real(wp), dimension(nveclen,nveclen) :: xjacinv !Jacobianxjac,

      real(wp) :: al,rnorm,rnormt
      real(wp), dimension(nveclen) :: dxi

!     real(wp) :: delta_k, delta_km1
!     real(wp) ::   eta_k,   eta_km1

      
      do k = 1,150
        icount = icount + 1
   
        Rnewton=Nonlinear_Residual(ep,dt,time,aI,iprob,L)
    
        rnorm = sqrt(dot_product(Rnewton(:),Rnewton(:)))

        call Build_Jac(ep,dt,time,aI,iprob,L)

        select case(Jac_case)
        case('DENSE')
        xjacinv=Mat_invert(xjac)

        dxi = MatMul(xjacinv,Rnewton(:))
        case('SPARSE')
        dxi = LU_solver(Rnewton)
        end select
        
        al = 1.0_wp
        do j = 1,10    !   under-relax the value of the parameter alpha
        
          ustor(:)=uvec(:)!becuause uvec is global, it needs to be temp
                          !stored so that calculations are done correctly
          uvec(:) = uvec(:) - al*dxi
             
          Rnewton=Nonlinear_Residual(ep,dt,time,aI,iprob,L)
    
          rnormt = sqrt(dot_product(Rnewton(:),Rnewton(:)))

          uvec(:)=ustor(:)!reset uvec from temp storage

           if((rnormt >= rnorm) .and. (rnorm >= 1.0e-4_wp)) then
             al = 0.5_wp * al
           else
             exit
           endif
        enddo

        uvec(:) = uvec(:) - al*dxi

        rnorm = rnormt
        if (k >= 140) print*,'L',L,'k',k,'tmp',rnorm!,'j',j
        if(rnorm <= 1.0e-9_wp) then
          ierr = 0
          return
        endif
      enddo      
      end subroutine newt_line_search

!==============================================================================

      subroutine Inexact_Newton_Dogleg(iprob,L,ep,dt,nveclen,time,aI,icount,k)
      
      use control_variables, only: uvec,xjac,Jac_case
      use Jacobian_CSR_Mod,  only: iaJac,jaJac,aJac
      use matvec_module,     only: amux, atmux
      
      integer,  intent(in   ) :: iprob,L
      real(wp), intent(in   ) :: ep,dt
      integer,  intent(in   ) :: nveclen
      real(wp), intent(in   ) :: time,aI
      integer,  intent(inout) :: icount,k
      
      real(wp), dimension(nveclen,nveclen) :: xjacinv !Jacobianxjac,

      real(wp), dimension(nveclen) :: ustor
      real(wp), dimension(nveclen) :: S_k, S_InNewt, S_Cauchy, SteepD_k
      real(wp), dimension(nveclen) :: J_SteepD_k, J_S_Cauchy, J_S_k
 
      real(wp), dimension(nveclen) :: R_NonLin_k,        R_NonLin_predct_k
      real(wp), dimension(nveclen) :: R_Cauchy_predct_k, R_InNewt_predct_k

      real(wp)                     :: S_Cauchy_L2_k,        S_InNewt_L2
      real(wp)                     :: S_L2_k
      real(wp)                     :: R_NonLin_L2_k,        R_NonLin_L2_km1
      real(wp)                     :: R_NonLin_predct_L2_k, R_NonLin_predct_L2_km1
      real(wp)                     :: R_Cauchy_predct_L2_k
      real(wp)                     :: SteepD_L2_k, J_SteepD_L2_k

      real(wp)                     :: delta_k
      real(wp)                     ::   eta_k
      real(wp)                     :: gamma_k
      real(wp)                     :: ratio, small
      real(wp), parameter          :: T1 = -10000000.0_wp

      integer                      :: ierr
      
      R_NonLin_L2_k          = T1 ; R_NonLin_L2_km1        = T1 ;
      R_NonLin_predct_L2_k   = T1 ; R_NonLin_predct_L2_km1 = T1 ;

      S_Cauchy               = T1 ; J_S_Cauchy             = T1 ;
      S_InNewt               = T1 ; S_k                    = T1 ;
      SteepD_k               = T1 ; J_SteepD_k             = T1 ;
      gamma_k = 1.0_wp            ;  small = 1.0e-40_wp
      eta_k   = 0.01_wp           ;  Delta_k = 0.1 ;

      !     F(x_k)  &  ||F(x_k)||
      R_NonLin_k = Nonlinear_Residual(ep,dt,time,aI,iprob,L)
      R_NonLin_L2_k  = sqrt(dot_product(R_NonLin_k(:),R_NonLin_k(:))+small)

      do k = 1,25

      write(*,*)'k,R_NonLin_L2_k',k,R_NonLin_L2_k
      write(*,*)'k,Delta_k',k,Delta_k
      write(*,*)'k,Eta_k',k,Eta_k

        ustor(:)=uvec(:)  ! uvec is global and needs to be temp

      !     {\partial F}{\partial (x_k)}
        call Build_Jac(ep,dt,time,aI,iprob,L)

         
      !  First build Cauchy point: it does not require matrix inverse
        select case(Jac_case)
          case('DENSE')
              SteepD_k = + MatMul(Transpose(xjac),R_NonLin_k)
            J_SteepD_k = + MatMul(xjac,SteepD_k)
          case('SPARSE')
            call atmux(nveclen,R_NonLin_k,SteepD_k,aJac,jaJac,iaJac)
            call amux (nveclen,SteepD_k,J_SteepD_k,aJac,jaJac,iaJac)
        end select

          SteepD_L2_k = sqrt(dot_product(  SteepD_k,  SteepD_k)+small)
        J_SteepD_L2_k = sqrt(dot_product(J_SteepD_k,J_SteepD_k)+small)

        ratio = SteepD_L2_k**2/J_SteepD_L2_k**2 
        S_Cauchy(:)      = + ratio * SteepD_k(:)
        S_Cauchy_L2_k    = sqrt(dot_product(S_Cauchy,S_Cauchy)+small)

        if (k<=2) then
           delta_k = S_Cauchy_L2_k
           S_k(:)  = S_Cauchy(:)
!          write(*,*)'norm S_k: first step',l,k,sqrt(dot_product(S_k,S_k))
        elseif ( S_Cauchy_L2_k >= Delta_k) then
            ratio = Delta_k / S_Cauchy_L2_k
           S_k(:) = ratio * S_Cauchy(:)
!          write(*,*)'Cauchy path violates trust delta_k: S_Cauchy_L2_k', S_Cauchy_L2_k
           write(*,*)'norm S_k: reduced Cauchy',l,k,sqrt(dot_product(S_k,S_k))
        else
          select case(Jac_case)
            case('DENSE')
              J_S_Cauchy = MatMul(xjac,S_Cauchy)
            case('SPARSE')
              call amux(nveclen,S_Cauchy,J_S_Cauchy,aJac,jaJac,iaJac)
          end select
          R_Cauchy_predct_k(:) = R_NonLin_k(:) + J_S_Cauchy(:)
          R_Cauchy_predct_L2_k = sqrt(dot_product(R_Cauchy_predct_k, &
                                                  R_Cauchy_predct_k)+small)

          if(R_Cauchy_predct_L2_k <= Eta_k * R_NonLin_L2_k ) then
             S_k(:) = S_Cauchy(:)
!            write(*,*)'Cauchy path satisfies reduction tolerance'
             write(*,*)'norm S_k: Cauchy reduction',l,k,sqrt(dot_product(S_k,S_k))
          else
            select case(Jac_case)
              case('DENSE')
                xjacinv=Mat_invert(xjac) ; S_InNewt = MatMul(xjacinv,R_NonLin_k)
              case('SPARSE')
                S_InNewt = LU_solver(R_NonLin_k)
            end select

            S_InNewt_L2 = sqrt(dot_product(S_InNewt,S_InNewt)+small)

            if(S_InNewt_L2 <= Delta_k) then
              S_k(:) = S_InNewt(:)
              write(*,*)'norm S_k: newton',l,k,sqrt(dot_product(S_k,S_k))
            else
              
              call amux(nveclen,S_InNewt,R_InNewt_predct_k,aJac,jaJac,iaJac)

              call NL_Optimization_gamma(S_Cauchy, R_Cauchy_predct_k, &
                                         S_InNewt, R_InNewt_predct_k, &
                                         delta_k, gamma_k) 
              S_k(:) = (1.0_wp - gamma_k) * S_Cauchy(:)  &
                      + (         gamma_k) * S_InNewt(:)
              write(*,*)'norm S_k: dogleg',l,k,sqrt(dot_product(S_k,S_k))
            endif
          endif
            
        endif

        S_L2_k = sqrt(dot_product(S_k,S_k))

        call amux(nveclen,S_k,J_S_k,aJac,jaJac,iaJac)
        R_NonLin_predct_k      = R_NonLin_k + J_S_k
        R_NonLin_predct_L2_k   = sqrt(dot_product(R_NonLin_predct_k, &
                                                  R_NonLin_predct_k)+small)


        R_NonLin_predct_L2_km1 = R_NonLin_predct_L2_k
        R_NonLin_L2_km1        = R_NonLin_L2_k

        uvec(:) = ustor(:) - S_k(:)

        !     F(x_k)  
        R_NonLin_k = Nonlinear_Residual(ep,dt,time,aI,iprob,L)
        !   ||F(x_k)||
        R_NonLin_L2_k  = sqrt(dot_product(R_NonLin_k(:),R_NonLin_k(:))+small)


        call NL_Optimization_Delta_k(k, R_NonLin_L2_k, R_NonLin_L2_km1,   &
                                        R_NonLin_predct_L2_k, S_L2_k,     &
                                        Delta_k)

        call NL_Optimization_Eta_k  (k, R_NonLin_L2_k, R_NonLin_L2_km1,   &
                                        R_NonLin_predct_L2_k,         Eta_k)
        print*,'L',L,'k',k,'tmp',R_NonLin_L2_k,'||Uvec||',sqrt(dot_product(uvec(:),uvec(:)))
        if(R_NonLin_L2_k <= 1.0e-10_wp) then
          ierr = 0
          icount = icount + k
          return
        endif

      enddo      
      icount = icount + k - 1

      end subroutine Inexact_Newton_Dogleg
      
!==============================================================================
!******************************************************************************
! Subroutine to perform QR decomposition and solving
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
! QR_MODULE.F90             *CONTAINS QR ROUTINES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From control_variables:
!   uvec -> Array containing variables, real(wp), dimension(u-vector length), modified
! From QR_Module:
!   qrdcmp -> performs QR decomposition
!   qrsolv -> performs QR solving 
! MODULE ROUTINES:
! Nonlinear_Residual -> Calls problemsub to build RHS and RNewton
! 
!******************************************************************************
! INPUTS:
! iprob   -> defines problem number,                                    integer
! L       -> stage number,                                              integer
! ep      -> Stiffness epsilon value,                                   real(wp)
! dt      -> timestep: created and then used later in the problem case, real(wp)
! nveclen -> total length of u-vector used for problem                  integer
! time    -> current solution time,                                     integer
! aI      -> Diagonal term from RK scheme,                              real(wp)
! INOUTS:
! icount  -> counter used to determine average iterations,              integer
! OUTPUTS:
! k       -> total number of newton iterations used for max iteration,  integer
!****************************************************************************** 
      subroutine QR_decomp(mat,nveclen,Rnewton)
      
      use control_variables, only: uvec
      use QR_Module,         only: qrdcmp,qrsolv
      
      real(wp), dimension(:,:), intent(inout) :: mat
      integer,                  intent(in   ) :: nveclen
      real(wp), dimension(:),   intent(inout) :: Rnewton
      
      real(wp), dimension(nveclen) :: wrk_c,wrk_d
      logical                      :: snglr
      
      call qrdcmp(mat,nveclen,nveclen,wrk_c,wrk_d,snglr)
      if (snglr) then
        write(*,*)'matrix is singular in Newton Iteration: Stopping'
        stop
      else
        call qrsolv(mat,nveclen,nveclen,wrk_c,wrk_d,Rnewton)
      endif
      uvec(:)=uvec(:)-Rnewton(:) 
      end subroutine QR_decomp       
      
!==============================================================================
!******************************************************************************
! Function to invert a dense matrix of rank<=4
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! 
!******************************************************************************
! INPUTS:
! mat -> input dense matrix, real(wp), dimension(:,:) 
! OUTPUTS:
! Mat_invert -> output inverted matrix, real(wp), dimension(:,:) (same as input)
!****************************************************************************** 
      function Mat_invert(mat)

      real(wp), dimension(:,:), intent(in   ) :: mat
      integer  :: mat_size
      real(wp), dimension(size(mat(:,1)),size(mat(1,:))) :: xinv
      real(wp), dimension(size(mat(:,1)),size(mat(1,:))) :: Mat_invert
      real(wp) :: x11,x12,x13,x14
      real(wp) :: x21,x22,x23,x24
      real(wp) :: x31,x32,x33,x34
      real(wp) :: x41,x42,x43,x44
      real(wp) :: det,detI   

      real(wp), parameter :: tol = 1.0e-20_wp

      mat_size=size(mat(:,1))
           
      if(mat_size==2)then

        det = (mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1))

        if(abs(det) <= tol)write(*,*)'determinant is nearly singular',det

        xinv(1,1) =  mat(2,2)/det
        xinv(1,2) = -mat(1,2)/det
        xinv(2,1) = -mat(2,1)/det
        xinv(2,2) =  mat(1,1)/det

      elseif(mat_size==3)then

        x11 = mat(1,1)
        x12 = mat(1,2)
        x13 = mat(1,3)
        x21 = mat(2,1)
        x22 = mat(2,2)
        x23 = mat(2,3)
        x31 = mat(3,1)
        x32 = mat(3,2)
        x33 = mat(3,3)

        det = - x13*x22*x31 + x12*x23*x31 +  x13*x21*x32& 
     &        - x11*x23*x32 - x12*x21*x33 +  x11*x22*x33

        if(abs(det) <= tol)write(*,*)'determinant is nearly singular',det

        detI = 1.0_wp/det

        xinv(1,1) =  (- x23*x32 + x22*x33) * detI
        xinv(1,2) =  (+ x13*x32 - x12*x33) * detI
        xinv(1,3) =  (- x13*x22 + x12*x23) * detI
        xinv(2,1) =  (+ x23*x31 - x21*x33) * detI
        xinv(2,2) =  (- x13*x31 + x11*x33) * detI
        xinv(2,3) =  (+ x13*x21 - x11*x23) * detI
        xinv(3,1) =  (- x22*x31 + x21*x32) * detI
        xinv(3,2) =  (+ x12*x31 - x11*x32) * detI
        xinv(3,3) =  (- x12*x21 + x11*x22) * detI

      elseif(mat_size==4)then

        x11 = mat(1,1)
        x12 = mat(1,2)
        x13 = mat(1,3)
        x14 = mat(1,4)
        x21 = mat(2,1)
        x22 = mat(2,2)
        x23 = mat(2,3)
        x24 = mat(2,4)
        x31 = mat(3,1)
        x32 = mat(3,2)
        x33 = mat(3,3)
        x34 = mat(3,4)
        x41 = mat(4,1)
        x42 = mat(4,2)
        x43 = mat(4,3)
        x44 = mat(4,4)

        det =& 
     &  (x14*x23*x32*x41 - x13*x24*x32*x41 - &
     &   x14*x22*x33*x41 + x12*x24*x33*x41 + &
     &   x13*x22*x34*x41 - x12*x23*x34*x41 - &
     &   x14*x23*x31*x42 + x13*x24*x31*x42 + &
     &   x14*x21*x33*x42 - x11*x24*x33*x42 - &
     &   x13*x21*x34*x42 + x11*x23*x34*x42 + &
     &   x14*x22*x31*x43 - x12*x24*x31*x43 - &
     &   x14*x21*x32*x43 + x11*x24*x32*x43 + &
     &   x12*x21*x34*x43 - x11*x22*x34*x43 - &
     &   x13*x22*x31*x44 + x12*x23*x31*x44 + &
     &   x13*x21*x32*x44 - x11*x23*x32*x44 - &
     &   x12*x21*x33*x44 + x11*x22*x33*x44)

        if(abs(det) <= tol)write(*,*)'determinant is nearly singular',det

        detI = 1.0_wp/det

        xinv(1,1) = (&
     & -(x24*x33*x42) + x23*x34*x42 + x24*x32*x43 - &
     &   x22*x34*x43  - x23*x32*x44 + x22*x33*x44  ) * detI
        xinv(1,2) = (&
     &   x14*x33*x42  - x13*x34*x42 - x14*x32*x43 + &
     &   x12*x34*x43  + x13*x32*x44 - x12*x33*x44  ) * detI
        xinv(1,3) = (&
     & -(x14*x23*x42) + x13*x24*x42 + x14*x22*x43 - &
     &   x12*x24*x43  - x13*x22*x44 + x12*x23*x44  ) * detI
        xinv(1,4) = (&
     &   x14*x23*x32  - x13*x24*x32 - x14*x22*x33 + &
     &   x12*x24*x33  + x13*x22*x34 - x12*x23*x34  ) * detI
        xinv(2,1) = (&
     &   x24*x33*x41  - x23*x34*x41 - x24*x31*x43 + &
     &   x21*x34*x43  + x23*x31*x44 - x21*x33*x44  ) * detI
        xinv(2,2) = (&
     & -(x14*x33*x41) + x13*x34*x41 + x14*x31*x43 - &
     &   x11*x34*x43  - x13*x31*x44 + x11*x33*x44  ) * detI
        xinv(2,3) = (&
     &   x14*x23*x41  - x13*x24*x41 - x14*x21*x43 + &
     &   x11*x24*x43  + x13*x21*x44 - x11*x23*x44  ) * detI
        xinv(2,4) = (&
     & -(x14*x23*x31) + x13*x24*x31 + x14*x21*x33 - &
     &   x11*x24*x33  - x13*x21*x34 + x11*x23*x34  ) * detI
        xinv(3,1) = (&
     & -(x24*x32*x41) + x22*x34*x41 + x24*x31*x42 - &
     &   x21*x34*x42  - x22*x31*x44 + x21*x32*x44  ) * detI
        xinv(3,2) = (&
     &   x14*x32*x41  - x12*x34*x41 - x14*x31*x42 + &
     &   x11*x34*x42  + x12*x31*x44 - x11*x32*x44  ) * detI
        xinv(3,3) = (&
     & -(x14*x22*x41) + x12*x24*x41 + x14*x21*x42 - &
     &   x11*x24*x42  - x12*x21*x44 + x11*x22*x44  ) * detI
        xinv(3,4) = (&
     &   x14*x22*x31  - x12*x24*x31 - x14*x21*x32 + &
     &   x11*x24*x32  + x12*x21*x34 - x11*x22*x34  ) * detI
        xinv(4,1) = (&
     &   x23*x32*x41  - x22*x33*x41 - x23*x31*x42 + &
     &   x21*x33*x42  + x22*x31*x43 - x21*x32*x43  ) * detI
        xinv(4,2) = (&
     & -(x13*x32*x41) + x12*x33*x41 + x13*x31*x42 - &
     &   x11*x33*x42  - x12*x31*x43 + x11*x32*x43  ) * detI
        xinv(4,3) = (&
     &   x13*x22*x41  - x12*x23*x41 - x13*x21*x42 + &
     &   x11*x23*x42  + x12*x21*x43 - x11*x22*x43  ) * detI
        xinv(4,4) = (&
     & -(x13*x22*x31) + x12*x23*x31 + x13*x21*x32 - &
     &   x11*x23*x32  - x12*x21*x33 + x11*x22*x33  ) * detI

      endif

      Mat_invert=xinv

      end function Mat_invert

!==============================================================================
!******************************************************************************
! Subroutine to convert a dense matrix to CSR and store it as the global CSR
! Jacobian
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! JACOBIAN_CSR_MOD.F90      *ALLOCATE AND STORE CSR JACOBIAN VARIABLES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From Jacobian_CSR_Mod:
!   iaJac  -> ia matrix for global storage of Jacobian, integer,  dimension(u-vector length + 1), set
!   jaJac  -> ja matrix for global storage of Jacobian, integer,  dimension(u-vector length **2), set
!    aJac  ->  a matrix for global storage of Jacobian, real(wp), dimension(u-vector length **2), set
!   Allocate_Jac_CSR_Storage -> Subroutine to create Jacobian and LU decomposition arrays for CSR problems
!******************************************************************************
! INPUTS:
! a        -> input dense matrix, real(wp), dimension(:,:) 
!****************************************************************************** 
      subroutine Convert_to_CSR(a)
 
      use Jacobian_CSR_Mod, only: Allocate_Jac_CSR_Storage,iaJac,jaJac,aJac      
      
      real(wp), dimension(:,:), intent(in   ) :: a
      integer                                 :: i,j,icnt,jcnt,nnz,dimen
      real(wp),    parameter              ::   toljac = 1.0e-13_wp     
      dimen=size(a(:,1))
      nnz=size(a)

      call Allocate_Jac_CSR_Storage(dimen,nnz)

!     U_t = F(U);  Jac = \frac{\partial F(U)}{\partial U};  xjac = I - akk dt Jac

      ! Initialize CSR
      iaJac(:) = 0
      iaJac(1) = 1
      jaJac(:) = 0 
      aJac(:) = 0.0_wp
      
      ! Store dense matrix into CSR format
      icnt = 0   
      do i = 1,dimen
        jcnt = 0   
        do j = 1,dimen
          if(abs(a(i,j)) >= tolJac) then
            icnt = icnt + 1 
            jcnt = jcnt + 1 
             
            jaJac(icnt) = j
             aJac(icnt) = a(i,j)
          endif
        enddo
        iaJac(i+1) = iaJac(i) + jcnt
      enddo

      end subroutine Convert_to_CSR

!==============================================================================      

      subroutine NL_Optimization_Delta_k(k, R_NonLin_L2_k, R_NonLin_L2_km1,   &
                                            R_NonLin_predct_L2_k, S_InNewt_L2,&
                                            Delta_k)

!============================================================================

!   Uses the current and previous nonlinear residuals and the current predicted residual
!
!   Inexact Newton Dogleg Methods
!   R. PAWLOWSKI, J. SIMONIS, H. WALKER, AND J. SHADID
!
!   Worcester Polytechnic Institute:     DigitalCommons@WPI
!   SIAM J. Numer. Anal, Vol. 46, No. 4, pp. 2112-2132,   (2008)
!
!
!   Integer:
!      k                                                    :  Nonlinear iteration

!   Real:
!      R_NonLin_L2_k        =  || F(x_k  ) ||                    :  Nonlinear residual at iteration k
!      R_NonLin_L2_km1      =  || F(x_km1) ||                    :  Nonlinear residual at iteration km1
!      R_NonLin_Predct_L2_k =  || F(x_km1) + F'(x_km1) S_km1 ||  :  Predicted nonlinear residual at k
!      S_InNewt_L2          =  || S_Inexact_Newton ||            :  Update provided by GMRES iteration
!      Delta_k                                                   :  
!

      integer,  intent(in)    :: k
      real(wp), intent(in)    :: R_NonLin_L2_k, R_NonLin_L2_km1
      real(wp), intent(in)    :: R_NonLin_Predct_L2_k, S_InNewt_L2
      real(wp), intent(inout) :: Delta_k
  
      real(wp)                :: act_red_k, pred_red_k
      real(wp), parameter     ::   rho_s   = 0.10_wp   ,   rho_e   = 0.75_wp
      real(wp), parameter     ::  beta_s   = 0.25_wp   ,  beta_e   = 4.00_wp
      real(wp), parameter     :: delta_min = 1.0e-03_wp, delta_max = 1.0e+08_wp
      real(wp), parameter     ::     eps   = 1.0e-10_wp

      continue

      if(k == 0) then

        if( S_InNewt_L2 < Delta_min) then
          Delta_k = 2.0_wp * Delta_min
        else
          Delta_k = S_InNewt_L2
        endif

      else

         act_red_k = R_NonLin_L2_km1 - R_NonLin_L2_k
        pred_red_k = R_NonLin_L2_km1 - R_NonLin_Predct_L2_k
  
        if(act_red_k / pred_red_k < rho_s) then
           if(S_InNewt_L2 <= Delta_k) then
              Delta_k = max(S_InNewt_L2,Delta_min)
           else
              Delta_k = max(beta_s*Delta_k,Delta_min)
           endif
        else
           if( (act_red_k / pred_red_k > rho_e) .and. (S_InNewt_L2 - Delta_k <= eps) ) then
             Delta_k = min(beta_e*Delta_k, Delta_max)
           else
             Delta_k = Delta_k
           endif
  
        endif
      endif

      end subroutine NL_Optimization_Delta_k

!============================================================================

      subroutine NL_Optimization_Eta_k(k, F_norm_k, F_norm_km1, &
                                       F_norm_pred_km1, Eta_k)

!============================================================================
!       Algorithm IN. Inexact newton method 
!           Let x0 be given.
!           For k = 0, 1, . . . (until convergence) do:
!           Choose η_k ∈ [0, 1) and s^IN_k such that
!           ||F(x_k) + F'(x_k) s^IN_k|| ≤ \eta_k || F(x_k) || .
!           Set x_{k+1} = x_k + s^IN_k.

!   Adjust the forcing term \eta_k region  
!
!   Uses the current and previous nonlinear residuals and the current predicted residual
!
!   Inexact Newton Dogleg Methods
!   R. PAWLOWSKI, J. SIMONIS, H. WALKER, AND J. SHADID
!
!   Worcester Polytechnic Institute:     DigitalCommons@WPI
!   SIAM J. Numer. Anal, Vol. 46, No. 4, pp. 2112-2132,   (2008)
!
!
!   Integer:
!      k                                                    :  Nonlinear iteration

!   Real:
!      F_norm_k        =  || F(x_k  ) ||                    :  Nonlinear residual at iteration k
!      F_norm_km1      =  || F(x_km1) ||                    :  Nonlinear residual at iteration km1
!      F_norm_pred_k   =  || F(x_km1) + F'(x_km1) S_km1 ||  :  Predicted nonlinear residual at k
!      Eta_k                                                :  Trust region radius
!
!============================================================================

      implicit none

      integer,  intent(in)    :: k
      real(wp), intent(in)    :: F_norm_k,F_norm_km1,F_norm_pred_km1
      real(wp), intent(inout) :: Eta_k
  
      real(wp)                :: Eta_k_tmp, expo
      real(wp), parameter     :: Eta_max = 0.9_wp, Eta_0 = 0.01_wp

      real(wp), save          :: Eta_km1
  
      continue

      Eta_k_tmp = min(Eta_max, abs( F_norm_k - F_norm_pred_km1 ) / F_norm_km1)
      expo = (1.0_wp+sqrt(5.0_wp))/2.0_wp
      if( Eta_km1**expo  > 0.1_wp ) then
        Eta_k = max(Eta_k_tmp, Eta_km1**expo)
      else
        Eta_k = Eta_k_tmp
      endif

      Eta_km1 = Eta_k

      end subroutine NL_Optimization_Eta_k

!============================================================================

      subroutine NL_Optimization_gamma(S_Cauchy, R_Cauchy_predct, &
                                       S_InNewt, R_InNewt_predct, &
                                       delta_k, gamma_k) 

      implicit none

      real(wp), dimension(:), intent(in   )    :: S_Cauchy, R_Cauchy_predct
      real(wp), dimension(:), intent(in   )    :: S_InNewt, R_InNewt_predct
      real(wp),               intent(in   )    :: delta_k
      real(wp),               intent(  out)    :: gamma_k

      real(wp), dimension(size(S_Cauchy))      :: wrk

      real(wp)                                 :: gam_min, gam_P
      real(wp)                                 :: T1, T2, T3, T4

  
      continue

      wrk(:) = R_Cauchy_predct(:) - R_InNewt_predct(:)

      gam_min = dot_product(R_Cauchy_predct,wrk) /  &
              & dot_product(            wrk,wrk) 

      wrk(:) = S_Cauchy(:) - S_InNewt(:)

      T1 = dot_product(S_Cauchy, wrk     )
      T2 = dot_product(S_Cauchy, S_Cauchy)
      T3 = dot_product(    wrk , wrk     )
      T4 = delta_k**2

      if (T4 >= T2) then
        gam_P = (T1 + sqrt( T1 * T1 + (T4 - T2) * T3 )) / T3
      else
        write(*,*)'somethings wrong in NL_Optimization_gamma' ; stop ;
      endif

      gamma_k = min(gam_min, gam_P)

      end subroutine NL_Optimization_gamma

      end module Newton
