!******************************************************************************
! Module to output data for test_cases.f90 to data files in the form:
! /"probname"/"casename"/"probname"_"casename"_"variable number".dat
! /"probname"/"casename"/"probname"_"casename"_"variable number"P.dat
! /"probname"/"casename"/"probname"_"casename"_conv.dat
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90            *DEFINES PRECISION FOR ALL VARIABLES
! POLY_FIT_MOD.F90              *CONTAINS SUBROUTINES TO OUTPUT DATA
! CONTROL_VARIABLES.F90         *ONTAINS VARIABLES USED IN THE PROGRAM
! RUNGE_KUTTA.F90               *CONTAINS RK CONSTANTS
!******************************************************************************

      module  output_module

      use precision_vars,    only : wp

      implicit none; save
      
      public :: create_file_paths, output_names, init_output_files
      public :: output_conv_stiff,output_terminal_iteration
      public :: output_terminal_final,output_conv_error,write_time_depen_sol
      public :: output_iDT_sol,write_exact_sol
      
      private
      character(len=80)            :: fileloc   !file  location
      character(len=4), parameter  :: ext='.dat'  !file extension

      contains
      
!==============================================================================
!******************************************************************************
! Set file location variable
!******************************************************************************
! REQUIRED FILES:
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
! RUNGE_KUTTA.F90           *CONTAINS RK CONSTANTS
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From control_variables:
!   temporal_splitting -> string used in choose_RHS_type and choose_Jac_type, character(len=80), not modified
!   probname           -> string used to set output file names,               character(len=9),  not modified
! From runge_kutta:
!   casename           -> string holding name of RK case,                     character(len=25), not modified
!
! MODULE VARIABLES/ROUTINES:
! fileloc -> string used to hold the filename of the output files,            character(len=80), set
!******************************************************************************
      subroutine set_file_loc()
      
      use control_variables, only : temporal_splitting,probname
      use runge_kutta,       only : casename
      
      fileloc=trim(temporal_splitting)//'/'// &
     &        trim(probname)//'/'//trim(casename)//'/' ! get file location  
      
      end subroutine set_file_loc           
      
!==============================================================================
!******************************************************************************
! Creates file paths for output files if they do not exsist
!******************************************************************************
! MODULE VARIABLES/ROUTINES:
! fileloc -> string used to hold the filename of the output files, character(len=80), not modified
! set_file_loc -> Subroutine for setting fileloc
!******************************************************************************
      subroutine create_file_paths()
      
      logical           :: dirExists   !check if directory exists
      
      call set_file_loc
      inquire( file=trim(fileloc)//'/.', exist=dirExists ) 
      if (.not.dirExists)call system('mkdir -p '//trim(fileloc)) 
       
      end subroutine create_file_paths
      
!==============================================================================  
!******************************************************************************
! Outputs RK case, IMEX/IMPLICIT/EXPLICIT, and problem names to the terminal
!******************************************************************************
! REQUIRED FILES:
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
! RUNGE_KUTTA.F90           *CONTAINS RK CONSTANTS
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From control_variables:
!   temporal_splitting -> string used in choose_RHS_type and choose_Jac_type, character(len=80), not modified
!   probname           -> string used to set output file names,               character(len=9),  not modified
! From runge_kutta:
!   casename           -> string holding name of RK case,                     character(len=25), not modified
!******************************************************************************
      subroutine output_names()

      use control_variables, only : probname,temporal_splitting
      use runge_kutta,       only : casename
       
      write(*,*)'Problem name: ',trim(probname)
      write(*,*)'Casename: ', trim(casename)
      write(*,*)'Splitting: ', trim(temporal_splitting)
       
      end subroutine output_names
       
!==============================================================================
!******************************************************************************
! Create zones in output files for Tecplot/Tec360
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
! RUNGE_KUTTA.F90           *CONTAINS RK CONSTANTS
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From control_variables:
!   probname -> string used to set output file names,              character(len=9),  not modified
! From runge_kutta:
!   casename -> string holding name of RK case,                    character(len=25), not modified
!
! MODULE VARIABLES/ROUTINES:
! fileloc -> string used to hold the filename of the output files, character(len=80), not modified
! ext     -> string used for file extension,                       character(len=4),  not modified
!******************************************************************************
! INPUTS:
! neq  -> number of equations in problem, integer
! ep   -> Stiffness epsilon value,        real(wp)
!******************************************************************************

      subroutine init_output_files(neq,ep)

      use control_variables, only : probname
      use runge_kutta,       only : casename
     
      integer,  intent(in) :: neq
      real(wp), intent(in) :: ep
      character(len=80)    :: filename   !file name and location
      character(len=1)     :: istr        !loop index placeholder
      integer              :: i
          
      do i = 1,neq
        write(istr,"(I1.1)")i
        filename=trim(probname)//'_'//trim(casename)//'_'//istr//ext
        open(49+i,file=trim(fileloc)//filename)
        write(49+i,*)'zone T = "ep = ',ep,'",'
        filename=trim(probname)//'_'//trim(casename)//'_'//istr//'P'//ext
        open(59+i,file=trim(fileloc)//filename)
        write(59+i,*)'zone T = "ep = ',ep,'",'
      enddo

      end subroutine init_output_files
      
!==============================================================================      
!******************************************************************************
! Output data to convergence vs. stiffness files
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
! RUNGE_KUTTA.F90           *CONTAINS RK CONSTANTS
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From control_variables:
!   temporal_splitting -> string used in choose_RHS_type and choose_Jac_type, character(len=80),                               not modified
!   probname           -> string used to set output file names,               character(len=9),                                not modified
!   jactual            -> actual number of epsilon values,                    integer,                                         not modified
!   b1save             -> convergence rate storage,                           real(wp), dimension(num. ep's, u-vector length), not modified
!   b1Psave            -> convergence rate storage, predicted,                real(wp), dimension(num. ep's, u-vector length), not modified
!   b1L2save           -> convergence rate storage, L2 norm ,                 real(wp), dimension(num. ep's, num. eq's),       not modified
!   var_names          -> string array used to label output graphs,           character(len=12), dimension(num. eq's),         not modified
! From runge_kutta:    
!   casename           -> string holding name of RK case,                     character(len=25),                               not modified
!
! MODULE VARIABLES/ROUTINES:
! fileloc -> string used to hold the filename of the output files, character(len=80), not modified
! ext     -> string used for file extension,                       character(len=9),  not modified
!******************************************************************************
! INPUTS:
! nveclen -> u-vector length,                 integer
! neq     -> number of equations,             integer
! epsave  -> array containing epsilon values, real(wp), dimension(num. ep's)
!******************************************************************************
      subroutine output_conv_stiff(nveclen,neq,epsave)
      
      use control_variables, only : temporal_splitting,probname,jactual, &
     &                              var_names,b1L2save!,b1save,b1Psave  
      use runge_kutta,       only : casename    
      
      integer,                intent(in) :: nveclen,neq
      real(wp), dimension(:), intent(in) :: epsave
      
      integer           :: i,j
      character(len=80) :: filename   !file name and location
      
      filename=trim(probname)//'_'//trim(casename)//'_conv'//ext
      
      open(35,file=trim(fileloc)//filename)
!      do i=1,nveclen
!        write(35,*)'zone T = "Var ',i,': ',trim(temporal_splitting),'",'
!        do j=1,jactual
!          write(35,50)epsave(j),b1save(j,i)
!        enddo
!          write(35,*)'zone T = "Var ',i,': Predicted",'
!        do j=1,jactual
!          write(35,50)epsave(j),b1Psave(j,i)
!        enddo
!      enddo
           
      do i=1,neq
        write(35,*)'zone T = "',trim(var_names(i)),' Variable - ',trim(temporal_splitting),'",'
        do j=1,jactual
          write(35,50)epsave(j),b1L2save(j,i)
        enddo
      enddo

      50 format( 10(e12.5,1x))
         
      end subroutine output_conv_stiff

!==============================================================================
!******************************************************************************
! Outputs data to terminal every epsilon iteration
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! POLY_FIT_MOD.F90          *CONTAINS SUBROUTINES TO OUTPUT DATA
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From poly_fit_Mod:
!   fit          -> Subroutine to get convergence slopes
! From control_variables:
!   isamp        -> number of dt's,                            integer,                                              not modified
!   dt_error_tol -> L2 error tolerance for convergence plots,  real(wp),                                             not modified
!   error        -> solution error,                            real(wp), dimension(num. dt's, u-vector length),      not modified
!   errorP       -> solution error, predicted,                 real(wp), dimension(num. dt's, u-vector length),      not modified
!   b            -> convergence rate,                          real(wp), dimension(2 * u-vector length + num. eq's),     set          
!   errorL2      -> solution error, L2 norm,                   real(wp), dimension(num. dt's, num. eq's),            not modified
!******************************************************************************
! INPUTS:
! cost    -> cost of a particular dt value,             real(wp)
! mwt     -> fit subroutine input,                      integer
! ep      -> Stiffness epsilon value,                   real(wp)
! nveclen -> total length of u-vector used for problem, integer
! neq     -> number of equations in problem,            integer
!******************************************************************************
      subroutine output_terminal_iteration(cost,mwt,ep,nveclen,neq)
      
      use poly_fit_Mod,      only : fit
      use control_variables, only : isamp,dt_error_tol,b,error,errorP,errorL2
            
      real(wp), dimension(:),   intent(in) :: cost

      integer,                  intent(in) :: mwt
      real(wp),                 intent(in) :: ep
      integer,                  intent(in) :: nveclen,neq
      
      real(wp), dimension(nveclen*2+neq) :: a,siga1,sigb1,chi2
      real(wp)                       :: q
      integer                        :: i,j
      integer                        :: jsamp
      real(wp), dimension(isamp)     :: sig
      
      sig(:)=0.0_wp


      !**GATHER OUTPUT VALUES

      do i = 1,nveclen
        do j=2,size(error(:,i))
          jsamp=j
          if (error(j,i)<=log10(dt_error_tol)) exit
        enddo
 
        call fit(cost,error(:,i),jsamp,sig,mwt,a(i),b(i),                  &
     &           siga1(i),sigb1(i),chi2(i),q)
        call fit(cost,errorP(:,i),jsamp,sig,mwt,a(i+nveclen),b(i+nveclen), &
     &           siga1(i+nveclen),sigb1(i+nveclen),chi2(i+nveclen),q)   
     enddo
     
     do i = 1,neq
        do j=2,size(errorL2(:,i))
          jsamp=j
          if (errorL2(j,i)<=log10(dt_error_tol)) exit
        enddo

        call fit(cost,errorL2(:,i),jsamp,sig,mwt,a(nveclen*2+i),b(nveclen*2+i),                  &
     &           siga1(nveclen*2+i),sigb1(nveclen*2+i),chi2(nveclen*2+i),q)
      enddo

      !**OUTPUT TO TERMINAL**
      write(*,60,advance="no")ep
! HACK - Every variables's convergence suppressed
!        do i = 1,nveclen*2+1
!          write(*,60,advance="no")a(i),b(i)
!        enddo
!      write(*,60)a(nveclen*2+2),b(nveclen*2+2)
! HACK
      ! L2 norm convergences for each variable 
      do i = 1,neq-1
        write(*,60,advance="no")a(nveclen*2+i),b(nveclen*2+i)
      enddo
      write(*,60)a(nveclen*2+neq),b(nveclen*2+neq)
   60 format( e12.5,1x,12(f8.3,1x))
   
      end subroutine output_terminal_iteration

!==============================================================================
!******************************************************************************
! Outputs data to terminal at end of RK case/problem
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp -> working precision
! From runge_kutta:
!   ns -> number of RK stage for this particular scheme, integer
!******************************************************************************
! INPUTS:
! icount  -> counter used to determine average iterations, integer
! jcount  -> total number of RK stages,                    integer
! stageE  -> ?                                             real(wp), dimension(max stages)
! stageI  -> ?                                             real(wp), dimension(max stages)
! maxiter -> maximum number of newton iterations           real(wp), dimension(max stages)
!******************************************************************************
      subroutine output_terminal_final(icount,jcount,stageE,stageI,maxiter)
      
      use runge_kutta, only : ns
      
      integer,                intent(in   ) :: icount,jcount
      real(wp), dimension(:), intent(inout) :: stageE,stageI,maxiter
      
      real(wp) :: tmp
      integer  :: i
      
      tmp = 1.0_wp*icount/jcount
      write(*,*)'average iterations per step',tmp
             
      stageE(2:) = (ns-1)*stageE(2:)/jcount
      stageI(2:) = (ns-1)*stageI(2:)/jcount

      write(*,*)'error of initial guess is '
      do i = 2,ns
        write(*,*) i,stageE(i),maxiter(i),stageI(i)
      enddo
      
      end subroutine output_terminal_final

!==============================================================================
!******************************************************************************
! Outputs data to files initilized earilier in the program
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From control_variables:
!   uvec    -> Array containing variables,                   real(wp), dimension(u-vector length), not modified
!   uexact  -> Array containing exact solution to variables, real(wp), dimension(u-vector length), not modified
!   errvecT -> error storage for each time step,             real(wp), dimension(u-vector length), not modified
!******************************************************************************
! INPUTS:
! cost    -> cost of a particular dt value,             real(wp)
! nveclen -> total length of u-vector used for problem, integer
! neq     -> number of equations in problem,            integer
!******************************************************************************
      subroutine output_conv_error(cost,nveclen,neq)
      
      use control_variables, only : uvec,uexact,errvecT
      
      real(wp),         intent(in) :: cost
      integer,          intent(in) :: nveclen,neq
      
      integer                      :: i
      real(wp), dimension(nveclen) :: tmp
     
      tmp=abs(uvec(:)-uexact(:))
      do i = 1,nveclen
       if (tmp(i)==0.0_wp)tmp(i)=1.0e-15_wp  
      enddo
      
      do i = 1,neq
        write(49+i,*)cost,log10(sqrt(dot_product(tmp(i:nveclen:neq), &
     &                                   tmp(i:nveclen:neq))/(nveclen/neq)))
        write(59+i,*)cost,log10(sqrt(dot_product(errvecT(i:nveclen:neq), &
     &                                errvecT(i:nveclen:neq))/(nveclen/neq)))
      enddo

      
      end subroutine output_conv_error
!==============================================================================
!******************************************************************************
! Outputs time dependant solution
!******************************************************************************
! REQUIRED FILES:
! PRECISION_VARS.F90        *DEFINES PRECISION FOR ALL VARIABLES
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
! RUNGE_KUTTA.F90           *CONTAINS RK CONSTANTS
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From precision_variables:
!   wp  -> working precision
! From control_variables:
!   uvec    -> Array containing variables,                   real(wp), dimension(u-vector length), not modified
!   probname           -> string used to set output file names,               character(len=9),  not modified
! From runge_kutta:
!   casename           -> string holding name of RK case,                     character(len=25), not modified
!
! MODULE VARIABLES/ROUTINES:
! fileloc -> string used to hold the filename of the output files, character(len=80), not modified
! ext     -> string used for file extension,                       character(len=9),  not modified
!******************************************************************************
! INPUTS:
! time          -> current solution time,                     real(wp)
!  ep           -> Stiffness epsilon value,                   real(wp)
! nveclen       -> total length of u-vector used for problem, integer
! neq           -> number of equations in problem,            integer
!******************************************************************************
      subroutine write_time_depen_sol(time,ep,nveclen,neq)
      
      use control_variables, only: uvec,probname
      use runge_kutta,       only: casename
      
      real(wp), intent(in) :: time,ep
      integer,  intent(in) :: nveclen,neq
      integer              :: i
      character(len=80)    :: filename   !file name
      character(len=1)     :: istr
      
      do i = 1,neq  
          
        if (time<=1e-16_wp) then ! check if new epsilon
          write(istr,"(I1.1)")i
          filename=trim(probname)//'_'//trim(casename)//'_time_dependant_solution_'//istr//ext
          open(849+i,file=trim(fileloc)//filename)   
          write(849+i,*)'zone T = "ep = ',ep,'",'
        endif
         
        write(849+i,*)time,uvec(i:nveclen:neq)
      enddo
      
      end subroutine write_time_depen_sol  
!==============================================================================
!******************************************************************************
! Outputs solution for each iDT value for use in finding exact solution
!******************************************************************************
! REQUIRED FILES:
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From control_variables:
!   uvec    -> Array containing variables,                   real(wp), dimension(u-vector length), not modified
! 
! MODULE VARIABLES/ROUTINES:
! fileloc -> string used to hold the filename of the output files, character(len=80), not modified
! ext     -> string used for file extension,                       character(len=9),  not modified
!******************************************************************************
! INPUTS:
! iDT     -> Timestep counter from timestep loop to define dt,          integer
!******************************************************************************
      subroutine output_iDT_sol(iDT)
      
      use control_variables, only: uvec
      
      integer,  intent(in) :: iDT
      character(len=80)    :: filename   !file name
      character(len=2)     :: istr
      
      write(istr,"(I2.1)")iDT
      filename='iDT_solution_'//trim(istr)//ext
      open(899+iDT,file=trim(fileloc)//filename)          
      write(899+iDT,*)uvec
      close(899+iDT)
      
      end subroutine output_iDT_sol 
!==============================================================================
!******************************************************************************
! Outputs exact solution to a file in the appropriate directory
!******************************************************************************
! REQUIRED FILES:
! CONTROL_VARIABLES.F90     *CONTAINS VARIABLES AND ALLOCATION ROUTINES
!******************************************************************************
! GLOBAL VARIABLES/ROUTINES:
! From control_variables:
!   uvec    -> Array containing variables,                   real(wp), dimension(u-vector length), not modified
!   probname           -> string used to set output file names,               character(len=9),  not modified
!
! MODULE VARIABLES/ROUTINES:
! fileloc -> string used to hold the filename of the output files, character(len=80), not modified
! ext     -> string used for file extension,                       character(len=9),  not modified
!******************************************************************************
      subroutine write_exact_sol()
      
      use control_variables, only: uvec,probname
      
      character(len=80)    :: filename   !file name
      
      filename='exact_'//trim(probname)//ext
      open(95,file=trim(fileloc)//filename)          
      write(95,*)uvec
      
      end subroutine write_exact_sol
!==============================================================================
      end module output_module
