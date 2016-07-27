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
      public :: output_terminal_final,output_conv_error
      
      private
      character(len=80)            :: fileloc   !file  location
      character(len=4), parameter  :: ext='.dat'  !file extension

      contains
      
!==============================================================================
!  SETS FILE LOCATION VARIABLE
      subroutine set_file_loc()
      
      use control_variables, only : temporal_splitting,probname
      use runge_kutta,       only : casename
      
      fileloc=trim(temporal_splitting)//'/'// &
     &        trim(probname)//'/'//trim(casename)//'/' ! get file location  
      
      end subroutine set_file_loc           
      
!==============================================================================
!  CHECKS AND CREATES FILE PATHS FOR OUTPUTS
      subroutine create_file_paths()
      
      logical           :: dirExists   !check if directory exists
      
      call set_file_loc
      inquire( file=trim(fileloc)//'/.', exist=dirExists ) 
      if (.not.dirExists)call system('mkdir -p '//trim(fileloc)) 
       
      end subroutine create_file_paths
      
!==============================================================================  
!  WRITES PROBLEM NAME AND CASE NAME TO THE TOP OF THE TERMIINAL OUTPUT
      subroutine output_names()

      use control_variables, only : probname,temporal_splitting
      use runge_kutta,       only : casename
       
      write(*,*)'Problem name: ',trim(probname)
      write(*,*)'Casename: ', trim(casename)
      write(*,*)'Splitting: ', trim(temporal_splitting)
       
      end subroutine output_names
       
!==============================================================================
!  INITIALIZE OUTPUT FILES FOR ERROR OUTPUTS
      subroutine init_output_files(nveclen,ep)

      use control_variables, only : probname
      use runge_kutta,       only : casename
     
      integer,  intent(in) :: nveclen
      real(wp), intent(in) :: ep
      character(len=80)    :: filename   !file name and location
      character(len=1)     :: istr        !loop index placeholder
      integer              :: i
          
      if(nveclen <= 9) then
        do i = 1,nveclen
          write(istr,"(I1.1)")i
          filename=trim(probname)//'_'//trim(casename)//'_'//istr//ext
          open(49+i,file=trim(fileloc)//filename)
          write(49+i,*)'zone T = "ep = ',ep,'",'
          filename=trim(probname)//'_'//trim(casename)//'_'//istr//'P'//ext
          open(59+i,file=trim(fileloc)//filename)
          write(59+i,*)'zone T = "ep = ',ep,'",'
        enddo
      endif

      end subroutine init_output_files
      
!==============================================================================      
!  OUTPUTS DATA TO CONVERGENCE VS. STIFFNESS FILE
      subroutine output_conv_stiff(nveclen,epsave)
      
      use control_variables, only : temporal_splitting,probname,jactual, &
     &                              b1save,b1Psave,b1L2save  
      use runge_kutta,       only : casename    
      
      integer,                intent(in) :: nveclen
      real(wp), dimension(:), intent(in) :: epsave
      
      integer           :: i,j
      character(len=80) :: filename   !file name and location
      
      filename=trim(probname)//'_'//trim(casename)//'_conv'//ext
      
      open(35,file=trim(fileloc)//filename)
      do i=1,nveclen
        write(35,*)'zone T = "Var ',i,': ',trim(temporal_splitting),'",'
        do j=1,jactual
          write(35,50)epsave(j),b1save(j,i)
        enddo
          write(35,*)'zone T = "Var ',i,': Predicted",'
        do j=1,jactual
          write(35,50)epsave(j),b1Psave(j,i)
        enddo
      enddo
      
      do i=1,2
        write(35,*)'zone T = "Var',i,': ',i,'",'
        do j=1,jactual
          write(35,50)epsave(j),b1L2save(j,i)
        enddo
      enddo

      50 format( 10(e12.5,1x))
         
      end subroutine output_conv_stiff

!==============================================================================
!  OUTPUTS ITERATION INFORMATION TO TERMINAL      
      subroutine output_terminal_iteration(cost,mwt,ep,nveclen)
      
      use poly_fit_mod,      only : fit
      use control_variables, only : isamp,dt_error_tol,b,error,errorP,errorL2
            
      real(wp), dimension(:),   intent(in) :: cost

      integer,                  intent(in) :: mwt
      real(wp),                 intent(in) :: ep
      integer,                  intent(in) :: nveclen
      
      real(wp), dimension(nveclen*2+2) :: a,siga1,sigb1,chi2
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
     
     do i = 1,2
        do j=2,size(errorL2(:,i))
          jsamp=j
          if (errorL2(j,i)<=log10(dt_error_tol)) exit
        enddo
        call fit(cost,errorL2(:,i),jsamp,sig,mwt,a(nveclen*2+i),b(nveclen*2+i),                  &
     &           siga1(nveclen*2+i),sigb1(nveclen*2+i),chi2(nveclen*2+i),q)
      enddo

      !**OUTPUT TO TERMINAL**
      write(*,60,advance="no")ep
! HACK
!        do i = 1,nveclen*2+1
!          write(*,60,advance="no")a(i),b(i)
!        enddo
!      write(*,60)a(nveclen*2+2),b(nveclen*2+2)
      write(*,60,advance="no")a(nveclen*2+1),b(nveclen*2+1)
      write(*,60)a(nveclen*2+2),b(nveclen*2+2)
   60 format( e12.5,1x,12(f8.3,1x))
   
      end subroutine output_terminal_iteration

!==============================================================================
!  OUTPUTS FINAL DATA BLOCK TO TERMINAL. ITERATION INFORMATION
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
!  OUTPUTS CONVERGENCE AND ERROR DATA TO FILES INITIALIZED EARILER
      subroutine output_conv_error(cost)
      
      use control_variables, only : uvec,uexact,errvecT
      
      real(wp),               intent(in) :: cost
      
      integer  :: i
      real(wp) :: tmp
  
      if (size(uvec)<=9) then
        do i = 1,size(uvec)
          tmp=abs(uvec(i)-uexact(i))
          if (tmp==0.0_wp)tmp=1.0e-15_wp
          write(49+i,*)cost,log10(tmp)
          write(59+i,*)cost,log10(errvecT(i))
        enddo
      endif
      
      end subroutine output_conv_error

!==============================================================================
      end module output_module
