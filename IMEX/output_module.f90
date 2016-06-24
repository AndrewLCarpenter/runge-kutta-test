      module  output_module

      use precision_vars     , only : wp
      use poly_fit_Mod

      implicit none
      
      public :: create_file_paths, output_names, init_output_files
      public :: output_conv_stiff,output_terminal_iteration
      public :: output_terminal_final,output_conv_error
      private

      contains
      
!==============================================================================

      subroutine create_file_paths(probname,casename)
      
      character(len=*), intent(in) :: probname,casename
      character(len=80) :: fileloc   !file name and location
      logical           :: dirExists   !check if directory exists
      
      fileloc=trim(probname)//'/'//trim(casename)//'/' ! create directory
      inquire( file=trim(fileloc)//'/.', exist=dirExists ) 
      if (.not.dirExists)call system('mkdir -p '//trim(fileloc)) 
       
      end subroutine create_file_paths
      
!==============================================================================  

      subroutine output_names(probname,casename)
       
      character(len=*), intent(in) :: probname,casename
       
      write(*,*)'Problem name: ',trim(probname)
      write(*,*)'Casename: ', trim(casename)
       
      end subroutine output_names
       
!==============================================================================

      subroutine init_output_files(probname,casename,nveclen,ep)
      
      character(len=*), intent(in) :: probname,casename
      integer,          intent(in) :: nveclen
      character(len=4)             :: ext='.dat'  !file extension
      character(len=1)             :: istr        !loop index placeholder
      integer                      :: i
      real(wp)                     :: ep
      character(len=80) :: fileloc,filename   !file name and location
      
      fileloc=trim(probname)//'/'//trim(casename)//'/' ! get file location      
      do i = 1,nveclen
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

      end subroutine init_output_files
      
!==============================================================================      

      subroutine output_conv_stiff(probname,casename,nveclen,jactual,epsave, &
     &                             b1save,b1Psave)
      
      integer, parameter                   :: jmax=81
      
      character(len=*),         intent(in) :: probname,casename
      integer,                  intent(in) :: nveclen,jactual
      real(wp), dimension(:,:), intent(in) :: b1save,b1Psave
      real(wp), dimension(:),   intent(in) :: epsave
      
      character(len=4)             :: ext='.dat'  !file extension
      character(len=1)             :: istr        !loop index placeholder
      integer                      :: i,j
      real(wp)                     :: ep
      character(len=80) :: fileloc,filename   !file name and location
     
      fileloc= trim(probname)//'/'//trim(casename)//'/' ! get file location   
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

      50 format( 10(e12.5,1x))
         
      end subroutine output_conv_stiff

!==============================================================================
      subroutine output_terminal_iteration(cost,error,errorP,jsamp,sig,mwt, &
     &           ep,nveclen,b)
            
      real(wp), dimension(:),   intent(in) :: cost      
      real(wp), dimension(:,:), intent(in) :: error,errorP    
      integer,                  intent(in) :: jsamp
      real(wp), dimension(:),   intent(in) :: sig
      integer,                  intent(in) :: mwt
      real(wp),                 intent(in) :: ep
      integer,                  intent(in) :: nveclen
      real(wp), dimension(nveclen*2), intent(out) :: b
      
      real(wp), dimension(nveclen*2) :: a,siga1,sigb1,chi2
      real(wp)                       :: q
      integer                        :: i

            !**GATHER OUTPUT VALUES
            do i = 1,nveclen
              call fit(cost,error(:,i),jsamp,sig,mwt,a(i),&
     &         b(i),siga1(i),sigb1(i),chi2(i),q)
              call fit(cost,errorP(:,i),jsamp,sig,mwt,a(i+nveclen),&
     &         b(i+nveclen),siga1(i+nveclen),sigb1(i+nveclen),&
     &         chi2(i+nveclen),q)
            enddo

            !**OUTPUT TO TERMINAL**
            write(*,60,advance="no")ep
            do i = 1,nveclen*2-1
              write(*,60,advance="no")a(i),b(i)
            enddo
            write(*,60)a(nveclen*2),b(nveclen*2)

   
   60 format( e12.5,1x,12(f8.3,1x))
      end subroutine output_terminal_iteration

!==============================================================================

      subroutine output_terminal_final(icount,jcount,nrk,stageE,stageI,maxiter)
      
      integer,                intent(in   ) :: icount,jcount,nrk
      real(wp), dimension(:), intent(inout) :: stageE,stageI,maxiter
      
      real(wp) :: tmp
      integer  :: i
      
      tmp = 1.0_wp*icount/jcount
      write(*,*)'average iterations per step',tmp
             
      stageE(2:) = (nrk-1)*stageE(2:)/jcount
      stageI(2:) = (nrk-1)*stageI(2:)/jcount

      write(*,*)'error of initial guess is '
      do i = 2,nrk
        write(*,*) i,stageE(i),maxiter(i),stageI(i)
      enddo
      
      end subroutine output_terminal_final

!==============================================================================

      subroutine output_conv_error(cost,uvec,uexact,errvecT)
      real(wp),               intent(in) :: cost
      real(wp), dimension(:), intent(in) :: uvec,uexact,errvecT
      
      integer  :: i
      real(wp) :: tmp
  
      do i = 1,size(uvec)
        tmp=abs(uvec(i)-uexact(i))
        if (tmp==0.0_wp)tmp=1.0e-15_wp
        write(49+i,*)cost,log10(tmp)
        write(59+i,*)cost,log10(errvecT(i))
      enddo
      
      end subroutine output_conv_error

!==============================================================================
      end module output_module
