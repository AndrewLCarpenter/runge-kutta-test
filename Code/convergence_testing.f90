      program test_convergence

      implicit none

      integer, parameter      :: wp = 8
      integer, parameter      :: vLen = 96, samps = 1

      real(wp), dimension(vLen,samps) :: uexact, tester

      integer                         :: i,j
      real(wp)                        :: norm_2
      real(wp), dimension(70)         :: norm_vec
      character(3) :: istr
      character(8) :: filename


      do i=121,190 !for 71 time steps
        write(istr,"(I3.1)") (i+1)

        open(unit=40,file='fort.'//istr)
        write(istr,"(I3.1)")i
        open(unit=50,file='fort.'//istr)

        norm_2 = 0.0_wp
        do j = 1,samps
         read(40,*) uexact(:,j)
         read(50,*) tester(:,j)
      
         norm_2 = norm_2 + dot_product(uexact(:,j)-tester(:,j),uexact(:,j)-tester(:,j))
        enddo
        write(*,*)'L2 norm of error =  ',sqrt(norm_2)
        norm_vec(i-120)=sqrt(norm_2)
       
        close(unit=40)
        close(unit=50)
      
      enddo
      write(*,*)'Minimum Norm = ',minval(norm_vec)
      write(*,*)'Location = ',minloc(norm_vec)
      
      
      end program
