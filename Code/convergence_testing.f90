      program test_convergence

      implicit none

      integer, parameter      :: wp = 8
      integer, parameter      :: vLen = 3, samps = 81

      real(wp), dimension(vLen,samps) :: uexact, tester

      integer                         :: i,j
      real(wp)                        :: norm_2


      open(unit=40,file='fort.exact.data')
      open(unit=50,file='fort.test.data')

      norm_2 = 0.0_wp
      do j = 1,samps
         read(40,*) uexact(:,j)
         read(50,*) tester(:,j)
        norm_2 = norm_2 + dot_product(uexact(:,j)-tester(:,j),uexact(:,j)-tester(:,j))
      enddo
      write(*,*)'L2 norm of error =  ',sqrt(norm_2)

      end program
