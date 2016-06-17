module precision_vars
  implicit none
  
  ! Generates the kind number for double precision
  integer, parameter :: dp = kind(1.d0)
  ! Generate the kind number for single precision
  integer, parameter :: sp = kind(1.e0)
  
  ! Set the working precision
  integer, parameter :: wp = dp

  integer, parameter :: prec = dp

  ! Common real numbers needed
  real(wp), parameter :: dzero = 0.0d0
  real(wp), parameter :: zero = 0.0_wp
  real(wp), parameter :: one = 1.0_wp
  real(wp), parameter :: two = 2.0_wp
  real(wp), parameter :: three = 3.0_wp
  real(wp), parameter :: four = 4.0_wp
  real(wp), parameter :: five = 5.0_wp
  real(wp), parameter :: six = 6.0_wp
  real(wp), parameter :: seven = 7.0_wp
  real(wp), parameter :: eight = 8.0_wp
  real(wp), parameter :: nine = 9.0_wp
  real(wp), parameter :: ten = 10._wp
  real(wp), parameter :: half = 0.5_wp
  real(wp), parameter :: third = 1.0_wp/3.0_wp
  real(wp), parameter :: twothird = 2.0_wp/3.0_wp
  real(wp), parameter :: fourth = 1.0_wp/4.0_wp
  real(wp), parameter :: fifth = 1.0_wp/5.0_wp
  real(wp), parameter :: sixth = 1.0_wp/6.0_wp
  real(wp), parameter :: seventh = 1.0_wp/7.0_wp
  real(wp), parameter :: eigth = 1.0_wp/8.0_wp
  real(wp), parameter :: ninth = 1.0_wp/9.0_wp
  real(wp), parameter :: tenth = 1.0_wp/10._wp
  real(wp), parameter :: small = 1.0e-10_wp
  real(wp), parameter :: large = 1.0e05_wp
  real(wp), parameter :: realsmall = 1.d-20, big = 1.d10, realbig = 1.d20
  real(wp), parameter :: pi = 3.1415926535897932384626433832795028841971693993751058209749_wp

contains

  pure function magnitude(xin)
    real(wp), intent(in) :: xin(:)

    real(wp) :: magnitude

    magnitude = sqrt(dot_product(xin,xin))

  end function magnitude

  pure function determinant3(Pmat)
    real(wp), intent(in) :: Pmat(3,3)
    real(wp) determinant3

    determinant3 = Pmat(1,1)*(Pmat(2,2)*Pmat(3,3) &
      & - Pmat(2,3)*Pmat(3,2)) &
      & - Pmat(1,2)*(Pmat(2,1)*Pmat(3,3) &
      & - Pmat(3,1)*Pmat(2,3)) &
      & + Pmat(1,3)*(Pmat(2,1)*Pmat(3,2) &
      & - Pmat(3,1)*Pmat(2,2))

  end function determinant3 

  subroutine get_unit ( iunit )

    !*****************************************************************************80
    !
    !! GET_UNIT returns a free FORTRAN unit number.
    !
    !  Discussion:
    !
    !    A "free" FORTRAN unit number is an integer between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5, 6 and 9, which
    !    are commonly reserved for console I/O).
    !
    !    Otherwise, IUNIT is an integer between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    !  Modified:
    !
    !    18 September 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer IUNIT, the free unit number.
    !
    implicit none

    integer i
    integer ios
    integer iunit
    logical lopen

    iunit = 0

    do i = 1, 99

      if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

        inquire ( unit = i, opened = lopen, iostat = ios )

        if ( ios == 0 ) then
          if ( .not. lopen ) then
            iunit = i
            return
          end if
        end if

      end if
    end do

  end subroutine get_unit

end module precision_vars
