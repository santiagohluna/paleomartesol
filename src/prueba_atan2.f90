program prueba_atan2

    use,intrinsic :: iso_fortran_env, only : dp=>real64

    implicit none

    real(dp), parameter :: pi = 4.d0*datan(1.d0)
    real(dp) :: x,y,ang,radio,ang_calc,ang_rad,atg

    ang = 0.d0

    radio = 1.d0

    do while (ang < 360.d0)

        ang_rad = ang*pi/180d0

        x = radio*dcos(ang_rad)
        y = radio*dsin(ang_rad)

        ang_calc = atg(y,x)*180.d0/pi

        print '(3(1x,F7.2))',ang,datan2(y,x)*180.d0/pi,ang_calc

        ang = ang + 10.d0

    end do
    
end program prueba_atan2

real(8) function atg(y,x) result(retval)

    use,intrinsic :: iso_fortran_env, only : dp=>real64
    
    real(dp), parameter :: pi = 4.d0*datan(1.d0)
    real(dp), intent(in) :: x,y

    retval = atan(y,x)

    if (y < 0.0d0) retval = retval + 2.d0*pi
    
end function atg