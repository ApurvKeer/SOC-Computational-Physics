program main
    implicit none

    integer :: a
    a = euler(100, 0.0, 0.0, 0.1)

contains

    real function f2(x, y)
        implicit none

        real :: x, y

        f2 = sin(x) + cos(y)
    end function f2

    integer function euler(n, x0, y0, h)
        implicit none

        integer :: n, i
        real :: x0, y0, h
        real, dimension(n) :: x, y

        x(1) = x0
        y(1) = y0

        do i = 2,n-1
            x(i+1) = x(i) + h
            y(i+1) = y(i) + h*f2(x(i), y(i))
        end do

        open (unit=1, file = 'euler.dat')
            do i = 1,n
                write(1,*) x(i), y(i)
            end do
        close(unit=1)

        euler = 0
        return

    end function euler




end program main