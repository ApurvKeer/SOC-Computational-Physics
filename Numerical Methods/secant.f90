program main
    implicit none

    real :: x
    x = secant(-5.0, 5.0)

    print '(3x, "The root is ", f12.7)', x

contains

    real function f(x)
        implicit none

        real :: x

        f = x*x*x + 3*x*x + 12*x + 8

    end function f

    recursive function secant(x_0, x_1) result (sec)
        implicit none

        real :: x_0, x_1, x_2, f0, f1, f2, sec

        f0 = f(x_0)
        f1 = f(x_1)
        x_2 = x_1 - f1*(x_1 - x_0)/(f1 - f0)
        f2 = f(x_2)

        print '(3x, f12.7, 3x, f12.7, 3x, f12.7, 3x, f12.7)', x_0, x_1, x_2, f2

        if (abs(x_2 - x_1) < 0.0001) then
            sec = x_2
        else 
            sec = secant(x_1, x_2)
        end if

    end function secant
        

end program main





