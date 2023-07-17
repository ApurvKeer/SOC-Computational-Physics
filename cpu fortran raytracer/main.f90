program main
    use KerrEquations
    use RK_process
    use ray
    implicit none

    integer, parameter :: width = 300 
    integer, parameter :: height = 400
    integer :: x, y, color(3), grid(height*width,3)

    N = 6
    size = 300
    a = 0.99
    inclination = 85
    r0 = 1000.0
    theta0 = (PI/180.0) * inclination
    a2 = a*a
    r_horizon = 1.0 + sqrt(1.0-a2) + 1.0e-5
	r_disk = 20.0 !times schwarzschild radius
 	r_stable_orbit = inner_orbit()

    do y = height-1,0,-1
        do x = 0, width-1
            call fire_ray(color, x, y)
        end do
    end do

    open (unit=1, file = 'plane_accretion.ppm')
        write "(A)", "P3"
        write "(i3, i4)", height, width 
        write "(i3)", 255
        do y = height-1,0,-1
            do x = 0, width-1
                write "(i3, i4, i4)", color(1), color(2), color(3)
            end do
        end do
    close(unit=1)


contains

    real function inner_orbit()
        real :: z1, z2 
        z1 = 1+cbrt(1-a2)*(cbrt(1+a)+cbrt(1-a))
	    z2 = sqrt(3*a2+z1*z1)
	    inner_orbit = 3+z2-sqrt((3-z1)*(3+z1+2*z2))
        return 
    end function inner_orbit

    real function cbrt(a)
        real :: a
        cbrt = a**(1.0/3.0)
    end function cbrt 

end program main 