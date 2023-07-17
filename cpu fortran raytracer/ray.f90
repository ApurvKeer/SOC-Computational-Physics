module ray
    use KerrEquations
    use RK_process 
    implicit none

    public :: binarysearch, fire_ray
    !constants
    real, parameter :: PI = 4.0*atan(1.0)

contains

    subroutine binarysearch(y, dydx, hbig)
        real, intent(inout) :: y(6), dydx(6), hbig 

        real :: hsmall = 0.0
        real :: yout(N), yerr(N), hdiff, hmid
        integer :: side

        if (y(2) > PI/2.0) then
            side = 1
        else if (y(2) < PI/2.0) then
            side = -1
        else
            ! Already at the equator
            return
        end if 

        call geodesic(y,dydx)

        do while ((y(1) > r_horizon) .and. (y(1) < r0) .and. (side /= 0))
            hdiff = hbig - hsmall

            if (hdiff < 1e-7) then 
                call rkstep(y, dydx, hbig, yout, yerr)
                y = yout
                return
            end if 

            hmid = (hbig + hsmall) / 2

            call rkstep(y, dydx, hmid, yout, yerr)

            if (side * (yout(2) - PI/2.0) > 0) then
                hsmall = hmid
            else
                hbig = hmid
            end if 
        end do 
    end subroutine binarysearch 

    subroutine fire_ray(rgb, x1, y1)
        integer, intent(out) :: rgb(3)
        integer, intent(in) :: x1, y1 
        real ::  htry, escal, hdid, hnext, range

        real :: y(N), dydx(N), yscal(N), ylaststep(N)
        integer :: side, i

        htry = 0.5
        escal = 1.0e11
        hdid = 0.0
        hnext = 0.0
        range = 0.0025 * r_disk / (size - 1.0)

        call initial(y, dydx, (x1 - (size + 1.0) / 2) * range, (y1 - (size + 1.0) / 2) * range)

        do while (.true.)
            ylaststep = y

            call geodesic(y, dydx)

            do i = 1,N
                yscal(i) = abs(y(i)) + abs(dydx(i) * htry) + 1.0e-3
            end do

            if (y(2) > PI/2.0) then
                side = 1
            else if (y(2) < PI/2.0) then
                side = -1
            else 
                side = 0
            end if 

            hnext = rkqs(y, dydx, htry, escal, yscal, hdid)

            if ((y(2)-PI/2.0)*side < 0) then
                y = ylaststep

                call binarysearch(y, dydx, hdid)

                ! Did we hit the disk? 
                if ((y(1) <= r_disk) .and. (y(1) >= r_stable_orbit)) then 
                    rgb(1) = 255
                    rgb(2) = 215
                    rgb(3) = 0

                    return
                end if 
            end if 

            ! Inside the hole, or escaped to infinity */
            if ((y(1) < r_horizon) .or. (y(1) > r0)) then 
                rgb(1) = 0
                rgb(2) = 0
                rgb(3) = 0

                return
            end if 

            htry = hnext
        end do 
    end subroutine fire_ray

end module ray

