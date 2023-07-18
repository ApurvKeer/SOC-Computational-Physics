module RK_process
    use KerrEquations
    implicit none

    public :: rkstep, rkqs
    !constants

contains

    subroutine rkstep(y, dydx, h, yout, yerr)
        real, intent(in) :: y(6), dydx(6), h
        real, intent(out) :: yout(6), yerr(6)
        integer :: i

        real :: ak(N), ytemp1(N), ytemp2(N), ytemp3(N), ytemp4(N), ytemp5(N), hdydx, yi, yt

        do i = 1,N
            hdydx = h * dydx(i)
            yi = y(i)
            ytemp1(i) = yi + 0.2 * hdydx
            ytemp2(i) = yi + (3.0/40.0) * hdydx
            ytemp3(i) = yi + 0.3 * hdydx
            ytemp4(i) = yi -(11.0/54.0) * hdydx
            ytemp5(i) = yi + (1631.0/55296.0) * hdydx
            yout(i) = yi + (37.0/378.0) * hdydx
            yerr(i) = ((37.0/378.0)-(2825.0/27648.0)) * hdydx
        end do

        call geodesic(ytemp1, ak)

        do i = 1,N
            yt = h * ak(i)
            ytemp2(i) = ytemp2(1) + (9.0/40.0) * yt
            ytemp3(i) = ytemp3(i) - 0.9 * yt
            ytemp4(i) = ytemp4(i) + 2.5 * yt
            ytemp5(i) = ytemp5(i) + (175.0/512.0) * yt
        end do 

        call geodesic(ytemp2, ak)

        do i = 1,N
            yt = h * ak(i)
            ytemp3(i) = ytemp3(i) + 1.2 * yt
            ytemp4(i) = ytemp4(i) - (70.0/27.0) * yt
            ytemp5(i) = ytemp5(i) + (575.0/13824.0) * yt
            yout(i) = yout(i) + (250.0/621.0) * yt
            yerr(i) = yerr(i) + ((250.0/621.0)-(18575.0/48384.0)) * yt
        end do

        call geodesic(ytemp3, ak)

        do i = 1,N
            yt = h * ak(i)
            ytemp4(i) = ytemp4(i) + (35.0/27.0) * yt
            ytemp5(i) = ytemp5(i) + (44275.0/110592.0) * yt
            yout(i) = yout(i) + (125.0/594.0) * yt
            yerr(i) = yerr(i) + ((125.0/594.0)-(13525.0/55296.0)) * yt
        end do

        call geodesic(ytemp4, ak)

        do i = 1,N
            yt = h * ak(i)
            ytemp5(i) = ytemp5(i) + (253.0/4096.0) * yt
            yerr(i) = yerr(i) - (277.0/14336.0) * yt
        end do

        call geodesic(ytemp5, ak)

        do i = 1,N
            yt = h * ak(i)
            yout(i) = yout(i) + (512.0/1771.0) * yt
            yerr(i) = yerr(i) + ((512.0/1771.0)-0.25) * yt
        end do
    end subroutine rkstep 

    real function rkqs(y, dydx, htry, escal, yscal, hdid) 
        real, intent(inout) :: y(6), hdid
        real, intent(in) :: dydx(6)
        real, intent(in) :: htry, escal, yscal(6)
        real :: hnext

        integer :: i
        real :: errmax, htemp
        real :: h 
        real :: yerr(N), ytemp(N)
        real :: temp
        h = htry

        wloop: do while(.true.)
            call rkstep(y, dydx, h, ytemp, yerr)

            errmax = 0.0
            do i = 1, N 
                temp = abs(yerr(i)/yscal(i))
                if (temp > errmax) then
                    errmax = temp
                end if 
            end do

            errmax = errmax * escal
            if (errmax <= 1.0) then
                exit wloop
            end if 

            htemp = 0.9 * h / sqrt(sqrt(errmax))

            h = h * 0.1

            if (h >= 0.0) then 
                if (htemp > h) then
                    h = htemp
                end if
            else
                if (htemp < h) then
                    h = htemp
                end if
            end if 
        end do wloop 

        if (errmax > 1.89e-4) then
            hnext = 0.9 * h * errmax**(-0.2)
        else
            hnext = 5.0 * h
        end if 

        hdid = h

        y = ytemp 

        rkqs = hnext !hnext with different name
        return

    end function rkqs

end module RK_process

