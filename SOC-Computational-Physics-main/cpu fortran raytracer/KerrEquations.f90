module KerrEquations
    implicit none

    public
    !constants
    integer :: N, size
    real :: a, inclination, r0, theta0, a2, r_horizon, r_stable_orbit, r_disk, L, kappa
    public :: geodesic, initial

contains

    subroutine geodesic(y, dydx)
        real, intent(in) :: y(6)
        real, intent(inout) :: dydx(6)
        real :: r, theta, pr, ptheta, r2, twor, sintheta, costheta, cos2, sin2, sigma, delta, sd, siginv, bot

        r = y(1)
        theta = y(2)
        pr = y(5)
        ptheta = y(6)

        r2 = r*r
        twor = 2*r 
        sintheta = Sin(theta)
        costheta = Cos(theta)
        cos2 = costheta*costheta
        sin2 = sintheta*sintheta

        sigma = r2+a2*cos2
        delta = r2-twor+a2
        sd = sigma*delta
        siginv = 1.0/sigma
	    bot = 1.0/sd

        if (sintheta < 1e-8) then
		    sintheta = 1e-8
		    sin2 = 1e-16
        end if

        dydx(1) = -pr*delta*siginv
        dydx(2) = -ptheta*siginv
        dydx(3) = -(twor*a+(sigma-twor)*L/sin2)*bot
        dydx(4) = -(1.0+(twor*(r2+a2)-twor*a*L)*bot)
        dydx(5) = -(((r-1.0)*(-kappa)+twor*(r2+a2)-2.0*a*L)*bot-2.0*pr*pr*(r-1.0)*siginv)
        dydx(6) = -sintheta*costheta*(L*L/(sin2*sin2)-a2)*siginv

    end subroutine geodesic

    subroutine initial(y0, ydot0, x, y)
        real, intent(inout) :: y0(6), ydot0(6)
        real, intent(in) :: x, y 

        real :: sintheta, costheta, sin2, cos2, rdot0, thetadot0, phidot0, r2, sigma, delta, s1, energy, energy2

        y0(1) = r0
        y0(2) = theta0
        y0(3) = 0
        y0(4) = 0
        y0(5) = Cos(y)*Cos(x)
        y0(6) = Sin(y)/r0

        sintheta = Sin(theta0)
        costheta = Cos(theta0)
        cos2 = costheta*costheta
        sin2 = sintheta*sintheta

        rdot0 = y0(5)
	    thetadot0 = y0(6)
        phidot0 = ydot0(3)

        r2 = r0 * r0
        sigma = r2 + a2*cos2
        delta = r2 - 2.0 * r0 + a2
        s1 = sigma - 2.0 * r0

        energy2 = s1*(rdot0*rdot0/delta+thetadot0*thetadot0) + delta*sin2*phidot0*phidot0
        energy = sqrt(energy2)

        y0(5) = rdot0*sigma/(delta*energy)
        y0(6) = thetadot0*sigma/energy 

        ydot0(1) = rdot0
        ydot0(2) = thetadot0
        ydot0(3) = Cos(y)*Sin(x)/(r0*Sin(theta0))

        !Angular Momentum with E = 1
        L = ((sigma*delta*phidot0-2.0*a*r0*energy)*sin2/s1)/energy

        kappa = y0(6)*y0(6)+a2*sin2+L*L/sin2

        !Hack - make sure everything is normalized correctly by a call to geodesic 
        call geodesic(y0, ydot0)
    end subroutine initial

end module KerrEquations

