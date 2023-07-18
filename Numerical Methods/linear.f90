program linear
    implicit none

    integer, parameter :: n = 3, kk = 5
    integer :: i, j, k
    real :: A(n, n), b(n), x_0(n), res(n), x_i(n), s, d(n), x_ii(n), xi, xii, x_g(n), xg ,xd

    A = transpose(reshape((/ 4,-1,-1,-2,6,1,-1,1,7 /), shape(A)))
    b = (/ 3,9,-6 /)
    x_0 = (/ 0,0,0 /)

    x_i(:) = x_0(:)

    kloop : do k = 0, kk
        iloop : do i = 1, n
            s = 0.0
            jloop : do j = 1, n 
                if (i /= j) then
                    s = s + A(i, j)*x_i(j)
                end if
            end do jloop 
            x_ii(i) = (b(i) - s) / A(i, i)
        end do iloop 

        xi = sqrt(sum(x_i**2)/n)
        xii = sqrt(sum(x_ii**2)/n)

        if (abs(xi - xii) < 1.0E-5) then
            exit kloop 
        else
            x_i(:) = x_ii(:)
        end if         
    end do kloop 

    print *, "Jacobi Method"
    print *, x_ii 

    x_g(:) = x_0(:)

    kkloop : do k = 0, kk
        d(:) = x_g(:)
        iiloop : do i = 1, n
            s = 0.0
            jjloop : do j = 1, n 
                if (i /= j) then
                    s = s + A(i, j)*x_g(j)
                end if
            end do jjloop 
            x_g(i) = (b(i) - s) / A(i, i)
        end do iiloop 

        xd = sqrt(sum(d**2)/n)
        xg = sqrt(sum(x_g**2)/n)

        if (abs(xd - xg) < 1.0E-5) then
            exit kkloop 
        end if         
    end do kkloop 

    print *, "Gauss Seidel Method"
    print *, x_g 

end program linear


