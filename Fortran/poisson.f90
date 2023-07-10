program poisson
    implicit none

    integer, parameter :: n = 100
    integer :: i, j, a
    real :: cd_grid(n,n), cd_flat(n**2), laplacian(n**2, n**2), field_flat(n**2), field_grid(n,n)

    cd_grid = 0.0

    cd_grid(n/2 + 10,n/2) = 1.0
    cd_grid(n/2 - 10,n/2) = -1.0

    cd_grid(1,:) = 1.0
    cd_grid(n,:) = -1.0

    open (unit=1, file = 'charge.txt')
        do i = 1,n
            do j = 1,n 
                write(1,*) i, j, cd_grid(i, j) 
            end do 
        end do
    close(unit=1) 

    a = flatten(cd_grid, cd_flat)

    laplacian = 0.0

    do i = 1, n
        do j = 1, n 
            if (i < n) then
                laplacian(n*(i-1) + j, n*(i) + j) = 1.0
            end if 
            if (j < n) then
                laplacian(n*(i-1) + j, n*(i-1) + j+1) = 1.0
            end if 
            if (i > 1) then
                laplacian(n*(i-1) + j, n*(i-2) + j) = 1.0
            end if 
            if (j > 1) then
                laplacian(n*(i-1) + j, n*(i-1) + j-1) = 1.0
            end if 
            laplacian(n*(i-1) + j, n*(i-1) + j) = -4.0 
        end do
    end do 

    a = solve_system(laplacian, cd_flat, field_flat)
    a = flat2grid(field_flat, field_grid)

    open (unit=2, file = 'field.txt')
        do i = 1,n
            do j = 1,n 
                write(2,*) i, j, field_grid(i, j) 
            end do 
        end do
    close(unit=2) 


contains

    integer function flatten(grid, flat)
        implicit none

        integer :: i, j
        real :: grid(n, n), flat(n**2)

        do i = 1, n 
            do j = 1, n
                flat(n*(i-1) + j) = grid(i, j)
            end do
        end do

        flatten = 0

    end function flatten 

    integer function solve_system(a_mat, b_vec, x)
        implicit none

        integer :: max_iterations, i, j 
        real :: a_mat(n**2,n**2), b_vec(n**2), x(n**2), eps, x_new(n**2), s1, s2 

        max_iterations = 1000
        eps = 1.0E-7
        x = 0.0

        iloop : do i = 1, max_iterations
            x_new = 0.0
            jloop: do j = 1, n**2
                s1 = dot_product(a_mat(j, :j), x_new(:j))
                s2 = dot_product(a_mat(j, j+1:), x(j+1:))
                x_new(j) = (b_vec(j) -s1 - s2)/a_mat(j,j)
            end do jloop
            if (abs(sqrt(sum((x - x_new)**2)/n**2)) < eps) then
                exit iloop 
            end if
            x(:) = x_new(:)
        end do iloop
        solve_system = 0

    end function solve_system

    integer function flat2grid(flat, grid)
        implicit none

        integer :: k
        real :: grid(n, n), flat(n**2)

        do k = 1, n**2 
            grid((k-1)/n + 1, mod(k-1, n) + 1) = flat(k)
        end do

        flat2grid = 0
    end function flat2grid

end program poisson
