program schrodinger
    implicit none

    double precision, parameter :: hbar = 1.054571817e-34
    double precision, parameter :: e = 1.602176634e-19
    double precision, parameter :: m_e = 9.1093837015e-31
    double precision, parameter :: m_p = 1.67262192369e-27
    double precision, parameter :: epsilon_0 = 8.8541878128e-12
    double precision, parameter :: pi = 4.0 * atan(1.0)
    double precision, parameter :: mu = m_e * m_p / (m_e + m_p)
    integer, parameter :: N = 100
    double precision :: Delta_x
    double precision :: T_factor
    double precision :: V_factor
    real :: proton(0:2)

    integer :: rows(0:7*N**3 - 6*N**2-1), cols(0:7*N**3 - 6*N**2-1)
    real(kind = 8) :: data(0:7*N**3 - 6*N**2-1)
    real(kind = 8) :: e_val1, e_vec1(0:N**3 -1), e_val2, e_vec2(0:N**3 -1), eigenvector_list(0:N**3 - 1, 0:5), eigenvalue_list(0:5)

    Delta_x = 2e-9 / N
    T_factor = -(hbar ** 2) / (2.0 * mu * (Delta_x ** 2))
    V_factor = e ** 2 / (4.0 * pi * epsilon_0)
    eigenvector_list = 0.0
    eigenvalue_list = 0.0

    call index3d_to_coords((N - 1) / 2.0, (N - 1) / 2.0, (N - 1) / 2.0, Delta_x, proton(0), proton(1), proton(2))

    call gethamiltonian(rows, cols, data)

    call power_iteration(rows, cols, data, max_iterations = 1000, eigenvalue = e_val1, x = e_vec1)   
    call power_iteration(rows, cols, data, mu = e_val1, max_iterations = 1000, eigenvalue = e_val2, x = e_vec2)
    
    if (e_val1 < e_val2) then
        eigenvalue_list(0) = e_val1
        eigenvalue_list(1) = e_val2
        eigenvector_list(:, 0) = e_vec1
        eigenvector_list(:, 1) = e_vec2
    else 
        eigenvalue_list(0) = e_val2
        eigenvalue_list(1) = e_val1
        eigenvector_list(:, 0) = e_vec2
        eigenvector_list(:, 1) = e_vec1
    end if

    print *, eigenvalue_list

    call plot(eigenvector_list(:, 0), N/2)
    
contains

    integer function flat_index(i, j, k)
        implicit none
        integer :: i, j, k
        flat_index = i + j * N + k * N * N
        return
    end function flat_index

    subroutine flat_to_index3d(flat_index, i, j, k)
        implicit none
        integer, intent(in) :: flat_index
        integer, intent(out) :: i, j, k
        integer :: temp
        i = mod(flat_index, N)
        temp = flat_index / N
        j = mod(temp, N)
        k = temp / N
    end subroutine flat_to_index3d

    subroutine index3d_to_coords(i, j, k, Delta_x, x, y, z)
        implicit none
        real, intent(in) :: i, j, k
        real(kind = 8), intent(in) :: Delta_x
        real, intent(out) :: x, y, z
        real :: dx
        
        dx = real(Delta_x, 4)
        x = (i + 1) * dx
        y = (j + 1) * dx
        z = (k + 1) * dx
    end subroutine index3d_to_coords

    real function distance(i, j, k, proton)
        implicit none
        real, intent(in) :: i, j, k
        real, intent(in) :: proton(0:2)
        real :: x, y, z, d
        call index3d_to_coords(i, j, k, Delta_x, x, y, z)
        d = (x - proton(0)) ** 2 + (y - proton(1)) ** 2 + (z - proton(2)) ** 2
        distance = sqrt(d)
    end function distance

    subroutine gethamiltonian(rows, cols, data)
        implicit none
        integer, intent(out):: rows(0:7*N**3 - 6*N**2-1), cols(0:7*N**3 - 6*N**2-1)
        real(kind = 8), intent(out):: data(0:7*N**3 - 6*N**2-1)
        integer :: ii, jj, kk
        integer :: index

        index = 0
        do ii = 0, N - 1
            do jj = 0, N - 1
                do kk = 0, N - 1
                    if (ii < N - 1) then
                        rows(index) = flat_index(ii, jj, kk)
                        cols(index) = flat_index(ii + 1, jj, kk)
                        data(index) = T_factor
                        index = index + 1
                    end if
                    if (jj < N - 1) then
                        rows(index) = flat_index(ii, jj, kk)
                        cols(index) = flat_index(ii, jj + 1, kk)
                        data(index) = T_factor
                        index = index + 1
                    end if
                    if (kk < N - 1) then
                        rows(index) = flat_index(ii, jj, kk)
                        cols(index) = flat_index(ii, jj, kk + 1)
                        data(index) = T_factor
                        index = index + 1
                    end if
                    if (ii > 0) then
                        rows(index) = flat_index(ii, jj, kk)
                        cols(index) = flat_index(ii - 1, jj, kk)
                        data(index) = T_factor
                        index = index + 1
                    end if
                    if (jj > 0) then
                        rows(index) = flat_index(ii, jj, kk)
                        cols(index) = flat_index(ii, jj - 1, kk)
                        data(index) = T_factor
                        index = index + 1
                    end if
                    if (kk > 0) then
                        rows(index) = flat_index(ii, jj, kk)
                        cols(index) = flat_index(ii, jj, kk - 1)
                        data(index) = T_factor
                        index = index + 1
                    end if
                    rows(index) = flat_index(ii, jj, kk)
                    cols(index) = flat_index(ii, jj, kk)
                    data(index) = -6 * T_factor - V_factor / distance(real(ii), real(jj), real(kk), proton)
                    index = index + 1
                end do
            end do
        end do

    end subroutine gethamiltonian

    subroutine sparse_dot_vec(rows, cols, data, vec, dot)
        implicit none

        integer, intent(in) :: rows(0:7*N**3 - 6*N**2-1), cols(0:7*N**3 - 6*N**2-1)
        real(kind = 8) :: data(0:7*N**3 - 6*N**2-1), vec(0:N**3 -1) 
        real(kind = 8), intent(out) :: dot(0:N**3 -1)
        integer :: index, row_index, col_index

        dot(:) = 0.0

        !A[x,y].b[y] = c[x]
        do index = 0, 7*N**3 - 6*N**2-1
            row_index = rows(index)
            col_index = cols(index)
            dot(row_index) = dot(row_index) + vec(col_index) * data(index)
        end do 

    end subroutine sparse_dot_vec

    subroutine power_iteration(rows, cols, data, x_o, mu, max_iterations, excluded_space, eigenvalue, x)
        implicit none

        integer, intent(in) :: rows(0:7*N**3 - 6*N**2 -1), cols(0:7*N**3 - 6*N**2 -1)
        real(kind = 8) :: data(0:7*N**3 - 6*N**2 -1)
        real(kind = 8), intent(in), optional :: x_o(0:N**3 - 1)
        real(kind = 8), optional, intent(in) :: mu
        integer, intent(in) :: max_iterations
        real(kind = 8), intent(in), optional :: excluded_space(0:N**3 - 1, 0:5)

        real(kind = 8), intent(out) :: eigenvalue
        real(kind = 8), intent(out) :: x(0:N**3 - 1)

        integer :: index, iter
        real(kind = 8) :: list(0:N**3 -1) , dummy(0:N**3 -1)

        list = 0.0

        if (present(x_o)) then
            x = x_o
        else 
            call random_number(x)
        end if

        x = x/norm2(x)

        if (present(mu)) then
            do index = 0, 7*N**3 - 6*N**2-1
                if (rows(index) == cols(index)) then
                    data(index) = data(index) - mu
                end if 
            end do
        end if      

        do iter = 0, max_iterations-1
            if (present(excluded_space)) then
                do index = 0, 5
                    list(:) = excluded_space(:,index)
                    x(:) = x(:) - (dot_product(x, list)/dot_product(list,list))*list(:)
                end do
            end if

            call sparse_dot_vec(rows, cols, data, x, dummy)

            x = dummy/norm2(dummy)
        end do

        call sparse_dot_vec(rows, cols, data, x, dummy)

        eigenvalue = dot_product(x, dummy)

        if (present(mu)) then
            eigenvalue = eigenvalue + mu 
        end if   

    end subroutine power_iteration

    subroutine plot (eigenvector, slice_index)
        implicit none

        real(kind = 8), intent(in) :: eigenvector(0:N**3 -1)
        integer, intent(in) :: slice_index
        real(kind = 8) :: grid(0:N-1, 0:N-1, 0:N-1)
        integer :: i, j, k, flat_index

        do flat_index = 0, N**3 -1
            call flat_to_index3d(flat_index, i, j, k)

            grid(i, j, k) = eigenvector(flat_index)
        end do

        open (unit=1, file = 'ground.txt')
            do i = 1,N
                do j = 1,N 
                    write(1,*) i, j, grid(i-1, j-1, slice_index) 
                end do 
            end do
        close(unit=1) 

    end subroutine plot

end program schrodinger
