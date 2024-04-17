program icc_33_interpolation
  
    implicit none
    integer, parameter :: N = 20, Nint = 101
    real(8), parameter :: L = 10.0
    real(8) :: x(N), y(N), xint(Nint), yint(Nint)
    integer :: i, j
    real(8) :: dx, dxint
    open(unit=10, file="original_function.dat", status="replace")
    open(unit=20, file="grid_points.dat", status="replace")
    open(unit=30, file="nearest_neighbor_interepolation.dat", status="replace")
    open(unit=40, file="linear_interpolation.dat", status="replace")
    open(unit=50, file="lagrange_interpolation.dat", status="replace")

    !dx = L / real(N)
    dx = L / real(N-1)
    !dxint = L / real(Nint)
    dxint = L / real(Nint-1)

    ! Generate grid points and evaluate the original function
    do i = 1, N
        x(i) = real(i - 1) * dx
        y(i) = cos(x(i))
        write(10, *) x(i), y(i)
    end do

    ! Generate interpolating points using nearest-neighbor interpolation
    do j = 1, Nint
        xint(j) = real(j - 1) * dxint
        yint(j) = nearest_neighbor_interpolate(xint(j), x, y, N)
        write(30, *) xint(j), yint(j)
    end do

    ! Generate interpolating points using linear interpolation
    do j = 1, Nint
        xint(j) = real(j - 1) * dxint
        yint(j) = linear_interpolate(xint(j), x, y, N)
        write(40, *) xint(j), yint(j)
    end do

    ! Generate interpolating points using Lagrange polynomial interpolation
    do j = 1, Nint
        xint(j) = real(j - 1) * dxint
        yint(j) = lagrange_interpolate(xint(j), x, y, N)
        write(50, *) xint(j), yint(j)
    end do

    ! Write grid points to a file
    do i = 1, N
        write(20, *) x(i), y(i)
    end do

    close(10)
    close(20)
    close(30)
    close(40)
    close(50)

contains

    ! Nearest-neighbor interpolation function
    real(8) function nearest_neighbor_interpolate(x, xi, yi, N)
        real(8), intent(in) :: x
        real(8), intent(in) :: xi(N), yi(N)
        integer, intent(in) :: N
        integer :: i
        real(8) :: min_dist, dist

        min_dist = abs(x - xi(1))
        nearest_neighbor_interpolate = yi(1)

        do i = 2, N
            dist = abs(x - xi(i))
            if (dist < min_dist) then
                min_dist = dist
                nearest_neighbor_interpolate = yi(i)
            end if
        end do
    end function nearest_neighbor_interpolate

    ! Linear interpolation function
    real(8) function linear_interpolate(x, xi, yi, N)
        real(8), intent(in) :: x
        real(8), intent(in) :: xi(N), yi(N)
        integer, intent(in) :: N
        integer :: i
        real(8) :: x1, x2, y1, y2

        do i = 1, N - 1
            if (x >= xi(i) .and. x <= xi(i + 1)) then
                x1 = xi(i)
                x2 = xi(i + 1)
                y1 = yi(i)
                y2 = yi(i + 1)
                linear_interpolate = y1 + (x - x1) * (y2 - y1) / (x2 - x1)
                return
            end if
        end do

        ! If x is outside the range, return 0
        linear_interpolate = 0.0
    end function linear_interpolate

    ! Lagrange polynomial interpolation function
    real(8) function lagrange_interpolate(x, xi, yi, N)
        real(8), intent(in) :: x
        real(8), intent(in) :: xi(N), yi(N)
        integer, intent(in) :: N
        integer :: i, j
        real(8) :: result, term

        result = 0.0

        do i = 1, N
            term = yi(i)
            do j = 1, N
                if (i /= j) then
                    term = term * (x - xi(j)) / (xi(i) - xi(j))
                end if
            end do
            result = result + term
        end do

        lagrange_interpolate = result
    end function lagrange_interpolate

end program
