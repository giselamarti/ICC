program icc_332

    implicit none
    integer, parameter :: N = 20, Nint = 100
    real(8), parameter :: L = 10.0
    real(8) :: x(N), y(N), xint(Nint), yint(Nint)
    integer :: i, j
    real(8) :: dx, dxint
    open(unit=4, file="332_original_function.dat", status="replace")
    open(unit=5, file="332_grid_points.dat", status="replace")
    open(unit=6, file="332_linear_interpolation.dat", status="replace")

    dx = L / real(N)
    dxint = L / real(Nint)

    ! Generate grid points and evaluate the original function
    do i = 1, N
        x(i) = real(i - 1) * dx
        y(i) = cos(x(i))
        write(4, *) x(i), y(i)
    end do

    ! Generate interpolating points using linear interpolation
    do j = 1, Nint
        xint(j) = real(j - 1) * dxint
        yint(j) = linear_interpolate(xint(j), x, y, N)
        write(6, *) xint(j), yint(j)
    end do

    ! Write grid points to a file
    do i = 1, N
        write(5, *) x(i), y(i)
    end do

    close(4)
    close(5)
    close(6)

contains

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

end program