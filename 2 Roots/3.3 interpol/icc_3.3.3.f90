program icc_333
  implicit none

    integer, parameter :: N = 20, Nint = 100
    real(8), parameter :: L = 10.0
    real(8) :: x(N), y(N), xint(Nint), yint(Nint)
    integer :: i, j
    real(8) :: dx, dxint
    open(unit=7, file="333_original_function.dat", status="replace")
    open(unit=8, file="333_grid_points.dat", status="replace")
    open(unit=9, file="333_lagrange_interpolation.dat", status="replace")

    dx = L / real(N)
    dxint = L / real(Nint)

    ! Generate grid points and evaluate the original function
    do i = 1, N
        x(i) = real(i - 1) * dx
        y(i) = cos(x(i))
        write(7, *) x(i), y(i)
    end do

    ! Generate interpolating points using Lagrange polynomial interpolation
    do j = 1, Nint
        xint(j) = real(j - 1) * dxint
        yint(j) = lagrange_interpolate(xint(j), x, y, N)
        write(9, *) xint(j), yint(j)
    end do

    ! Write grid points to a file
    do i = 1, N
        write(8, *) x(i), y(i)
    end do

    close(7)
    close(8)
    close(9)

contains

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
