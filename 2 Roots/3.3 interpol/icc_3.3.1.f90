program icc_331
  
    implicit none
    integer, parameter :: N = 20, Nint = 101
    real(8), parameter :: L = 10.0
    real(8) :: x(N), y(N), xint(Nint), yint(Nint)
    integer :: i, j
    real(8) :: dx, dxint
    open(unit=1, file="original_function.dat", status="replace")
    open(unit=2, file="grid_points.dat", status="replace")
    open(unit=3, file="nearest_neighbor_approximation.dat", status="replace")

    dx = L / real(N)
    dxint = L / real(Nint)

    ! Generate grid points and evaluate the original function
    do i = 1, N
        x(i) = real(i - 1) * dx
        y(i) = cos(x(i))
        write(1, *) x(i), y(i)
    end do

    ! Generate interpolating points and perform nearest-neighbor interpolation
    do j = 1, Nint
        xint(j) = real(j - 1) * dxint
        yint(j) = nearest_neighbor_interpolate(xint(j), x, y, N)
        write(3, *) xint(j), yint(j)
    end do

    ! Write grid points to a file
    do i = 1, N
        write(2, *) x(i), y(i)
    end do

    close(1)
    close(2)
    close(3)

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

end program
