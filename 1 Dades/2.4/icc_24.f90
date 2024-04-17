program icc_24

    implicit none
    integer, parameter :: n = 61  ! Number of data points along each axis
    real(8) :: x, y, f, dx, dy
    integer :: i, j
    open(unit=10, file='icc_2.4_data.txt', status='replace')

    dx = 6.0d0 / (n-1)
    dy = 6.0d0 / (n-1)

    do i = 1, n
        x = -3.0d0 + (i-1) * dx
        do j = 1, n
            y = -3.0d0 + (j-1) * dy
            f = exp(- (x*x + y*y))
            write(10, *) x, y, f
        end do
    end do

    close(10)

end program