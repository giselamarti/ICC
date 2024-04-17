program icc_324
    implicit none

    real(8) :: a, b, c
    real(8) :: x1, x2, x, x_prev, epsilon
    integer :: max_iter, iter

    a = 1.0d0
    b = 4.0d0
    c = -1.0d0

    max_iter = 100
    epsilon = 1.0d-6

    ! Initial guess for the roots
    x1 = 0.0d0
    x2 = -b/a

    do iter = 1, max_iter
        !x_prev = x1
        x = -c / (b + x2)
        x1 = x2
        x2 = x

        if (abs(x2 - x1) < epsilon) exit
    end do

    print *, "Root 1: ", x1
    print *, "Root 2: ", x2

end program
