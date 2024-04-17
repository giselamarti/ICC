program icc_412

    implicit none
    integer, parameter :: n_max = 100
    real(8) :: eps, h, x, df_exact, df_forward, df_backward, df_central
    integer :: i, n
    character(50) :: filename
    real(8), dimension(n_max) :: error_forward, error_backward, error_central
    real(8), dimension(n_max) :: linear_fit, quadratic_fit
    real(8), parameter :: x_target = 1.0
    real(8), parameter :: pi = 3.14159265358979323846

    open(unit=2, file="4.1.2_output.dat", status="replace")
  
    ! Machine epsilon
    eps = 1.
    do while(eps + 1 > 1.)
        eps = eps/2.
    end do
  
    ! Loop over different grid spacings (h) from h=1 to 0.01*epsilon
    n = 0
    h = 1.0d0
    do while (h .ge. 0.01d0 * eps)
        n = n + 1
        
        ! Exact and numerical derivatives
        x = x_target
        df_exact = cos(x)
        df_forward = (sin(x + h) - sin(x)) / h
        df_backward = (sin(x) - sin(x - h)) / h
        df_central = (sin(x + h) - sin(x - h)) / (2.0d0 * h)
        
        ! Error for each method
        error_forward(n) = abs(df_exact - df_forward)
        error_backward(n) = abs(df_exact - df_backward)
        error_central(n) = abs(df_exact - df_central)

        ! Linear and quadratic fits
        linear_fit(n) = h
        quadratic_fit(n) = h**2
        
        write(2, *) h, error_forward(n), error_backward(n), error_central(n), linear_fit(n), quadratic_fit(n)
        
        ! Reduce h for the next iteration
        h = h / 2.0d0
    end do

    close(2)

end program
