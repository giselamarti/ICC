program icc_413

    implicit none
    real(8) :: x, df_exact, df_backward, df_forward, df_central, error_backward, error_forward, error_central
    integer :: i
    real(8), parameter :: x_min = 0.11, x_max = 1.0, h = 0.1
    integer, parameter :: n = 50
  
    open(unit=3, file="4.1.3_output.dat", status="replace")

    do i = 1, n
        x = x_min + (i - 1) * (x_max - x_min) / (n - 1)
        
        ! Exact derivative
        df_exact = -1.0d0 / (x**2)
        
        ! Derivatives using the three methods
        df_backward = (1.0d0 / (x - h) - 1.0d0 / x) / h
        df_forward = (1.0d0 / x - 1.0d0 / (x + h)) / h
        df_central = (1.0d0 / (x + h) - 1.0d0 / (x - h)) / (2.0d0 * h)
        
        ! Errors for each method
        error_backward = abs(df_exact - df_backward)
        error_forward = abs(df_exact - df_forward)
        error_central = abs(df_exact - df_central)
        
        write(3, *) x, error_backward, error_forward, error_central
    end do

    close(3)

end program
