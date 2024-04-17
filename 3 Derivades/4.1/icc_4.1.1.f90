program icc_411

    implicit none
    integer :: i
    integer, parameter :: n = 10
    real(8), parameter :: pi = 3.14159265358979323846
    real(8) :: h, x, df_backward, df_forward, df_central
    real(8), dimension(n) :: xi, fi
    
    open(unit=1, file="4.1.1_output.dat", status="replace")

    ! Grid spacing
    h = (2.0 * pi) / real(n-1)

    ! Grid and function values
    do i = 1, n
        x = (i - 1) * h
        xi(i) = x
        fi(i) = sin(x)
    end do

  ! Calculate the derivatives using different methods
    do i = 2, n-1
        df_backward = (fi(i) - fi(i-1)) / h
        df_forward = (fi(i+1) - fi(i)) / h
        df_central = (fi(i+1) - fi(i-1)) / (2.0 * h)
        
        ! Print the results to the console
        write(1, '(4f12.6)') xi(i), df_backward, df_forward, df_central
    end do

    close(1)

end program