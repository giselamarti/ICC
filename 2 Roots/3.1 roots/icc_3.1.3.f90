program icc_31_secant_method

    implicit none
    real(8) :: x0, x1, x, error, residue
    integer :: n
    integer, parameter :: max_iter = 1000
    real(8), parameter :: tol = 1.0d-8

    open(3, file='3.1.3_secant.txt', status='replace')
    write(3,*) "Secant method:"
    write(3,*) " "

    ! Initial guesses for x0 and x1
    x0 = 0.0d0
    x1 = 2.0d0

    ! Iterate using the secant method
    do n = 1, max_iter
        x = x1 - (cos(x1) - x1) * (x1 - x0) / ((cos(x1) - x1) - (cos(x0) - x0))
        error = abs(x - x1)
        residue = cos(x) - x

        ! Check for convergence
        if (error .lt. tol) then
            write(3,*) "Root found:", x
            write(3,*) "Number of iterations:", n
            write(3,*) "Final error:", error
            write(3,*) "Final residue:", residue
            exit
        end if

        ! Update x0 and x1
        x0 = x1
        x1 = x
    end do

    ! If max_iterations is reached, print an error message
    if (n .gt. max_iter) then
        write(3,*) "Root not found within maximum iterations."
    end if   

end program