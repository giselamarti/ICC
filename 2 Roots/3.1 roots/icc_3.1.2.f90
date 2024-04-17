program icc_31_dichotomy_method

    implicit none
    real(8) :: a, b, x, error, residue
    integer :: n ! iteration
    integer, parameter :: max_iter = 1000
    real(8), parameter :: tol = 1.0d-8

    open(2, file='3.1.2_dichotomy.txt', status='replace')
    write(2,*) "Dichotomy method:"
    write(2,*) " "

    ! Initial guesses for the interval [a, b]
    a = 0.0d0
    b = 2.0d0

    ! Check that the function values have opposite signs at a and b
    if (sign(cos(a) - a, cos(b) - b) .gt. 0.0d0) then
        write(2,*) "Initial interval [a, b] does not have opposite signs of f(a) and f(b)."
        stop
    end if

    ! Iterate using the dichotomy method
    do n = 1, max_iter
        x = (a + b) / 2.0d0
        error = abs(b - a) / 2.0d0
        residue = cos(x) - x

        ! Check for convergence
        if (error .lt. tol) then
            write(2,*) "Root found:", x
            write(2,*) "Number of iterations:", n
            write(2,*) "Final error:", error
            write(2,*) "Final residue:", residue
            exit
        end if

        ! Update the interval [a, b] based on the sign of f(x)
        if (sign(cos(a) - a, cos(x) - x) .gt. 0.0d0) then
            a = x
            else
            b = x
        end if
    end do

    ! If max_iterations is reached, print an error message
    if (n .gt. max_iter) then
        write(2,*) "Root not found within maximum iterations."
    end if   

    close(2)

end program