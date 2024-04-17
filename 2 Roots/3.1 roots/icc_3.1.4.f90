program icc_31_false_position_method

    implicit none
    real(8) :: a, b, x, error, residue
    integer :: n
    integer, parameter :: max_iter = 1000
    real(8), parameter :: tol = 1.0d-8

    open(4, file='3.1.4_false_position.txt', status='replace')
    write(4,*) "False position method:"
    write(4,*) " "

    ! Initial guesses for the interval [a, b]
    a = 0.0d0
    b = 2.0d0

    ! Check that the function values have opposite signs at a and b
    if (sign(cos(a) - a, cos(b) - b) > 0.0d0) then
        write(4,*) "Initial interval [a, b] does not have opposite signs of f(a) and f(b)."
        stop
    end if

    ! Iterate using the false position method
    do n = 1, max_iter
        x = (a * (cos(b) - b) - b * (cos(a) - a)) / (cos(b) - b - cos(a) + a)
        error = abs(b - a) / 2.0d0
        residue = cos(x) - x

        ! Check for convergence
        if (error .lt. tol) then
            write(4,*) "Root found:", x
            write(4,*) "Number of iterations:", n
            write(4,*) "Final error:", error
            write(4,*) "Final residue:", residue
            exit
        end if

        ! Update the interval [a, b] based on the sign of f(x)
        if (sign(cos(a) - a, cos(x) - x) > 0.0d0) then
            a = x
            else
            b = x
        end if
    end do

    ! If max_iterations is reached, print an error message
    if (n .gt. max_iter) then
        write(3,*) "Root not found within maximum iterations."
    end if  

    close(4)

end program