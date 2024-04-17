program icc_31_mullers_method

    implicit none
    real(8) :: x0, x1, x2, x, h1, h2, a, b, c, d, error, residue
    integer :: n
    integer, parameter :: max_iter = 1000
    real(8), parameter :: tol = 1.0d-8

    open(6, file='3.1.6_muller.txt', status='replace')
    write(6,*) "Muller's method:"
    write(6,*) " "

    ! Initial guesses
    x0 = 0.0d0
    x1 = 1.0d0
    x2 = 2.0d0

    ! Iterate using Muller's method
    do n = 1, max_iter
        h1 = x1 - x0
        h2 = x2 - x1
        a = (h2 - h1) / (h2 * h1)
        b = (h1 + h2) / h2 - (h2 - h1) / h1 * (1.0d0 + a * h2)
        c = x2 * (1.0d0 + a * h2) * (1.0d0 + b * h1)
        d = sqrt(c * c - 4.0d0 * a * x2)

        ! Determine the better root
        if (abs(c - d) > abs(c + d)) then
            x = x2 + (-2.0d0 * x2) / (c + d)
        else
            x = x2 + (-2.0d0 * x2) / (c - d)
        end if

        error = abs(x - x2)
        residue = cos(x) - x

        ! Check for convergence
        if (error .lt. tol) then
            write(6,*) "Root found:", x
            write(6,*) "Number of iterations:", n
            write(6,*) "Final error:", error
            write(6,*) "Final residue:", residue
            exit
        end if

        ! Update the values of x0, x1, and x2
        x0 = x1
        x1 = x2
        x2 = x
     end do

    ! If max_iterations is reached, print an error message
    if (n .gt. max_iter) then
        write(6,*) "Root not found within maximum iterations."
    end if  

    close(6)

end program