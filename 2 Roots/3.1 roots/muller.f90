program icc_31_mullers_method

    implicit none
    real(8) :: x0, x1, x2, x, error, residue, h0, h1, a, b, c, discriminant
    integer :: n
    integer, parameter :: max_iter = 1000
    real(8), parameter :: tol = 1.0d-8

    open(10, file='jiuefnwoe', status='replace')

    ! Initial guesses for x0, x1, and x2
    x0 = 0.0d0
    x1 = 0.5d0
    x2 = 1.0d0

    ! Iterate using Müller's method
    do n = 1, max_iter

        h0 = x1 - x0
        h1 = x2 - x1
        a = (h1 * (cos(x1) - x1) - h0 * (cos(x2) - x2)) / (h0 * h1 * (h0 - h1))
        b = ((cos(x2) - x2) - h1 * a * (h1 + h0)) / h0
        c = cos(x2) - x2

        discriminant = b**2 - 4.0d0 * a * c
        x = x2 - (2.0d0 * c) / (b + sign(sqrt(abs(discriminant)), b))

        error = abs(x - x2)
        residue = cos(x) - x

        ! Check for convergence
        if (error < tol) then
            write(10,*) "Root found:", x
            write(10,*) "Number of iterations:", n
            write(10,*) "Final error:", error
            write(10,*) "Final residue:", residue
            exit
        end if

        ! Update x0, x1, and x2
        x0 = x1
        x1 = x2
        x2 = x
    end do

    ! If max_iterations is reached, print an error message
    if (n .gt. max_iter) then
        write(10,*) "Root not found within maximum iterations."
    end if  

    close(10)

end program
