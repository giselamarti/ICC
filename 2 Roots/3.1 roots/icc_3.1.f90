program icc_31

    implicit none
    integer :: i ! iteration
    integer, parameter :: max_iter = 1000 !maximum iteration
    real(8), parameter :: tol = 1.0d-8 !tolerance (epsilon)
    real(8) :: x0, x1, x, error, residue, a, b

    open(7, file="3.1_results.txt", status="replace")

    write(7,*) "Sequential method:"
    write(7,*) " "
    call sequential(x0, x, error, residue, i)
    write(7,*) " "
    write(7,*) "-------------------------------------"
    write(7,*) " "
    write(7,*) "Dichotomy method:"
    write(7,*) " "
    call dichotomy(a, b, x, error, residue, i)
    write(7,*) " "
    write(7,*) "-------------------------------------"
    write(7,*) " "
    write(7,*) "Secant method:"
    write(7,*) " "
    call secant(x0, x1, x, error, residue, i)
    write(7,*) " "
    write(7,*) "-------------------------------------"
    write(7,*) " "
    write(7,*) "False position method:"
    write(7,*) " "
    call false_position(a, b, x, error, residue, i)
    write(7,*) " "
    write(7,*) "-------------------------------------"
    write(7,*) " "
    write(7,*) "Newton's method:"
    write(7,*) " "
    call newton(x0, x, error, residue, i)
    
    close(7)

contains

    subroutine sequential(x0, x, error, residue, i)

        real(8) :: x0, x, error, residue
        integer :: i ! iteration

        ! Initial guess
        x0 = 0.5d0

        ! Iterate to find the root
        do i = 1, max_iter
            x = cos(x0)
            error = abs(x - x0)
            residue = cos(x) - x
            x0 = x

            ! Check for convergence
            if (error .lt. 1.0d-8) then
                write(7,*) "Root found:", x
                write(7,*) "Number of iterations:", i
                write(7,*) "Final error:", error
                write(7,*) "Final residue:", residue
                exit
            end if
        end do

        ! If max_iterations is reached, print an error message
        if (i .gt. max_iter) then
            write(7,*) "Root not found within maximum iterations."
        end if 

    end subroutine sequential

    subroutine dichotomy(a, b, x, error, residue, i)

        real(8) :: a, b, x, error, residue
        integer :: i

        a = 0.0d0
        b = 2.0d0

        ! Check that the function values have opposite signs at a and b
        if (sign(cos(a) - a, cos(b) - b) .gt. 0.0d0) then
            write(7,*) "Initial interval [a, b] does not have opposite signs of f(a) and f(b)."
            stop
        end if

        do i = 1, max_iter
            x = (a + b) / 2.0d0
            error = abs(b - a) / 2.0d0
            residue = cos(x) - x

            if (error .lt. tol) then
                write(7,*) "Root found:", x
                write(7,*) "Number of iterations:", i
                write(7,*) "Final error:", error
                write(7,*) "Final residue:", residue
                exit
            end if

            ! Update the interval [a, b] based on the sign of f(x)
            if (sign(cos(a) - a, cos(x) - x) .gt. 0.0d0) then
                a = x
            else
                b = x
            end if
        end do

        if (i .gt. max_iter) then
            write(7,*) "Root not found within maximum iterations."
        end if 

    end subroutine dichotomy

    subroutine secant(x0, x1, x, error, residue, i)

        real(8) :: x0, x1, x, error, residue
        integer :: i

        x0 = 0.0d0
        x1 = 2.0d0

        do i = 1, max_iter
            x = x1 - (cos(x1) - x1) * (x1 - x0) / ((cos(x1) - x1) - (cos(x0) - x0))
            error = abs(x - x1)
            residue = cos(x) - x

            if (error .lt. tol) then
                write(7,*) "Root found:", x
                write(7,*) "Number of iterations:", i
                write(7,*) "Final error:", error
                write(7,*) "Final residue:", residue
                exit
            end if

            x0 = x1
            x1 = x
        end do

        if (i .gt. max_iter) then
            write(7,*) "Root not found within maximum iterations."
        end if

    end subroutine secant

    subroutine false_position(a, b, x, error, residue, i)

        real(8) :: a, b, x, error, residue
        integer :: i

        a = 0.0d0
        b = 2.0d0

        if (sign(cos(a) - a, cos(b) - b) > 0.0d0) then
            write(7,*) "Initial interval [a, b] does not have opposite signs of f(a) and f(b)."
            stop
        end if

        do i = 1, max_iter
            x = (a * (cos(b) - b) - b * (cos(a) - a)) / (cos(b) - b - cos(a) + a)
            error = abs(b - a) / 2.0d0
            residue = cos(x) - x

            if (error .lt. tol) then
                write(7,*) "Root found:", x
                write(7,*) "Number of iterations:", i
                write(7,*) "Final error:", error
                write(7,*) "Final residue:", residue
                exit
            end if

            if (sign(cos(a) - a, cos(x) - x) > 0.0d0) then
                a = x
            else
                b = x
            end if
        end do

        if (i .gt. max_iter) then
            write(7,*) "Root not found within maximum iterations."
        end if

    end subroutine false_position

    subroutine newton(x0, x, error, residue, i)

        real(8) :: x0, x, error, residue
        integer :: i

        x0 = 0.0d0

        do i = 1, max_iter
            x = x0 - (cos(x0) - x0) / (-sin(x0) - 1.0d0)
            error = abs(x - x0)
            residue = cos(x) - x

            if (error .lt. tol) then
                write(7,*) "Root found:", x
                write(7,*) "Number of iterations:", i
                write(7,*) "Final error:", error
                write(7,*) "Final residue:", residue
                exit
            end if

            x0 = x
        end do

        if (i .gt. max_iter) then
            write(7,*) "Root not found within maximum iterations."
        end if

    end subroutine newton

end program