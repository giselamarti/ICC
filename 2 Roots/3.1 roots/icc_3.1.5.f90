program icc_31_newtons_method

    implicit none
    real(8) :: x0, x, error, residue
    integer :: n
    integer, parameter :: max_iter = 1000
    real(8), parameter :: tol = 1.0d-8

    open(5, file='3.1.5_newton.txt', status='replace')
    write(5,*) "Newton's method:"
    write(5,*) " "

    ! Initial guess
    x0 = 0.0d0

    ! Iterate using Newton's method
    do n = 1, max_iter
        x = x0 - (cos(x0) - x0) / (-sin(x0) - 1.0d0)
        error = abs(x - x0)
        residue = cos(x) - x

        ! Check for convergence
        if (error .lt. tol) then
            write(5,*) "Root found:", x
            write(5,*) "Number of iterations:", n
            write(5,*) "Final error:", error
            write(5,*) "Final residue:", residue
            exit
        end if

        ! Update the current guess
        x0 = x
    end do

    ! If max_iterations is reached, print an error message
    if (n .gt. max_iter) then
        write(3,*) "Root not found within maximum iterations."
    end if  

    close(5)

end program