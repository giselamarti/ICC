program icc_31_sequential_method

    implicit none
    integer, parameter :: max_iter = 1000
    real(8) :: x0, x, error, residue
    integer :: n ! iteration

    open(1, file='3.1.1_sequential.txt', status='replace')
    write(1,*) "Sequential method:"
    write(1,*) "    "

    ! Initial guess
    x0 = 0.5d0

    ! Iterate to find the root
    do n = 1, max_iter
        x = cos(x0)
        error = abs(x - x0)
        residue = cos(x) - x
        x0 = x

        ! Check for convergence
        if (error .lt. 1.0d-8) then
            write(1,*) "Root found:", x
            write(1,*) "Number of iterations:", n
            write(1,*) "Final error:", error
            write(1,*) "Final residue:", residue
            exit
        end if
    end do

    ! If max_iterations is reached, print an error message
    if (n .gt. max_iter) then
        write(1,*) "Root not found within maximum iterations."
    end if    

    close(1)

end program