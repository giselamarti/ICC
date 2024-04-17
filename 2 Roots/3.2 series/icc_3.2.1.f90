program icc_321
  
    implicit none
    integer, parameter :: nmax = 20
    integer, parameter :: nmax_digits = 50
    integer :: n
    real(8) :: En, En_1, e_approx, term
    real(8), dimension(nmax) :: results

    open(30, file="3.2.1_output.txt", status="replace")

    ! --- Recursive relation --- (unstable result) ---

    write(30,*) "Starting from E1=1/e, the recursive relation En=1-n*En-1 leads to an unstable result:"
    write(30,*) " "
    En_1 = 1.0d0 / exp(1.0d0)
    results(1) = En_1

    do n = 2, nmax
        En = 1.0d0 - n * En_1
        results(n) = En
        En_1 = En
    end do

    write(30,*) "          n  ", "E_n"
    do n = 1, nmax
        write(30,*) n, results(n)
    end do

    write(30,*) " "
    write(30,*) "-------------------------------"
    write(30,*) " "

    ! --- Inverse recursion --- (stable result) ---

    write(30,*) "The inverse recursion, En-1=(1-En)/n, starting from n=10 with E10=0 is stable and permits to find E1:"
    write(30,*) " "

    results(10) = 0.0d0

    do n = 9, 1, -1
        En_1 = (1.0d0 - results(n+1)) / n
        results(n) = En_1
    end do

    write(30,*) "          n  ", "E_n"
    do n = 9, 1, -1
        write(30,*) n, results(n)
    end do

    write(30,*) " "
    write(30,*) "-------------------------------"
    write(30,*) " "

    ! --- Obtain as many digits of e as possible ---

    write(30,*) "Maximum number of digits of e:"
    write(30,*) " "

    e_approx = 1.0d0
    results(1) = e_approx

    do n = 1, nmax
        term = 1.0d0 / real(n, kind=8)
        e_approx = e_approx + term
        results(n) = e_approx
    end do

    write(30,*) "n", "Approximation of e"
    do n = 1, nmax
        write(30,*) n, results(n)
    end do

    !results(20) = 0.0d0

    !do n = 9, 1, -1
    !    En_1 = (1.0d0 - results(n+1)) / n
    !    results(n) = En_1
    !end do

    !write(30,*) "          n  ", "E_n"
    !do n = 1, 20
    !    write(30,*) n, results(n)
    !end do

    close(30)


end program