program icc_323
  
    implicit none
    real(8) x, result, term
    integer :: i
    integer,parameter :: n = 20 ! number of terms

    write(*,*) "Enter the value of x:"
    read(*,*) x

    result = 1.0d0
    term = 1.0d0

    do i = 1, n
        term = term * x / i
        result = result + term
    end do

    write(*,*) "exp(", x, ") =", result

end program