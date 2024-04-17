program icc_22c

    implicit none
    real*16 :: eps

    open(3, file="2.2c_output.txt")
    write(3,*) "Real*16 results:"
    eps = 1.
    do while(eps + 1 > 1.)
        eps = eps/2.
    end do

    write(3,*) "Machine epsilon (mantissa length) =", eps
    write(3,*) "Exponent lower bound (min number) =", tiny(eps)
    write(3,*) "Exponent upper bound (max number) =", huge(eps)

    close(3)

end program