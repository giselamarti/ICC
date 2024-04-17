program icc_22a

    implicit none
    real*4 :: eps

    open(1, file="2.2a_output.txt")
    write(1,*) "Real*4 results:"
    eps = 1.
    do while(eps + 1 > 1.)
        eps = eps/2.
    end do

    write(1,*) "Machine epsilon (mantissa length) =", eps
    write(1,*) "Exponent lower bound (min number) =", tiny(eps)
    write(1,*) "Exponent upper bound (max number) =", huge(eps)

    close(1)

end program