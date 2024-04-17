program icc_22b

    implicit none
    real*8 :: eps

    open(2, file="2.2b_output.txt")
    write(2,*) "Real*8 results:"
    eps = 1.
    do while(eps + 1 > 1.)
        eps = eps/2.
    end do

    write(2,*) "Machine epsilon (mantissa length) =", eps
    write(2,*) "Exponent lower bound (min number) =", tiny(eps)
    write(2,*) "Exponent upper bound (max number) =", huge(eps)

    close(2)

end program