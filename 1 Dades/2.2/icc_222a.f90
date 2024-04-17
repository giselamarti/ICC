program icc_222a

    real*4 :: eps
    eps = 1.
    do while(eps + 1 > 1.)
        eps = eps/2.
    end do
    print*, "eps =", tiny(eps)

end program