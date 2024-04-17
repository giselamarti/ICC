program icc_221c

    real*16 :: eps
    eps = 1.
    do while(eps + 1 > 1.)
        eps = eps/2.
    end do
    print*, "eps =",eps

end program