program icc_23

    real, dimension(100) :: x, f
    open(unit=1, file="icc_23_data.txt", status="replace")
    
    do i = 1, 100
        x(i) = i
        f(i) = sin(x(i))/x(i)
        write(1,*) x(i), f(i)
    end do

    close(1)

end program