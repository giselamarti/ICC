program icc_422

    implicit none
    integer :: n
    real(8) :: f, rect_left, rect_right, trapezoidal, simpson, exact_integral
    real(8), parameter :: pi = 3.14159265358979323846
    real(8), parameter :: a = 0.0d0, b = 1.0d0
    external :: f

    open(unit=12, file='4.2.2_output.dat', status='replace')
    open(unit=13, file='4.2.2_convergence.dat', status='replace')
    write(12,*) "  N   ", "rect_left        ", "rect_right      ", "trapezoidal      ", "simpson"

    ! Exact integral of f(x) from 0 to 1 = 3/2 - cos(1)
    exact_integral = 3.0d0/2.0d0-cos(1.0d0)

    do n = 1, 100, 1
        if (mod(n,2).eq.0) then
            write(12,*) n, abs(pi/4.0d0 - rect_left(f,a,b,n)), abs(pi/4.0d0 - rect_right(f,a,b,n)), & 
            abs(pi/4.0d0 - trapezoidal(f,a,b,n)), abs(pi/4.0d0 - simpson(f,a,b,n))
            write(13,*) n, 1.0d0/n, 1.0d0/(n**2), 1.0d0/(n**3), 1.0d0/(n**4)
        end if
    end do
    
    close(12)
    close(13)

end program

! f(x) = x + sin(x)
real(8) function f(x)
    implicit none
    real(8) :: x

    f = sqrt(1-x**2)
    return
end function

! Rectangular left
real(8) function rect_left(f,a,b,n)
    implicit none
    integer :: i, n
    real(8) :: f, h, a, b, x

    rect_left = 0
    h = (b-a)/n

    do i = 0, n-1
        x = a + i*h
        rect_left = rect_left + h*f(x)
    end do

    return
end function

! Rectangular right
real(8) function rect_right(f,a,b,n)
    implicit none
    integer :: i, n
    real(8) :: f, h, a, b, x

    rect_right = 0
    h = (b-a)/n

    do i = 1, n
        x = a + i*h
        rect_right = rect_right + h*f(x)
    end do

    return
end function

! Trapezoidal
real(8) function trapezoidal(f,a,b,n)
    implicit none
    integer :: i, n
    real(8) :: f, h, a, b, x

    h = (b-a)/n
    trapezoidal = h*(f(a)/2.0d0+f(b)/2.0d0)

    do i = 1, n-1
        x = a + i*h
        trapezoidal = trapezoidal + h*f(x)
    end do

    return
end function

! Simpson
real(8) function simpson(f,a,b,n)
    implicit none
    integer :: i, n
    real(8) :: f, h, a, b, x

    h = (b-a)/n
    simpson = h*(f(a)/3.0d0+f(b)/3.0d0)

    do i = 1, n-1
        x = a + i*h
        if (mod(i,2).eq.0) then
            simpson = simpson + 2.0d0/3.0d0*h*f(x)
        else
            simpson = simpson + 4.0d0/3.0d0*h*f(x)
        end if
    end do

    return
end function