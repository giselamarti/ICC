program icc_432

    implicit none

    real(8) :: f, w, FT
    real(8), parameter :: a = -20.0d0, b = 20.0d0
    integer :: i
    external :: f

    open(unit=21, file='4.3.2_output.dat', status='replace')
    write(21,*) "  w                        ", "Fourier Transform          ", "Real function"

    do i = -100, 100
        w = i/10.0d0
        write(21,*) w, FT(f,w,a,b), f(w)
    end do

    close(21)

end program

! f(t) = 1 if abs(t) < 0.5, 0 otherwise
real(8) function f(t)
    implicit none
    real(8) :: t

    if (abs(t) .lt. 0.5d0) then
        f = 1.0d0
    else
        f = 0.0d0
    end if

    return
end function

! Fourier transform using Simpson's method
real(8) function FT(f, w, a, b)
    implicit none
    integer :: i
    integer, parameter :: n = 100
    real(8) :: f, w, t, a, b, h, x

    h = (b-a)/n
    FT = h*(f(a)*cos(w*a)/3.0d0 + f(b)*cos(w*b)/3.0d0)

    do i = 1, n-1
        x = a + i*h
        if (mod(i,2).eq.0) then
            FT = FT + 2.0d0/3.0d0*h*f(x)*cos(w*x)
        else
            FT = FT + 4.0d0/3.0d0*h*f(x)*cos(w*x)
        end if
    end do

    return
end function
