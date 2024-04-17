program icc_431

    implicit none

    real(8) :: f, w, FT
    real(8), parameter :: a = -10.0d0, b = 10.0d0
    integer :: i
    external :: f

    open(unit=20, file='4.3.1_output.dat', status='replace')
    write(20,*) "  w                        ", "Fourier Transform          ", "Real function"

    do i = -100, 100
        w = i/10.0d0
        write(20,*) w, FT(f,w,a,b), f(w)
    end do

    close(20)

end program

! f(t) = exp(-t^2)
real(8) function f(t)
    implicit none
    real(8) :: t

    f = exp(-t**2)
    return
end function

! Fourier transform using Simpson's method
real(8) function FT(f,w,a,b)
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