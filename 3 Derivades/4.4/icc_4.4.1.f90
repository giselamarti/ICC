program icc_441

    implicit none

    integer :: i, j
    real(8), parameter :: pi = 3.141592653589793238462643383279502884D0
    real(8), parameter :: A = 0.1
    real(8), parameter :: w1 = 20.0 !6.0 * pi
    real(8), parameter :: w2 = 25.0 !8.0 * pi
    real(8) :: t, dt, w, dw
    real(8) :: t_min, t_max
    integer, parameter :: N = 1000
    real(8), dimension(N) :: time, signal, noise, signalnoise
    real(8), dimension(N) :: FT_signal, FT_noise, FT_signalnoise
    real(8), dimension(N) :: IFT_signal, IFT_noise, IFT_signalnoise

    open(unit=3, file="4.4.1_funcions.dat", status="replace")
    open(unit=4, file="4.4.1_FT.dat", status="replace")
    open(unit=5, file="4.4.1_InvFT.dat", status="replace")

    ! Define time parameters
    t_min = -10.0
    t_max = 10.0
    dt = (t_max - t_min) / real(N - 1)
    !dw = (60.0) / real(N - 1)

    ! Generate time array and calculate functions
    do i = 1, N
        t = t_min + (i - 1) * dt
        time(i) = t
        signal(i) = exp(-t**2)
        noise(i) = A * (cos(w1 * t) + cos(w2 * t))
        signalnoise(i) = signal(i) + noise(i)
        write(3,*) t, signal(i), noise(i), signalnoise(i)
    end do

    ! Calculate the Fourier transforms
    do j = 1, N
        w = -30.0 + (j-1) * (60.0/(N-1))
        FT_signal(j) = 0.0
        FT_noise(j) = 0.0
        FT_signalnoise(j) = 0.0
        do i = 1, N
            FT_signal(j) = FT_signal(j) + signal(i) * cos(w * time(i)) * dt
            FT_noise(j) = FT_noise(j) + noise(i) * cos(w * time(i)) * dt
            FT_signalnoise(j) = FT_signalnoise(j) + signalnoise(i) * cos(w * time(i)) * dt
        end do
        write(4,*) w, FT_signal(j), FT_noise(j), FT_signalnoise(j)
    end do

    ! Calculate the inverse Fourier transforms
    do j = 1, N
        t = time(j)
        !dw = 2.0 * pi / (N * dt)
        dw = (60.0) / real(N - 1)
        IFT_signal(j) = 0.0
        IFT_noise(j) = 0.0
        IFT_signalnoise(j) = 0.0
        do i = 1, N
            IFT_signal(j) = IFT_signal(j) + FT_signal(i) * cos(dw * (i - 1) * t) * dw
            IFT_noise(j) = IFT_noise(j) + FT_noise(i) * cos(dw * (i - 1) * t) * dw
            IFT_signalnoise(j) = IFT_signalnoise(j) + FT_signalnoise(i) * cos(dw * (i - 1) * t) * dw
        end do
        IFT_signal(j) = IFT_signal(j) / (2.0 * pi)
        IFT_noise(j) = IFT_noise(j) / (2.0 * pi)
        IFT_signalnoise(j) = IFT_signalnoise(j) / (2.0 * pi)
        write(5,*) t, IFT_signal(j), IFT_noise(j), IFT_signalnoise(j)
    end do

  close(3)
  close(4)
  close(5)

end program
