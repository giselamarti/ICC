program icc_45

    implicit none

    integer, parameter :: iter = 21
    integer, parameter :: n_iter = 200
    integer :: i, j, k_index
    real(8), dimension(iter) :: x_old, x_new
    integer, dimension(iter) :: iterations
    real(8) :: x_current
    real(8), dimension(7) :: k_values = (/0.5d0, 1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, 3.5d0/)
    real(8), dimension(4) :: k2_val = (/2.9d0, 3.5d0, 3.56d0, 3.57d0/)
    real(8), parameter :: k1_min = 0.0d0
    real(8), parameter :: k1_max = 4.0d0
    real(8), parameter :: k2_min = 2.8d0
    real(8), parameter :: k3_min = 3.84d0
    real(8), parameter :: k3_max = 3.856d0
    real(8) :: k1, k2, k3, x1, x2, x3

    open(unit=1, file="4.5_population_evolution.dat", status="replace")
    open(unit=2, file="4.5_poincare_map.dat", status="replace")
    open(unit=3, file="4.5_bifurcation1.dat", status="replace")
    open(unit=4, file="4.5_bifurcation2.dat", status="replace")
    open(unit=5, file="4.5_bifurcation3.dat", status="replace")

    ! Population evolution
    do j = 1, size(k_values)
        call population_evolution(0.5d0, k_values(j), x_old, x_new, iterations)
        write(1, *) "k =", k_values(j)
        write(1, *) "Iteration     x_old"
        do i = 1, iter
            write(1, *) iterations(i), x_old(i)
        end do
    end do

    ! Pointcare map
    do j = 1, size(k2_val)
        call population_evolution(0.5d0, k2_val(j), x_old, x_new, iterations)

        write(2, *) "k =", k2_val(j)
        write(2, *) "x(t)", "x(t+1)"

        do i = 1, iter
            x_current = x_old(i)
            if (x_current >= 0.d0 .and. x_current <= 1.d0) then
                write(2, *) x_current, x_new(i)
            end if
        end do
    end do

    ! Bifurcation diagrams (one output file for each diagram)

    ! k from 0 to 4
    do i = 0, 100
        k1 = k1_min + (k1_max - k1_min) * real(i) / 100.0d0
        x1 = 0.5d0 ! Initial value of x
        do j = 1, n_iter
            x1 = k1 * x1 * (1.0d0 - x1)
            if (j > 100) then
                write(3,*) k1, x1
            end if
        end do
    end do

    ! k from 2.8 to 4
    do i = 0, 100
        k2 = k2_min + (k1_max - k2_min) * real(i) / 100.0d0
        x2 = 0.5d0 ! Initial value of x
        do j = 1, n_iter
            x2 = k2 * x2 * (1.0d0 - x2)
            if (j > 100) then
                write(4,*) k2, x2
            end if
        end do
    end do

    ! k from 3.84 to 3.856
    do i = 0, 100
        k3 = k3_min + (k3_max - k3_min) * real(i) / 100.0d0
        x3 = 0.5d0 ! Initial value of x
        do j = 1, n_iter
            x3 = k3 * x3 * (1.0d0 - x3)
            if (j > 100) then
                write(5,*) k3, x3
            end if
        end do
    end do

    close(1)
    close(2)
    close(3)
    close(3)
    close(5)

contains

    subroutine population_evolution(xo, k, x_old, x_new, iterations)
        real(8), intent(in) :: xo, k
        real(8), dimension(:), intent(out) :: x_old, x_new
        integer, dimension(:), intent(out) :: iterations
        integer :: i
        real(8) :: xn
        real(8) :: xo_local ! local variable to update xo

        xo_local = xo
        do i = 1, size(x_old)
            xn = k * xo_local * (1.0d0 - xo_local)
            x_old(i) = xo_local
            x_new(i) = xn
            iterations(i) = i - 1
            xo_local = xn  ! update xo_local
        end do
    end subroutine population_evolution

end program