program icc_322

  implicit none

  integer, parameter :: n = 10 ! number of terms in the Taylor series
  real(8) :: e_approx

  open(20, file="3.2.2_output.txt", status="replace")

  e_approx = calc_e(n)
  write(20,*) "Approximation of e using", n, "terms:", e_approx

  close(20)

contains

  recursive function factorial(n) result(fact)
    integer, intent(in) :: n
    real(8) :: fact

    if (n == 0) then
      fact = 1.0
    else
      fact = real(n, 8) * factorial(n - 1)
    end if
  end function factorial

  function calc_e(n) result(e_approx)
    integer, intent(in) :: n
    real(8) :: e_approx
    integer :: i
    real(8) :: term

    e_approx = 0.0

    do i = 0, n
      term = 1.0d0 / real(Factorial(i), 8)
      e_approx = e_approx + term
    end do
  end function calc_e

end program

