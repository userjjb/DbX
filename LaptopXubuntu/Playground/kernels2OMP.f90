subroutine K6S7(t, d, x, y, z, tgt, res, n) ! output
    use omp_lib
    implicit none
    integer, intent(in) :: n, t
    real(8), intent(in) :: d
    real(8), intent(in), dimension(n) :: x, y, z
    real(8), intent(in) :: tgt(3)
    real(8), intent(out) :: res(n)
!f2py depend(n*m) res
    real(8) :: r, k
    integer :: i

    call omp_set_num_threads(t)
    !$OMP PARALLEL DO PRIVATE(r, k) SHARED(res)
    DO i = 1, n
        r = (x(i) - tgt(1))**2 + (y(i) - tgt(2))**2 + (z(i) - tgt(3))**2
        k = sqrt(d**2 + r)

        res(i) = (15*d**6)/(8*k**7) + (15*d**2*(2*d**6 - 5*d**4*r))/(16*k**9) + (15*d**4 + 20*d**2*r + 8*r**2)/(8*k**5) + &
        (5*d**3*(6*d**7 - 37*d**5*r + 20*d**3*r**2))/(16*k**11) + &
        (15*d**4*(8*d**8 - 88*d**6*r + 115*d**4*r**2 - 20*d**2*r**3))/(64*k**13) + &
        (3*d**5*(40*d**9 - 680*d**7*r + 1563*d**5*r**2 - 680*d**3*r**3 + 40*d*r**4))/(64*k**15) + &
        (5*d**6*(48*d**10 - 1160*d**8*r + 4078*d**6*r**2 - 3195*d**4*r**3 + 520*d**2*r**4 - 8*r**5))/(128*k**17) + &
        (15*d**7*(16*d**11 - 520*d**9*r + 2578*d**7*r**2 - 3143*d**5*r**3 + 980*d**3*r**4 - 56*d*r**5))/(128*k**19)
    END DO
    !$OMP END PARALLEL DO

end subroutine