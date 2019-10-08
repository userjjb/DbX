subroutine K6S7(d, r, res, n, m) ! output
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: d
    real(8), intent(in) :: r(n, m)
    real(8), intent(out) :: res(n, m)
!f2py depend(n,m) res
    real(8) :: k(n, m)

    k = sqrt(d**2 + r)
    res = (15*d**6)/(8*k**7) + (15*d**2*(2*d**6 - 5*d**4*r))/(16*k**9) + (15*d**4 + 20*d**2*r + 8*r**2)/(8*k**5) + &
    (5*d**3*(6*d**7 - 37*d**5*r + 20*d**3*r**2))/(16*k**11) + &
    (15*d**4*(8*d**8 - 88*d**6*r + 115*d**4*r**2 - 20*d**2*r**3))/(64*k**13) + &
    (3*d**5*(40*d**9 - 680*d**7*r + 1563*d**5*r**2 - 680*d**3*r**3 + 40*d*r**4))/(64*k**15) + &
    (5*d**6*(48*d**10 - 1160*d**8*r + 4078*d**6*r**2 - 3195*d**4*r**3 + 520*d**2*r**4 - 8*r**5))/(128*k**17) + &
    (15*d**7*(16*d**11 - 520*d**9*r + 2578*d**7*r**2 - 3143*d**5*r**3 + 980*d**3*r**4 - 56*d*r**5))/(128*k**19)

end subroutine

subroutine K4S6(d, r, res, n, m) ! output
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: d
    real(8), intent(in) :: r(n, m)
    real(8), intent(out) :: res(n, m)
!f2py depend(n,m) res
    real(8) :: k(n, m)

    k = sqrt(d**2 + r)
    res = (3*d**4)/(2*k**5) + (3*d**2 + 2*r)/(2*k**3) + (3*d**2*(2*d**4 - 3*d**2*r))/(4*k**7) + &
    (d**4*(6*d**4 - 23*d**2*r + 6*r**2))/(4*k**9) + &
    (3*d**6*(8*d**6 - 88*d**4*r + 115*d**2*r**2 - 20*r**3))/(16*k**13) + &
    (3*d**4*(8*d**6 - 56*d**4*r + 39*d**2*r**2 - 2*r**3))/(16*k**11) + &
    (d**6*(48*d**8 - 760*d**6*r + 1590*d**4*r**2 - 585*d**2*r**3 + 20*r**4))/(32*k**15)

end subroutine

subroutine K4S5(d, r, res, n, m) ! output
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: d
    real(8), intent(in) :: r(n, m)
    real(8), intent(out) :: res(n, m)
!f2py depend(n,m) res
    real(8) :: k(n, m)

    k = sqrt(d**2 + r)
    res = (3*d**4)/(2*k**5) + (3*d**2 + 2*r)/(2*k**3) + (3*d**2*(2*d**4 - 3*d**2*r))/(4*k**7) + &
    (d**4*(6*d**4 - 23*d**2*r + 6*r**2))/(4*k**9) + (3*d**6*(8*d**6 - 88*d**4*r + 115*d**2*r**2 - 20*r**3))/(16*k**13) + &
    (3*d**4*(8*d**6 - 56*d**4*r + 39*d**2*r**2 - 2*r**3))/(16*k**11)

end subroutine

subroutine K4S4(d, r, res, n, m) ! output
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: d
    real(8), intent(in) :: r(n, m)
    real(8), intent(out) :: res(n, m)
!f2py depend(n,m) res
    real(8) :: k(n, m)

    k = sqrt(d**2 + r)
    res = (3*d**4)/(2*k**5) + (3*d**2 + 2*r)/(2*k**3) + (3*d**2*(2*d**4 - 3*d**2*r))/(4*k**7) - &
    (d**3*(-6*d**5 + 23*d**3*r - 6*d*r**2))/(4*k**9) + (3*d**4*(8*d**6 - 56*d**4*r + 39*d**2*r**2 - 2*r**3))/(16*k**11)

end subroutine

subroutine K2S6(d, r, res, n, m) ! output
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: d
    real(8), intent(in) :: r(n, m)
    real(8), intent(out) :: res(n, m)
!f2py depend(n,m) res
    real(8) :: k(n, m)

    k = sqrt(d**2 + r)
    res = d**2/k**3 + 1/k + (d**4*(2*d**2 - 3*r))/(2*k**7) + (d**2*(2*d**2 - r))/(2*k**5) + &
    (d**4*(8*d**4 - 24*d**2*r + 3*r**2))/(8*k**9) + (d**6*(8*d**4 - 40*d**2*r + 15*r**2))/(8*k**11) + &
    (d**6*(16*d**6 - 120*d**4*r + 90*d**2*r**2 - 5*r**3))/(16*k**13)

end subroutine

subroutine K2S3E(d, r, res, n, m) ! output
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: d
    real(8), intent(in) :: r(n, m)
    real(8), intent(out) :: res(n, m)
!f2py depend(n,m) res
    real(8) :: k(n, m)

    k = sqrt(d**2 + r)
    res = d**2/k**3 + 1/k + (d**2*(2*d**2 - r))/(2*k**5) - (d**3*(-2*d**3 + 3*d*r))/(2*k**7)

end subroutine

subroutine K2S3(d, r, res, n, m) ! output
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: d
    real(8), intent(in) :: r(n, m)
    real(8), intent(out) :: res(n, m)
!f2py depend(n,m) res

    res = (8*d**6 + 8*d**4*r + 7*d**2*r**2 + 2*r**3)/(2*sqrt(d**2 + r)**7)

end subroutine



