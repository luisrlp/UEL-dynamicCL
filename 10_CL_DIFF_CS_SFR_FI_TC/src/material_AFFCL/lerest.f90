SUBROUTINE lerest(le_rest, cb, N0, aratio, r0c, etac, mu0str, beta, b0, lambda0, kactin)
!> A1 (python_scripts/now-i-want2.md): elastic stretch lam_e^rest at which a direction with bound content cb
!> carries the rest direction tension N0:   n(cb) * f(lam_e^rest; cb) * lambda0 * r0(cb) * lam_e^rest = N0.
!> Solved for the force f with the explicit inverse of the beta-model force law (as in xmech):
!>   lam_f*lam_0f = aratio*rhs(f),   lam_e = 1 + (lam_f - 1)/kappa,   kappa = etac*r0/r0f,
!> i.e. f*lam_e(f) = N0/(n*lambda0*r0), increasing in f. Newton from the target force, bisection as fallback.
use global
IMPLICIT NONE
DOUBLE PRECISION, INTENT(OUT) :: le_rest
DOUBLE PRECISION, INTENT(IN)  :: cb, N0, aratio, r0c, etac, mu0str, beta, b0, lambda0, kactin
INTEGER :: it
DOUBLE PRECISION :: pi, r0f, r0, ll, n, kap, l0f, c, kk, ftarg, f, f_new, fl, fh, le, dle

pi  = four * ATAN(one)
r0f = 1.6d0 * (cb * 1.0d3)**(-two / 5.0d0)
r0  = r0f + r0c
ll  = aratio * r0f
n   = kactin / ll
kap = etac * r0 / r0f
l0f = kap * (lambda0 - one) + one
c   = ll * ll / (pi * pi * b0)
kk  = one - one / aratio
ftarg = N0 / (n * lambda0 * r0)                  ! f*lam_e at the target

f = ftarg
DO it = 1, 50
  CALL le_of_f(f, le, dle)
  IF (ABS(f * le - ftarg) <= 1.0d-14 * ftarg) THEN
    le_rest = le
    RETURN
  END IF
  f_new = f - (f * le - ftarg) / (le + f * dle)
  IF (.NOT. (f_new > zero)) EXIT
  f = f_new
END DO
fl = zero                                         ! fallback: bisection on f
fh = two * ftarg
CALL le_of_f(fh, le, dle)
DO WHILE (fh * le - ftarg <= zero)
  fh = two * fh
  CALL le_of_f(fh, le, dle)
END DO
DO it = 1, 200
  f = half * (fl + fh)
  CALL le_of_f(f, le, dle)
  IF (f * le - ftarg > zero) THEN
    fh = f
  ELSE
    fl = f
  END IF
  IF (fh - fl <= 1.0d-14 * MAX(one, fh)) EXIT
END DO
CALL le_of_f(half * (fl + fh), le, dle)
le_rest = le
RETURN
CONTAINS
  SUBROUTINE le_of_f(ff, lee, dlee)              ! lam_e(f) and dlam_e/df
    DOUBLE PRECISION, INTENT(IN)  :: ff
    DOUBLE PRECISION, INTENT(OUT) :: lee, dlee
    DOUBLE PRECISION :: u, v, dd, q, rhs, drhs
    u  = one + ff / mu0str
    v  = one + two * ff / mu0str
    dd = one + ff * c * u
    q  = kk * u**beta * v * dd**(-beta)
    rhs  = u - q
    drhs = one / mu0str - q * (beta / (mu0str * u) + two / (mu0str * v) - beta * c * (u + ff / mu0str) / dd)
    lee  = one + (aratio * rhs / l0f - one) / kap
    dlee = aratio * drhs / (l0f * kap)
  END SUBROUTINE le_of_f
END SUBROUTINE lerest

