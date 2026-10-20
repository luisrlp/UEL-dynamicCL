SUBROUTINE koffcs(koff, dkoffdf, fi, koff0, kcatch0, dxs_kT, dxc_kT)

!>    CATCH-SLIP (TWO-PATHWAY BELL-EVANS) UNBINDING RATE AND ITS FORCE DERIVATIVE
!>      koff    = kcatch0*exp(-dxc_kT*fi) + koff0*exp(dxs_kT*fi)
!>      dkoffdf = -dxc_kT*kcatch0*exp(-dxc_kT*fi) + dxs_kT*koff0*exp(dxs_kT*fi)
!>    koff0 / dxs_kT: slip pathway, kcatch0 / dxc_kT: catch pathway (reactive distance / kB*theta).
!>    kcatch0 = 0 recovers the slip-only Bell law.
!>    Exponents are clipped to [-500, 500] to avoid overflow (-ffpe-trap=overflow); the derivative
!>    of a clipped pathway is set to zero (as material.py::_koff). 500 (koff <= ~1e216 1/s, i.e.
!>    instantaneous unbinding) instead of 700: koff is later multiplied by thetab/(1-thetab) and
!>    1/(1-thetab)**2 (up to ~1e15 near cbmax) and by dt in kineticsFunc, which must not overflow.
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: koff
DOUBLE PRECISION, INTENT(OUT)            :: dkoffdf
DOUBLE PRECISION, INTENT(IN)             :: fi
DOUBLE PRECISION, INTENT(IN)             :: koff0
DOUBLE PRECISION, INTENT(IN)             :: kcatch0
DOUBLE PRECISION, INTENT(IN)             :: dxs_kT
DOUBLE PRECISION, INTENT(IN)             :: dxc_kT

DOUBLE PRECISION, PARAMETER :: argmax = 500.d0
DOUBLE PRECISION :: arg_s, arg_c, k_s, k_c, dk_s, dk_c

arg_s = dxs_kT * fi
arg_c = -dxc_kT * fi

k_s = koff0 * EXP(MIN(MAX(arg_s, -argmax), argmax))
k_c = kcatch0 * EXP(MIN(MAX(arg_c, -argmax), argmax))

dk_s = zero
dk_c = zero
IF (ABS(arg_s) < argmax) dk_s = dxs_kT * k_s
IF (ABS(arg_c) < argmax) dk_c = -dxc_kT * k_c

koff = k_s + k_c
dkoffdf = dk_s + dk_c

RETURN
END SUBROUTINE koffcs
