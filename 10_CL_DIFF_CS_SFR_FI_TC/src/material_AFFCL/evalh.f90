SUBROUTINE evalh(h,cb0,cabp,cfmax,cbmax,chi,Keq,eqprops)



!>     ESTABLISHMENT OF H(F)=LHS-RHS(F) THAT RELATES
!>       KEQ - CB0 RELATIONSHIP OF A SINGLE EXTENSIBLE FILAMENT
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: h
DOUBLE PRECISION, INTENT(IN)             :: cb0
DOUBLE PRECISION, INTENT(IN)             :: Keq
DOUBLE PRECISION, INTENT(IN)             :: cabp
DOUBLE PRECISION, INTENT(IN)             :: cfmax
DOUBLE PRECISION, INTENT(IN)             :: cbmax
DOUBLE PRECISION, INTENT(IN)             :: chi
DOUBLE PRECISION, INTENT(IN)             :: eqprops(11)

DOUBLE PRECISION :: lhs,rhs


DOUBLE PRECISION :: aux0,aux1,aux2,aux3,aux4
DOUBLE PRECISION :: r0f, r0, ll, lambdaf, lambda0f, mu0str, beta, b0, f0, koff_f0, dkoffdf, keq_eff

! Resting force of the prestretched filament
f0 = zero
IF (cb0 > zero) THEN
    r0f = 1.6d0 * (1.0d3 * cb0)**(-two/5.0d0)
    r0  = r0f + eqprops(2)
    ll  = eqprops(1) * r0f
    IF ((eqprops(3) > zero) .AND. (eqprops(3) .LE. one)) THEN
        lambda0f = eqprops(3) * (r0/r0f) * (eqprops(7) - one) + one
    ELSE
        lambda0f = eqprops(7)
    END IF
    lambdaf = one
    IF (lambdaf * lambda0f > one) THEN
        mu0str = eqprops(4)
        beta   = eqprops(5)
        b0     = eqprops(6)
        CALL pullforce(f0, zero, 1.0d9, 2.2204d-16, 1.0d-14, lambdaf, lambda0f, ll, r0f, mu0str, beta, b0)
    END IF
END IF
keq_eff = Keq
IF (f0 > zero) THEN
    CALL koffcs(koff_f0, dkoffdf, f0, eqprops(8), eqprops(9), eqprops(10), eqprops(11))
    keq_eff = Keq * (eqprops(8) + eqprops(9)) / koff_f0
END IF

aux0 = cabp - cb0
aux1 = aux0 / cfmax
aux2 = 1.0d0 - aux1
aux3 = exp(- chi * (1.0d0 - 2.0d0 * aux1))
    
rhs = cb0 * aux2 * aux3
lhs = keq_eff * aux0 * (1.0d0 - cb0 / cbmax)

h = lhs-rhs

RETURN
END SUBROUTINE evalh

