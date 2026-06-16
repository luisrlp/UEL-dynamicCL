SUBROUTINE evalh(h,cb0,cabp,cfmax,cbmax,chi,Keq)



!>     ESTABLISHMENT OF H(F)=LHS-RHS(F) THAT RELATES
!>       KEQ - CB0 RELATIONSHIP OF A SINGLE EXNTESIBLE FILAMENT
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: h
DOUBLE PRECISION, INTENT(IN)             :: cb0
DOUBLE PRECISION, INTENT(IN)             :: Keq
DOUBLE PRECISION, INTENT(IN)             :: cabp
DOUBLE PRECISION, INTENT(IN)             :: cfmax
DOUBLE PRECISION, INTENT(IN)             :: cbmax
DOUBLE PRECISION, INTENT(IN)             :: chi

DOUBLE PRECISION :: lhs,rhs


DOUBLE PRECISION :: aux0,aux1,aux2,aux3,aux4

aux0 = cabp - cb0
aux1 = aux0 / cfmax
aux2 = 1.0d0 - aux1
aux3 = exp(- chi * (1.0d0 - 2.0d0 * aux1))
    
rhs = cb0 * aux2 * aux3
lhs = Keq * aux0 * (1.0d0 - cb0 / cbmax)

h = lhs-rhs

RETURN
END SUBROUTINE evalh
