SUBROUTINE vol(ssev,pv,ppv,k,det,phi_tau)

! Code converted using TO_F90 by Alan Miller
! Date: 2020-12-12  Time: 12:08:12

!>     VOLUMETRIC CONTRIBUTION :STRAIN ENERGY FUNCTION AND DERIVATIVES
use global
implicit none


DOUBLE PRECISION :: g, aux, detfe, detfs
DOUBLE PRECISION, INTENT(OUT)            :: ssev
DOUBLE PRECISION, INTENT(OUT)            :: pv
DOUBLE PRECISION, INTENT(OUT)            :: ppv
DOUBLE PRECISION, INTENT(IN)             :: k
DOUBLE PRECISION, INTENT(IN)             :: det
DOUBLE PRECISION, INTENT(IN)             :: phi_tau

! g=(one/four)*(det*det-one-two*LOG(det))

! ssev=k*g

! pv=k*(one/two)*(det-one/det)
! aux=k*(one/two)*(one+one/(det*det))
! ppv=pv+det*aux


DETFE = DET * PHI_TAU
DETFS = 1.0D0 / PHI_TAU

SSEV = 0.5D0 * k * DETFS * (DLOG(DETFE))**2

! Pressure: PV = (J^s * K * ln(J^e)) / J
PV = (DETFS * k * DLOG(DETFE)) / DET

! Derivative of Pressure w.r.t J: PPV = dp/dJ
PPV = (DETFS * k * (1.0D0 - DLOG(DETFE))) / (DET**2)

RETURN
END SUBROUTINE vol
