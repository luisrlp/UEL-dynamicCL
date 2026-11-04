SUBROUTINE vol(ssev,pv,ppv,k,det,Jc)

! Code converted using TO_F90 by Alan Miller
! Date: 2020-12-12  Time: 12:08:12

!>     VOLUMETRIC CONTRIBUTION :STRAIN ENERGY FUNCTION AND DERIVATIVES
use global
implicit none


DOUBLE PRECISION :: Je
DOUBLE PRECISION, INTENT(OUT)            :: ssev
DOUBLE PRECISION, INTENT(OUT)            :: pv
DOUBLE PRECISION, INTENT(OUT)            :: ppv
DOUBLE PRECISION, INTENT(IN)             :: k
DOUBLE PRECISION, INTENT(IN)             :: det
DOUBLE PRECISION, INTENT(IN)             :: Jc

Je = det / Jc

! Volumetric Strain Energy: Psi_vol = 0.5 * K * (ln(Je))^2
SSEV = 0.5D0 * k * (DLOG(Je))**2

! Cauchy Pressure: PV = (K * ln(Je)) / J
PV = (k * DLOG(Je)) / det

! Tangent Modulus Term (p + J*dp/dJ): PPV = K / J
PPV = k / det

RETURN
END SUBROUTINE vol
