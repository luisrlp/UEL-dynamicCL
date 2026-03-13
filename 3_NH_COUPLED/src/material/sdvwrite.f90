SUBROUTINE sdvwrite(det,statev,sigma,phi_tau)
!>    VISCOUS DISSIPATION: WRITE STATE VARS
use global
implicit none

INTEGER :: pos1, min_idx
!
DOUBLE PRECISION, INTENT(IN)             :: det
! DOUBLE PRECISION, INTENT(IN)             :: etac_sdv(nsdv-1)
DOUBLE PRECISION, INTENT(IN)             :: sigma(6)
DOUBLE PRECISION, INTENT(IN)             :: phi_tau
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)
!
pos1=1
statev(pos1)=phi_tau
statev(pos1+1)=det

! Find out how many stress components we actually have room for
min_idx = MIN(6, nsdv - 2)

IF (min_idx > 0) THEN
    statev(3 : 2 + min_idx) = sigma(1 : min_idx)
END IF

RETURN

END SUBROUTINE sdvwrite

