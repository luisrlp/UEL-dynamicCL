SUBROUTINE initialize(statev, phi_t)
use global
IMPLICIT NONE

!      DOUBLE PRECISION TIME(2),KSTEP
INTEGER :: pos1, i
DOUBLE PRECISION, INTENT(IN)             :: phi_t
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)


pos1=1
!     VOLUME FRACTION
statev(pos1)=phi_t
!       DETERMINANT
statev(pos1+1)=one
!       CL RELATIVE STIFFNESS
DO i = pos1+2, nsdv
    statev(i)=zero
END DO
!        CONTRACTION VARIANCE
!statev(pos1+2)=zero

RETURN

END SUBROUTINE initialize
