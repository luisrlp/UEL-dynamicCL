SUBROUTINE initialize(statev, phi_t, Vmol)
use global
IMPLICIT NONE

!      DOUBLE PRECISION TIME(2),KSTEP
INTEGER :: pos1, i
DOUBLE PRECISION, INTENT(IN)             :: phi_t, Vmol
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)


pos1=1
!     VOLUME FRACTION
statev(pos1)=phi_t
!       DETERMINANT
statev(pos1+1)=one
!      FLUID CONTENT
statev(pos1+2) = (1.0d0 - phi_t) / (Vmol * phi_t)
!       CL RELATIVE STIFFNESS
DO i = pos1+3, nsdv
    statev(i)=zero
END DO
!        CONTRACTION VARIANCE
!statev(pos1+2)=zero

RETURN

END SUBROUTINE initialize
