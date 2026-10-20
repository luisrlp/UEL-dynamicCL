SUBROUTINE initialize(statev, thetaf_t, Vmol, cb0, cfmax)
use global
IMPLICIT NONE

!      DOUBLE PRECISION TIME(2),KSTEP
INTEGER :: pos1, i
DOUBLE PRECISION, INTENT(IN)             :: thetaf_t, Vmol, cb0, cfmax
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)


pos1=1
!     VOLUME FRACTION
statev(pos1)=thetaf_t
!       DETERMINANT
statev(pos1+1)=one
!      CL CONTENT
statev(pos1+2) = thetaf_t * cfmax + cb0
!      TOTAL CB
statev(pos1+3) = cb0
!       STRESSES and CL FLUX
DO i = pos1+4, nsdv - 2*ndir
    statev(i)=zero
END DO

DO i = 1, ndir
    statev(nsdv - 2*ndir + i) = cb0 ! INITIAL CL CONCENTRATION
    statev(nsdv - ndir + i) = one  ! INITIAL REFERENCE STRETCH 
END DO
!        CONTRACTION VARIANCE
!statev(pos1+2)=zero

RETURN

END SUBROUTINE initialize
