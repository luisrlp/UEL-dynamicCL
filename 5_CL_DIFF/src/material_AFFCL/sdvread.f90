SUBROUTINE sdvread(statev,cb)
use global
implicit none
!>    VISCOUS DISSIPATION: READ STATE VARS
DOUBLE PRECISION, INTENT(IN)             :: statev(nsdv)
DOUBLE PRECISION, INTENT(OUT)            :: cb(ndir)
INTEGER :: IDIR

DO IDIR = 1, ndir
    cb(IDIR) = statev(NSDV - ndir + IDIR)
END DO

RETURN

END SUBROUTINE sdvread
