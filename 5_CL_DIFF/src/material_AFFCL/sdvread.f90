SUBROUTINE sdvread(statev,cb,cb_tot)
use global
implicit none
!>    VISCOUS DISSIPATION: READ STATE VARS
DOUBLE PRECISION, INTENT(IN)             :: statev(nsdv)
DOUBLE PRECISION, INTENT(OUT)            :: cb(ndir), cb_tot
INTEGER :: IDIR

DO IDIR = 1, ndir
    cb(IDIR) = statev(NSDV - ndir + IDIR)
END DO

cb_tot  = statev(4)

RETURN

END SUBROUTINE sdvread
