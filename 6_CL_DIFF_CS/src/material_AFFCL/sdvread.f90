SUBROUTINE sdvread(statev,thetaf,cb,cb_tot)
use global
implicit none
!>    VISCOUS DISSIPATION: READ STATE VARS
DOUBLE PRECISION, INTENT(IN)             :: statev(nsdv)
DOUBLE PRECISION, INTENT(OUT)            :: thetaf, cb(ndir), cb_tot
INTEGER :: IDIR

DO IDIR = 1, ndir
    cb(IDIR) = statev(NSDV - ndir + IDIR)
END DO

thetaf = statev(1)
cb_tot  = statev(4)

RETURN

END SUBROUTINE sdvread
