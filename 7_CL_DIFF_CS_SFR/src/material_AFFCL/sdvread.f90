SUBROUTINE sdvread(statev,thetaf,cb,cb_tot,lamn)
use global
implicit none
!>    VISCOUS DISSIPATION: READ STATE VARS
DOUBLE PRECISION, INTENT(IN)             :: statev(nsdv)
DOUBLE PRECISION, INTENT(OUT)            :: thetaf, cb(ndir), cb_tot, lamn(ndir)
INTEGER :: IDIR

DO IDIR = 1, ndir
    cb(IDIR) = statev(NSDV - 2*ndir + IDIR)
    lamn(IDIR) = statev(NSDV - ndir + IDIR)
END DO

thetaf = statev(1)
cb_tot  = statev(4)

RETURN

END SUBROUTINE sdvread
