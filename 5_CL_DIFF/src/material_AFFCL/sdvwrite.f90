SUBROUTINE sdvwrite(det,statev,sigma,phi_tau,dmudx,Vmol,jfluid)
!>    VISCOUS DISSIPATION: WRITE STATE VARS
use global
implicit none

INTEGER :: pos1, min_idx
!
DOUBLE PRECISION, INTENT(IN)             :: det
! DOUBLE PRECISION, INTENT(IN)             :: etac_sdv(nsdv-1)
DOUBLE PRECISION, INTENT(IN)             :: sigma(6)
DOUBLE PRECISION, INTENT(IN)             :: phi_tau
DOUBLE PRECISION, INTENT(IN)             :: dmudx(3,1), JFLUID(3,1)
DOUBLE PRECISION, INTENT(IN)             :: Vmol
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)
!
pos1=1
statev(pos1)=phi_tau
statev(pos1+1)=det
! Add fluid content cR to statev
! cR
! statev(pos1+2) = (1.0d0 - phi_tau) / (Vmol * phi_tau)
! c
statev(pos1+2) = (1.0d0 - phi_tau) / (Vmol * phi_tau * det)

! Find out how many stress components we actually have room for
min_idx = MIN(6, nsdv - 3)

IF (min_idx > 0) THEN
    statev(4 : 3 + min_idx) = sigma(1 : min_idx)
END IF

statev(pos1+3+min_idx : pos1+5+min_idx) = - dmudx(1:3,1)

statev(pos1+6+min_idx : nsdv) = JFLUID(1:3,1)

! write(*,*) 'nsdv = ', nsdv
! write(*,*) 'statev = ', statev

RETURN

END SUBROUTINE sdvwrite

