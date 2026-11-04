PROGRAM main_closed
!> Closed-system LAOS driver (homogeneous simple shear, no exchange of crosslinkers with the surroundings),
!> the Fortran counterpart of optimise_laos_v2.simulate / material_TC.evaluate_state(F, None, dt).
!> Each increment solves thetaf from the mass balance  cfmax*thetaf + cb_tot(thetaf) = c_tot  (RTSAFE with
!> kinetics-only calls of affclnetfic_discrete, imode = 1) and then evaluates the network stress (imode = 2).
!> In simple shear (det F = 1, phinet = 1) sigma_12 = sfic(1,2): the deviatoric projection does not change
!> off-diagonal entries and the pressure only acts on the diagonal.
!> Input (stdin):   line 1: PROPS(1:28) (as the UMAT)
!>                  line 2: gamma0, freq [Hz], dtime [s], ninc
!> Output (stdout): per increment  t, gamma, sigma_12, thetaf, cb_tot
use global
IMPLICIT NONE
INTEGER, PARAMETER :: ndi = 3, nprops = 28, maxit = 50
DOUBLE PRECISION, PARAMETER :: tol_th = 1.0d-14, tol_g = 1.0d-15
DOUBLE PRECISION :: props(nprops)
!     MODEL CONSTANTS (as _umat_.f90)
DOUBLE PRECISION :: a, r0c, etac, mu0str, beta, Lp, theta, dx, bb, lambda0, cactin, R, Rfmax, Rbmax
DOUBLE PRECISION :: CHI, D, MU0, VMOL, Koff0, Keq, Kcatch0, dxc, sfr, xi_t
DOUBLE PRECISION :: kb, b0, rgas, Mactin, rhoactin, NA, cfmax, cbmax, cabp, cb_upper, machep, tol, cb0
DOUBLE PRECISION :: x_rest, f_rest, g_rest, dum1, dum2, efi
DOUBLE PRECISION :: filprops(10), affprops(7), chemprops(12), eqprops(11), prefdir(nelem,4)
!     KINEMATICS AND NETWORK STATE
DOUBLE PRECISION :: dfgrd1(ndi,ndi), distgr(ndi,ndi), unit2(ndi,ndi), det
DOUBLE PRECISION :: sfic(ndi,ndi), cfic(ndi,ndi,ndi,ndi), dPK2ficdcb(ndi,ndi), dcbdc(ndi,ndi), dsficdthetaf(ndi,ndi)
DOUBLE PRECISION :: cb(ndir), lamn(ndir), cb_try(ndir), lamn_try(ndir)
DOUBLE PRECISION :: thetaf, cb_tot, c_tot, cbtot_try, dcbdthetaf, pnewdt
!     LOADING AND THETAF SOLVE
DOUBLE PRECISION :: gamma0, freq, dtime, time, gam, pi
DOUBLE PRECISION :: thl, thh, thx, gx, dgx, dth, dthold
INTEGER :: ninc, kinc, it, noel, nwarn, nfail

READ(*,*) props
READ(*,*) gamma0, freq, dtime, ninc

!     MATERIAL CONSTANTS (same block as _umat_.f90)
a        = props(5)
r0c      = props(6)
etac     = props(7)
mu0str   = props(8)
beta     = props(9)
Lp       = props(10)
theta    = props(11)
dx       = props(12)
bb       = props(13)
lambda0  = props(14)
cactin   = props(15)
R        = props(16)
Rfmax    = props(17)
Rbmax    = props(18)
CHI      = props(19)
D        = props(20)
MU0      = props(21)
VMOL     = props(22)
Koff0    = props(23)
Keq      = props(24)
Kcatch0  = props(25)
dxc      = props(26)
sfr      = props(27)
xi_t     = props(28)
kb = 1.380649d-5
b0 = Lp * theta * kb
rgas = 8.314462618d0
Mactin = 42.0d-3       ! [MDa]
rhoactin = 16.0d0      ! [MDa/microm]
NA = 6.022d5           ! [1/amol]
cfmax = Rfmax * cactin
cbmax = Rbmax * cactin
filprops = (/a, r0c, etac, mu0str, beta, Lp, theta, dx, kb, NA/)
affprops = (/bb, lambda0, cactin, Mactin, rhoactin, sfr, xi_t/)
eqprops = (/a, r0c, etac, mu0str, beta, b0, lambda0, Koff0, Kcatch0, dx/(kb*theta), dxc/(kb*theta)/)
!     REST EQUILIBRIUM AND C4 REST CONSTANTS
cabp = cactin * R
cb_upper = MIN(cabp, cbmax)
machep = 2.22d-16
tol = 1.0d-12
CALL pullchem(cb0, zero, cb_upper, machep, tol, cabp, cfmax, cbmax, CHI, Keq, eqprops)
CALL xmech(x_rest, f_rest, dum1, dum2, cb0, one, one, a, r0c, etac, mu0str, beta, b0, lambda0, &
           cactin*NA*Mactin/rhoactin, rgas*theta)
CALL koffcs(g_rest, dum1, f_rest, Koff0, Kcatch0, dx/(kb*theta), dxc/(kb*theta))
chemprops = (/cbmax, cfmax, CHI, D, MU0, VMOL, Koff0, Keq, Kcatch0, dxc, g_rest, x_rest/)
CALL erfi(efi, bb)
!     PREFERRED DIRECTION (as prefdir.inp: element 1, direction (1,0,0)); no uexternaldb needed
prefdir = zero
prefdir(1,1) = one
prefdir(1,2) = one
noel = 1
unit2 = zero
unit2(1,1) = one; unit2(2,2) = one; unit2(3,3) = one

!     INITIAL STATE (rest equilibrium, as initialize.f90)
cb = cb0
lamn = one
thetaf = (cabp - cb0) / cfmax
cb_tot = cb0
pi = four * ATAN(one)
nwarn = 0
nfail = 0

DO kinc = 1, ninc
  time = DBLE(kinc) * dtime
  gam = gamma0 * SIN(two * pi * freq * time)
  dfgrd1 = unit2
  dfgrd1(1,2) = gam
  CALL fslip(dfgrd1, distgr, det, ndi)
  pnewdt = one

  ! Closed system: g(thetaf) = cfmax*thetaf + cb_tot(thetaf) - c_tot = 0, conserved total content.
  ! Bracket as material_TC: g(1e-12) < 0 (only unbinding, cb_tot <= cb_tot^t < c_tot), g(c_tot/cfmax) >= 0.
  c_tot = cfmax * thetaf + cb_tot
  thl = 1.0d-12
  thh = MIN(one - 1.0d-8, c_tot / cfmax)
  thx = thetaf
  CALL gfun(thx, gx, dgx)
  dthold = thh - thl
  dth = dthold
  DO it = 1, maxit
    IF (ABS(gx) < tol_g * c_tot) EXIT
    IF (gx < zero) THEN
      thl = thx
    ELSE
      thh = thx
    END IF
    ! Newton if it stays inside the bracket and decreases fast enough, else bisection (RTSAFE)
    IF ((((thx - thh) * dgx - gx) * ((thx - thl) * dgx - gx) >= zero) .OR. &
        (ABS(two * gx) > ABS(dthold * dgx))) THEN
      dthold = dth
      dth = half * (thh - thl)
      thx = thl + dth
    ELSE
      dthold = dth
      dth = gx / dgx
      thx = thx - dth
    END IF
    CALL gfun(thx, gx, dgx)
    IF (ABS(dth) < tol_th) EXIT
  END DO
  IF (it > maxit) nfail = nfail + 1
  thetaf = thx

  ! Final evaluation at the converged thetaf: kinetics from the start-of-increment state + stress
  CALL affclnetfic_discrete(sfic, cfic, distgr, unit2, filprops, affprops, efi, noel, det, prefdir, ndi, &
       cb, dtime, chemprops, thetaf, cb_tot, dPK2ficdcb, dcbdc, pnewdt, dcbdthetaf, dsficdthetaf, lamn, 2)
  IF (pnewdt < one) nwarn = nwarn + 1

  WRITE(*,'(5ES24.15)') time, gam, sfic(1,2), thetaf, cb_tot
END DO

IF (nfail > 0) WRITE(0,*) 'WARNING: closed-system thetaf not converged in', nfail, 'increments'
IF (nwarn > 0) WRITE(0,*) 'WARNING: local kinetics not resolved (PNEWDT < 1) in', nwarn, 'increments'

CONTAINS

SUBROUTINE gfun(th, g, dg)
  ! Mass-balance residual and its derivative at thetaf = th (kinetics only, from the start-of-increment state)
  DOUBLE PRECISION, INTENT(IN)  :: th
  DOUBLE PRECISION, INTENT(OUT) :: g, dg
  cb_try = cb
  lamn_try = lamn
  CALL affclnetfic_discrete(sfic, cfic, distgr, unit2, filprops, affprops, efi, noel, det, prefdir, ndi, &
       cb_try, dtime, chemprops, th, cbtot_try, dPK2ficdcb, dcbdc, pnewdt, dcbdthetaf, dsficdthetaf, lamn_try, 1)
  g = cfmax * th + cbtot_try - c_tot
  dg = cfmax + dcbdthetaf
END SUBROUTINE gfun

END PROGRAM main_closed

SUBROUTINE getoutdir(outdir, lenoutdir)
!> Stub for the Abaqus utility (uexternaldb is linked in umat.f90 but not called by this driver)
CHARACTER(LEN=256), INTENT(OUT) :: outdir
INTEGER, INTENT(OUT) :: lenoutdir
outdir = '.'
lenoutdir = 1
END SUBROUTINE getoutdir