PROGRAM main_sa
!> Material-point driver for the sensitivity analysis: homogeneous simple shear F = I + gamma e1 (x) e2 with a
!> piecewise loading history, open or closed system, through the complete UMAT routine MATERIAL (_umat_.f90).
!> One source for both models: built against 10_CL_DIFF_CS_SFR_FI_TC (consistent kinetics, TC) or
!> 9_CL_DIFF_CS_SFR_FI (non-consistent kinetics) by python_scripts/build_sa.sh (xmech.f90 of folder 10 is added to
!> the folder-9 build: it is only used for the outputs below).
!>
!> Systems:
!>   open   (0): mu = mu_rest, the chemical potential of the unloaded rest equilibrium of the parameter set
!>               (exchange with a bath at the initial free-crosslinker concentration);
!>   closed (1): mu is solved every increment so that the total content c_f + c_b,tot keeps its initial value
!>               c_actin*R (Newton on mu with d(content)/dmu = DSRCDMU*DTIME, bracketed bisection as fallback).
!>
!> Outputs that use the free energy (the same for both models): the affinity of each direction,
!>   A_i/RT = (mu0_f - mu0_b)/RT + ln[thetaf/(1-thetaf)] + chi(1-2 thetaf) - ln[thetab_i/(1-thetab_i)] - x_i,
!>   (mu0_f - mu0_b)/RT = ln[keq (kc+ks) cfmax / (g_rest cbmax)] + x_rest   (paper Eq. keq_mu0),
!> the reaction dissipation  sum_i w_i A_i (cb_i - cb_i^t)/dt  and the remodelling dissipation
!>   sum_i w_i n_i / lamn_i [F(lame_i) - F(1)] (lamn_i - lamn_i^t)/dt,  F(le) = f(le) lambda0 r0 le
!> (paper Eqs. remodelling_force and remodelling_dissipation, with the rates of the Backward-Euler update), where
!> w_i = rho_i * area_i are the direction weights of the network routine (sum = 1). Units: stresses in Pa,
!> contents in amol/um^3, mu in J/mol, dissipation and power in Pa/s. Stretch-free rebinding on the whole network
!> (xi_t = 1) is assumed for lame = lambar/lamn.
!>
!> Input (stdin, list-directed):
!>   line 1          : PROPS(1:28), as the UMAT
!>   line 2          : system (0 open, 1 closed), number of breakpoints, number of snapshots, maximum time step,
!>                     initial time step [s]
!>   one per breakpoint: t gamma (t increasing from t = 0, gamma = 0 at rest); gamma is linear between breakpoints
!>                     and every breakpoint is reached exactly (snapshot times must be breakpoints)
!>   next            : the snapshot times [s] (only if there are snapshots)
!> Time stepping: adaptive, controlled by PNEWDT of MATERIAL (Backward-Euler error estimate of the kinetics):
!> PNEWDT < 1 repeats the increment with dt*PNEWDT, otherwise the next step is dt*min(PNEWDT, GROWTH_MAX), at most
!> the maximum time step. One output row per accepted increment.
!>   last line       : output base name -> <base>.dat (every increment) and <base>_snap.dat (snapshots)
!> Kinetics warnings of MATERIAL go to stdout.
use global
IMPLICIT NONE
INTEGER, PARAMETER :: ndi = 3, nshr = 3, ntens = 6, nprops = 28, maxit_mu = 40
DOUBLE PRECISION, PARAMETER :: GROWTH_MAX = 1.5d0, DT_MIN = 1.0d-9, PNEWDT_CALL = 1.0d10
INTEGER, PARAMETER :: MAX_INC_FACTOR = 20
DOUBLE PRECISION :: props(nprops)
!     MATERIAL INTERFACE
DOUBLE PRECISION :: sigma(ndi,ndi), statev(nsdv), statev_t(nsdv), ddsigdde(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: dfgrd0(ndi,ndi), dfgrd1(ndi,ndi), det, time(2), dtime, predef(1), coords(3), pnewdt
DOUBLE PRECISION :: thetaf_tau, dthetafdt, dthetafdmu, rmacro, mfluid, dmdmu, dmudx(3,1), dmdj, vmol_out, cfmax_out
DOUBLE PRECISION :: dsigdmu(ndi,ndi), dsrcdmu, dsrcdj, dsrcdgu(ndi,ndi)
INTEGER :: noel, npt, kstep, kinc
COMMON /kfilp/prefdir
DOUBLE PRECISION :: prefdir(nelem,4)
!     MODEL CONSTANTS (as _umat_.f90)
DOUBLE PRECISION :: a, r0c, etac, mu0str, beta, Lp, theta, dx, bb, lambda0, cactin, R, Rfmax, Rbmax
DOUBLE PRECISION :: chi, mu0, vmol, koff0, keq, kcatch0, dxc
DOUBLE PRECISION :: kb, b0, rgas, rt, Mactin, rhoactin, NA, kactin, cfmax, cbmax, cabp, cb_upper, machep, tol
DOUBLE PRECISION :: eqprops(11), cb0, thetaf0, mu_rest, x_rest, f_rest, g_rest, dum1, dum2, ln_affinity0
!     DIRECTIONS
DOUBLE PRECISION :: m0(3,ndir), wdir(ndir), wsum
!     LOADING
INTEGER :: system, nstep, nsnap, istep, isnap, n_reject
DOUBLE PRECISION, ALLOCATABLE :: thist(:), ghist(:), tsnap(:)
DOUBLE PRECISION :: t, gamma, dt_max, dt_init, h, h_try, t_new, gamma_new, statev_acc(nsdv), f_acc(ndi,ndi), mu_acc
LOGICAL :: limited
CHARACTER(LEN=512) :: outbase
!     CLOSED SYSTEM
DOUBLE PRECISION :: mu, mu_t, c_target, content_err, dcontent_dmu, mu_new, mu_low, mu_high, tol_content, tol_mu
LOGICAL :: has_low, has_high, converged
INTEGER :: it_mu, n_unconverged, n_kinetics_fail
!     OUTPUTS
DOUBLE PRECISION :: f_t(ndi,ndi), finv(ndi,ndi), velgrad(ndi,ndi), p_mech, d_remod, d_react, exchange
DOUBLE PRECISION :: stat_cb(4), stat_lamn(4), stat_lame(4), force_max
INTEGER :: i

READ(*,*) props
READ(*,*) system, nstep, nsnap, dt_max, dt_init
ALLOCATE(thist(0:nstep), ghist(0:nstep), tsnap(MAX(nsnap,1)))
thist(0) = zero
ghist(0) = zero
DO istep = 1, nstep
  READ(*,*) thist(istep), ghist(istep)
  IF (.NOT. (thist(istep) > thist(istep-1))) THEN
    WRITE(0,*) 'ERROR: times must increase, increment', istep
    STOP 1
  END IF
END DO
IF (nsnap > 0) READ(*,*) tsnap(1:nsnap)
READ(*,'(A)') outbase

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
chi      = props(19)
mu0      = props(21)
vmol     = props(22)
koff0    = props(23)
keq      = props(24)
kcatch0  = props(25)
dxc      = props(26)
kb = 1.380649d-5
b0 = Lp * theta * kb
rgas = 8.314462618d0
rt = rgas * theta
Mactin = 42.0d-3
rhoactin = 16.0d0
NA = 6.022d5
kactin = cactin * NA * Mactin / rhoactin
cfmax = Rfmax * cactin
cbmax = Rbmax * cactin
eqprops = (/a, r0c, etac, mu0str, beta, b0, lambda0, koff0, kcatch0, dx/(kb*theta), dxc/(kb*theta)/)

!     REST EQUILIBRIUM, REST CHEMICAL POTENTIAL AND REFERENCE AFFINITY
cabp = cactin * R
cb_upper = MIN(cabp, cbmax)
machep = 2.22d-16
tol = 1.0d-12
CALL pullchem(cb0, zero, cb_upper, machep, tol, cabp, cfmax, cbmax, chi, keq, eqprops)
thetaf0 = (cabp - cb0) / cfmax
mu_rest = mu0 + rt * (LOG(thetaf0 / (one - thetaf0)) + chi * (one - two * thetaf0))   ! J = Jc = 1 at rest
CALL xmech(x_rest, f_rest, dum1, dum2, cb0, one, one, a, r0c, etac, mu0str, beta, b0, lambda0, kactin, rt)
CALL koffcs(g_rest, dum1, f_rest, koff0, kcatch0, dx/(kb*theta), dxc/(kb*theta))
ln_affinity0 = LOG(keq * (koff0 + kcatch0) * cfmax / (g_rest * cbmax)) + x_rest
c_target = cabp

!     PREFERRED DIRECTION (as prefdir.inp: element 1, direction (1,0,0)) AND DIRECTIONS OF THE NETWORK ROUTINE
prefdir = zero
prefdir(1,1) = one
prefdir(1,2) = one
CALL directions(m0, wdir)
wsum = SUM(wdir)

!     INITIAL STATE (rest equilibrium, as initialize.f90)
statev = zero
CALL initialize(statev, thetaf0, vmol, cb0, cfmax)
noel = 1; npt = 1; kstep = 1; kinc = 0
predef = zero; coords = zero; dmudx = zero; time = zero
dfgrd1 = zero
dfgrd1(1,1) = one; dfgrd1(2,2) = one; dfgrd1(3,3) = one
dfgrd0 = dfgrd1
dtime = dt_init
n_reject = 0
tol_content = 1.0d-11 * c_target
tol_mu = 1.0d-10 * rt
n_unconverged = 0
n_kinetics_fail = 0

OPEN(10, FILE=TRIM(outbase)//'.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(11, FILE=TRIM(outbase)//'_snap.dat', STATUS='REPLACE', ACTION='WRITE')
WRITE(10,'(A,I2,A,I5,A,ES23.15,A,ES23.15,A,ES23.15,A,ES23.15)') '# system ', system, '  ndir ', ndir, &
  '  cb0 ', cb0, '  thetaf0 ', thetaf0, '  mu_rest ', mu_rest, '  sum_w ', wsum
WRITE(10,'(A,ES23.15,A,ES23.15,A,ES23.15,A,ES23.15)') '# x_rest ', x_rest, '  f_rest ', f_rest, &
  '  g_rest ', g_rest, '  ln_affinity0 ', ln_affinity0
WRITE(10,'(A)') '# t gamma s11 s22 s33 s12 s13 s23 thetaf mu cb_tot content cb_mean cb_std cb_min cb_max ' // &
  'lamn_mean lamn_std lamn_min lamn_max lame_mean lame_std lame_min lame_max f_max d_remod d_react p_mech ' // &
  'exchange it_mu converged pnewdt'
WRITE(11,'(A)') '# isnap t gamma idir m0x m0y m0z w cb lambar lamn lame f g x affinity rate'

!     t = 0: rest state (one increment at F = I; nothing changes at the rest equilibrium)
t = zero
gamma = zero
mu = mu_rest
statev_t = statev
f_t = dfgrd1
CALL evaluate(mu)
it_mu = 1
converged = .TRUE.
CALL write_increment()
isnap = 1
CALL write_snapshots()

!     LOADING HISTORY (adaptive steps between the breakpoints)
statev_acc = statev
f_acc = dfgrd1
mu_acc = mu
h = dt_init
DO istep = 1, nstep
  DO WHILE (t < thist(istep) - 1.0d-12 * MAX(one, thist(istep)))
    IF (kinc + n_reject >= MAX_INC_FACTOR * (nstep + NINT(thist(nstep) / dt_max))) THEN
      WRITE(0,*) 'ERROR: more than', MAX_INC_FACTOR * (nstep + NINT(thist(nstep) / dt_max)), &
                 ' increments, stopped at t =', t
      STOP 1
    END IF
    h_try = MIN(h, dt_max, thist(istep) - t)
    IF (thist(istep) - t - h_try < 0.1d0 * h_try) h_try = thist(istep) - t   ! no sliver before the breakpoint
    limited = (h_try < h)
    t_new = t + h_try
    gamma_new = ghist(istep-1) + (ghist(istep) - ghist(istep-1)) * (t_new - thist(istep-1)) &
                / (thist(istep) - thist(istep-1))
    IF (t_new >= thist(istep)) gamma_new = ghist(istep)
    statev_t = statev_acc
    f_t = f_acc
    dfgrd0 = f_acc
    dfgrd1 = f_acc
    dfgrd1(1,2) = gamma_new
    dtime = h_try
    time(1) = t
    time(2) = t
    mu = mu_acc
    mu_t = mu_acc
    IF (system == 0) THEN
      CALL evaluate(mu_rest)
      mu = mu_rest
      it_mu = 1
      converged = .TRUE.
    ELSE
      CALL solve_closed()
    END IF
    IF (pnewdt < one) THEN                       ! rejected: repeat from the accepted state with a smaller step
      n_reject = n_reject + 1
      h = h_try * pnewdt
      statev = statev_acc
      dfgrd1 = f_acc
      mu = mu_acc
      IF (h < DT_MIN) THEN
        WRITE(0,*) 'ERROR: time step below', DT_MIN, ' at t =', t
        STOP 1
      END IF
      CYCLE
    END IF
    kinc = kinc + 1
    t = t_new
    gamma = gamma_new
    IF (.NOT. converged) n_unconverged = n_unconverged + 1
    CALL write_increment()
    CALL write_snapshots()
    statev_acc = statev
    f_acc = dfgrd1
    mu_acc = mu
    IF (limited) THEN                            ! step shortened by a breakpoint or dt_max: keep h
      h = MAX(h, h_try * MIN(pnewdt, GROWTH_MAX))
    ELSE
      h = h_try * MIN(pnewdt, GROWTH_MAX)
    END IF
  END DO
END DO
CLOSE(10)
CLOSE(11)
IF (n_unconverged > 0) WRITE(0,*) 'WARNING: closed-system mu not converged in', n_unconverged, 'increments'
WRITE(0,'(A,I8,A,I8,A)') 'increments: ', kinc, ' accepted, ', n_reject, ' rejected'

CONTAINS

SUBROUTINE evaluate(mu_in)
  ! MATERIAL at chemical potential mu_in from the state at the start of the increment
  DOUBLE PRECISION, INTENT(IN) :: mu_in
  DOUBLE PRECISION :: mu_call
  statev = statev_t
  pnewdt = PNEWDT_CALL
  mu_call = mu_in
  CALL material(sigma, statev, ddsigdde, dfgrd0, dfgrd1, det, time, dtime, predef, ndi, nshr, ntens, nsdv, &
       props, nprops, coords, pnewdt, noel, npt, kstep, kinc, mu_call, thetaf_tau, dthetafdt, dthetafdmu, &
       rmacro, mfluid, dmdmu, dmudx, dmdj, vmol_out, cfmax_out, dsigdmu, dsrcdmu, dsrcdj, dsrcdgu)
END SUBROUTINE evaluate

SUBROUTINE solve_closed()
  ! mu such that content(mu) = c_f + c_b,tot = c_target; content increases with mu
  converged = .FALSE.
  has_low = .FALSE.
  has_high = .FALSE.
  mu = mu_t
  DO it_mu = 1, maxit_mu
    CALL evaluate(mu)
    content_err = statev(3) - c_target
    IF (ABS(content_err) <= tol_content) THEN
      converged = .TRUE.
      EXIT
    END IF
    IF (content_err < zero) THEN
      mu_low = mu
      has_low = .TRUE.
    ELSE
      mu_high = mu
      has_high = .TRUE.
    END IF
    dcontent_dmu = dsrcdmu * dtime
    IF (.NOT. (dcontent_dmu > zero)) dcontent_dmu = cfmax * dthetafdmu
    mu_new = mu - content_err / dcontent_dmu
    IF (has_low .AND. has_high) THEN
      IF ((mu_new <= mu_low) .OR. (mu_new >= mu_high)) mu_new = half * (mu_low + mu_high)
    END IF
    IF (ABS(mu_new - mu) <= tol_mu) THEN
      mu = mu_new
      CALL evaluate(mu)
      converged = .TRUE.
      EXIT
    END IF
    mu = mu_new
  END DO
  IF (it_mu > maxit_mu) it_mu = maxit_mu
END SUBROUTINE solve_closed

SUBROUTINE direction_state(idir, cb_i, lamn_i, lambar_i, lame_i, force_i, g_i, x_i, affinity_i, rate_i, &
                           d_remod_i, d_react_i)
  ! State of direction idir at the end of the increment and its dissipation rates (weighted by w_i)
  INTEGER, INTENT(IN) :: idir
  DOUBLE PRECISION, INTENT(OUT) :: cb_i, lamn_i, lambar_i, lame_i, force_i, g_i, x_i, affinity_i, rate_i
  DOUBLE PRECISION, INTENT(OUT) :: d_remod_i, d_react_i
  DOUBLE PRECISION :: distgr(ndi,ndi), detf, mvec(ndi), cb_old, lamn_old, force_1, x_1, d1, d2
  DOUBLE PRECISION :: r0f, r0, n_i, thetab, fcal_e, fcal_1, fcopy(ndi,ndi)
  fcopy = dfgrd1
  CALL fslip(fcopy, distgr, detf, ndi)
  cb_i = statev(nsdv - 2*ndir + idir)
  lamn_i = statev(nsdv - ndir + idir)
  cb_old = statev_t(nsdv - 2*ndir + idir)
  lamn_old = statev_t(nsdv - ndir + idir)
  CALL deffil(lambar_i, mvec, m0(:,idir), distgr, ndi)
  lame_i = lambar_i / lamn_i
  CALL xmech(x_i, force_i, d1, d2, cb_i, lame_i, lamn_i, a, r0c, etac, mu0str, beta, b0, lambda0, kactin, rt)
  CALL xmech(x_1, force_1, d1, d2, cb_i, one, lamn_i, a, r0c, etac, mu0str, beta, b0, lambda0, kactin, rt)
  CALL koffcs(g_i, d1, force_i, koff0, kcatch0, dx/(kb*theta), dxc/(kb*theta))
  r0f = 1.6d0 * (1.0d3 * cb_i)**(-two/5.0d0)
  r0 = r0f + r0c
  n_i = kactin / (a * r0f)
  fcal_e = force_i * lambda0 * r0 * lame_i
  fcal_1 = force_1 * lambda0 * r0
  thetab = cb_i / cbmax
  affinity_i = ln_affinity0 + LOG(thetaf_tau / (one - thetaf_tau)) + chi * (one - two * thetaf_tau) &
             - LOG(thetab / (one - thetab)) - x_i
  rate_i = (cb_i - cb_old) / dtime
  d_remod_i = wdir(idir) * n_i / lamn_i * (fcal_e - fcal_1) * (lamn_i - lamn_old) / dtime
  d_react_i = wdir(idir) * rt * affinity_i * rate_i
END SUBROUTINE direction_state

SUBROUTINE write_increment()
  ! One row of <base>.dat: macroscopic state, statistics over the directions (weights w_i), dissipation and power
  DOUBLE PRECISION :: cb_i, lamn_i, lambar_i, lame_i, force_i, g_i, x_i, affinity_i, rate_i, dr_i, dc_i
  DOUBLE PRECISION :: vcb(ndir), vlamn(ndir), vlame(ndir)
  INTEGER :: idir, j
  d_remod = zero
  d_react = zero
  force_max = zero
  DO idir = 1, ndir
    CALL direction_state(idir, cb_i, lamn_i, lambar_i, lame_i, force_i, g_i, x_i, affinity_i, rate_i, dr_i, dc_i)
    vcb(idir) = cb_i
    vlamn(idir) = lamn_i
    vlame(idir) = lame_i
    force_max = MAX(force_max, force_i)
    d_remod = d_remod + dr_i
    d_react = d_react + dc_i
  END DO
  CALL statistics(vcb, stat_cb)
  CALL statistics(vlamn, stat_lamn)
  CALL statistics(vlame, stat_lame)
  ! Mechanical power sigma : L, L = (F - F^t) F^-1 / dt, and exchange rate of the content with the bath
  CALL matinv3d(dfgrd1, finv, ndi)
  velgrad = MATMUL(dfgrd1 - f_t, finv) / dtime
  p_mech = zero
  DO i = 1, ndi
    DO j = 1, ndi
      p_mech = p_mech + sigma(i,j) * velgrad(i,j)
    END DO
  END DO
  exchange = (statev(3) - statev_t(3)) / dtime
  WRITE(10,'(29ES24.15,2I4,ES12.4)') t, gamma, sigma(1,1), sigma(2,2), sigma(3,3), sigma(1,2), sigma(1,3), &
    sigma(2,3), thetaf_tau, mu, statev(4), statev(3), stat_cb, stat_lamn, stat_lame, force_max, d_remod, d_react, &
    p_mech, exchange, it_mu, MERGE(1, 0, converged), pnewdt
END SUBROUTINE write_increment

SUBROUTINE write_snapshots()
  ! All directions at the increment that ends at each requested time (snapshot times are breakpoints, which the
  ! adaptive steps reach exactly; the tolerance only covers round-off of the accumulated time)
  DOUBLE PRECISION :: cb_i, lamn_i, lambar_i, lame_i, force_i, g_i, x_i, affinity_i, rate_i, dr_i, dc_i
  INTEGER :: idir
  DO WHILE (isnap <= nsnap)
    IF (t < tsnap(isnap) - 1.0d-9 * MAX(one, tsnap(isnap))) EXIT
    DO idir = 1, ndir
      CALL direction_state(idir, cb_i, lamn_i, lambar_i, lame_i, force_i, g_i, x_i, affinity_i, rate_i, dr_i, dc_i)
      WRITE(11,'(I4,2ES24.15,I6,13ES24.15)') isnap, t, gamma, idir, m0(:,idir), wdir(idir), cb_i, lambar_i, &
        lamn_i, lame_i, force_i, g_i, x_i, affinity_i, rate_i
    END DO
    isnap = isnap + 1
  END DO
END SUBROUTINE write_snapshots

SUBROUTINE statistics(v, s)
  ! Weighted mean, standard deviation, minimum and maximum over the directions
  DOUBLE PRECISION, INTENT(IN) :: v(ndir)
  DOUBLE PRECISION, INTENT(OUT) :: s(4)
  s(1) = SUM(wdir * v) / wsum
  s(2) = SQRT(MAX(SUM(wdir * (v - s(1))**2) / wsum, zero))
  s(3) = MINVAL(v)
  s(4) = MAXVAL(v)
END SUBROUTINE statistics

SUBROUTINE directions(mdir, wd)
  ! Reference directions and weights rho_i*area_i in the order of affclnetfic_discrete (same loops)
  DOUBLE PRECISION, INTENT(OUT) :: mdir(3,ndir), wd(ndir)
  INTEGER, PARAMETER :: nfacedir = 2
  INTEGER :: ifacedir, f1, f2, f3, face, face_num, face_order_max, edge_num, point_num, node, ia, ib, ic, noel_d
  INTEGER :: f3_start(nfacedir), f3_end(nfacedir), f2_start(nfacedir)
  INTEGER, DIMENSION(3, nfacedir) :: off_a, off_b, off_c
  INTEGER, ALLOCATABLE :: edge_point(:,:), face_order(:), face_point(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: point_coord(:,:)
  DOUBLE PRECISION :: a_xyz(3), b_xyz(3), c_xyz(3), a2_xyz(3), b2_xyz(3), c2_xyz(3), node_xyz(3), area
  DOUBLE PRECISION :: ang, rho, efi, bdisp, prefdir0(3), ident(ndi,ndi)
  CALL icos_size(point_num, edge_num, face_num, face_order_max)
  ALLOCATE(point_coord(1:3,1:point_num), edge_point(1:2,1:edge_num), face_order(1:face_num), &
           face_point(1:face_order_max,1:face_num))
  CALL icos_shape(point_num, edge_num, face_num, face_order_max, point_coord, edge_point, face_order, face_point)
  f3_start(1) = 1; f3_end(1) = 3 * factor - 2
  f2_start(1) = 1
  f3_start(2) = 2; f3_end(2) = 3 * factor - 4
  f2_start(2) = 2
  off_a(:,1) = [2, -1, -1];  off_b(:,1) = [-1, 2, -1];  off_c(:,1) = [-1, -1, 2]
  off_a(:,2) = [-2, 1, 1];   off_b(:,2) = [1, -2, 1];   off_c(:,2) = [1, 1, -2]
  bdisp = bb
  CALL erfi(efi, bdisp)
  prefdir0 = prefdir(1,2:4)
  ident = zero
  ident(1,1) = one; ident(2,2) = one; ident(3,3) = one
  noel_d = 1
  node = 0
  DO face = 1, face_num / 2
    ia = face_point(1,face)
    ib = face_point(2,face)
    ic = face_point(3,face)
    a_xyz = point_coord(1:3,ia)
    b_xyz = point_coord(1:3,ib)
    c_xyz = point_coord(1:3,ic)
    DO ifacedir = 1, nfacedir
      DO f3 = f3_start(ifacedir), f3_end(ifacedir), 3
        DO f2 = f2_start(ifacedir), 3 * factor - f3 - ifacedir, 3
          f1 = 3 * factor - f3 - f2
          node = node + 1
          CALL sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2, f3, node_xyz)
          CALL sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 + off_a(1,ifacedir), f2 + off_a(2,ifacedir), &
               f3 + off_a(3,ifacedir), a2_xyz)
          CALL sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 + off_b(1,ifacedir), f2 + off_b(2,ifacedir), &
               f3 + off_b(3,ifacedir), b2_xyz)
          CALL sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1 + off_c(1,ifacedir), f2 + off_c(2,ifacedir), &
               f3 + off_c(3,ifacedir), c2_xyz)
          CALL sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz, area)
          mdir(:,node) = node_xyz
          CALL bangle(ang, ident, node_xyz, noel_d, prefdir0, ndi)
          CALL density(rho, ang, bdisp, efi)
          wd(node) = rho * area
        END DO
      END DO
    END DO
  END DO
  IF (node /= ndir) THEN
    WRITE(0,*) 'ERROR: number of directions', node, ' differs from NDIR', ndir
    STOP 1
  END IF
  DEALLOCATE(point_coord, edge_point, face_order, face_point)
END SUBROUTINE directions

END PROGRAM main_sa

SUBROUTINE getoutdir(outdir, lenoutdir)
!> Stub for the Abaqus utility (uexternaldb is linked in the material sources but not called by this driver)
CHARACTER(LEN=256), INTENT(OUT) :: outdir
INTEGER, INTENT(OUT) :: lenoutdir
outdir = '.'
lenoutdir = 1
END SUBROUTINE getoutdir
