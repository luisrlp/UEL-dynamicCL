
    SUBROUTINE MATERIAL(SIGMA,STATEV,DDSIGDDE,DFGRD0,DFGRD1,DET, &
    TIME,DTIME,PREDEF,NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS, &
    PNEWDT,NOEL,NPT,KSTEP,KINC,MU_TAU,THETAF_TAU,DTHETAFDT, &
      DTHETAFDMU,RMACRO,MFLUID,DMDMU,DMUDX,DMDJ,VMOL,CFMAX,DSIGDMU, &
      DSRCDMU,DSRCDJ,DSRCDGU)
!
use global
IMPLICIT NONE
!----------------------------------------------------------------------
!--------------------------- DECLARATIONS -----------------------------
!----------------------------------------------------------------------
INTEGER :: NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT, &
            LAYER, KSPT, KSTEP, KINC

INTEGER, PARAMETER :: nargs = 11

REAL(KIND=8) :: STRESS(NTENS), STATEV(NSTATEV), &
                DDSDDE(NTENS,NTENS), DDSDDT(NTENS), DRPLDE(NTENS), &
                STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1), DPRED(1), &
                PROPS(NPROPS), COORDS(3), DROT(3,3), DFGRD0(3,3), DFGRD1(3,3), &
                FIBORI(NELEM,4), ARGS(NARGS)

REAL(8), INTENT(IN)      :: MU_TAU, DMUDX(3,1)
REAL(8), INTENT(OUT)     :: DSIGDMU(NDI,NDI)
REAL(8), INTENT(OUT)     :: THETAF_TAU, DTHETAFDT, DTHETAFDMU, RMACRO! DPHIDMU, DPHIDOTDMU
REAL(8), INTENT(OUT)     :: MFLUID, DMDMU, DMDJ, VMOL, CFMAX
! Derivatives of the source s = -(cfmax*dthetaf/dt + dcb_tot/dt) (sign: -ds/d(.)), Part C item 1:
!   DSRCDMU = (cfmax+G)/dt * dthetaf/dmu,  DSRCDJ = (cfmax+G)/dt * dthetaf/dJbar,
!   DSRCDGU = -dcbdc/dt (spatial, isochoric: contracts with grad(du) at the local point)
REAL(8), INTENT(OUT)     :: DSRCDMU, DSRCDJ, DSRCDGU(NDI,NDI)

REAL(8) :: THETAF_T
! DIFFUSION VARIABLES
REAL(8) :: CHI, D, MU0, RGAS
REAL(8) :: PHI_PER, PHI_M, dPdt_per, dPdt_m, DELTAMU, JFLUID(3,1)
REAL(8) :: DphiDJ, DmDphi

REAL(KIND=8) :: SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP, &
                DTEMP, PNEWDT, CELENT

COMMON /kfilp/prefdir
COMMON /kfile/etadir
DOUBLE PRECISION :: prefdir(nelem,4)
DOUBLE PRECISION :: etadir(nelem*ngp, ndir+2)
DOUBLE PRECISION :: etadir_array(ndir)

!
!     FLAGS
!      INTEGER FLAG1
!     UTILITY TENSORS
DOUBLE PRECISION :: unit2(ndi,ndi),unit4(ndi,ndi,ndi,ndi),  &
    unit4s(ndi,ndi,ndi,ndi), proje(ndi,ndi,ndi,ndi),projl(ndi,ndi,ndi,ndi)
!     KINEMATICS
DOUBLE PRECISION :: distgr(ndi,ndi),c(ndi,ndi),b(ndi,ndi),  &
    cbar(ndi,ndi),bbar(ndi,ndi),distgrinv(ndi,ndi),  &
    ubar(ndi,ndi),vbar(ndi,ndi),rot(ndi,ndi), dfgrd1inv(ndi,ndi), &
    cinv(ndi,ndi), cbarinv(ndi,ndi)
DOUBLE PRECISION :: det,detfe, detfs,cbari1,cbari2
!     VOLUMETRIC CONTRIBUTION
DOUBLE PRECISION :: pkvol(ndi,ndi),svol(ndi,ndi),  &
    cvol(ndi,ndi,ndi,ndi),cmvol(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: k,pv,ppv,ssev,ppv_chem
!     ISOCHORIC CONTRIBUTION
DOUBLE PRECISION :: siso(ndi,ndi),pkiso(ndi,ndi),pk2(ndi,ndi),  &
    ciso(ndi,ndi,ndi,ndi),cmiso(ndi,ndi,ndi,ndi),  &
    sfic(ndi,ndi),cfic(ndi,ndi,ndi,ndi), pkfic(ndi,ndi),cmfic(ndi,ndi,ndi,ndi)
!     ISOCHORIC ISOTROPIC CONTRIBUTION
DOUBLE PRECISION :: c10,c01,sseiso,diso(5),pkmatfic(ndi,ndi),  &
    smatfic(ndi,ndi),sisomatfic(ndi,ndi), cmisomatfic(ndi,ndi,ndi,ndi),  &
    cisomatfic(ndi,ndi,ndi,ndi)
!     FILAMENTS NETWORK CONTRIBUTION
DOUBLE PRECISION :: filprops(10), affprops(7), chemprops(10), eqprops(11)
DOUBLE PRECISION :: cactin,cabp,ll,lambda0,mu0str,beta,nn,b0,bb
DOUBLE PRECISION :: phinet,r0,r0c,r0f,a,p,etac,na,mactin,rhoactin
DOUBLE PRECISION :: pknetfic(ndi,ndi),cmnetfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: snetfic(ndi,ndi),cnetfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: pknetficaf(ndi,ndi),pknetficnaf(ndi,ndi)
DOUBLE PRECISION :: snetficaf(ndi,ndi),snetficnaf(ndi,ndi)
DOUBLE PRECISION :: cmnetficaf(ndi,ndi,ndi,ndi), cmnetficnaf(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: cnetficaf(ndi,ndi,ndi,ndi), cnetficnaf(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: efi, kb, dx, dxc, Lp, theta
DOUBLE PRECISION :: R, Rfmax, Rbmax, Keq, Koff0, Kon0, Kcatch0
DOUBLE PRECISION :: cb(ndir), cb0, cbmax, thetab, thetaf0
DOUBLE PRECISION :: lamn(ndir), sfr, xi_t
DOUBLE PRECISION :: cb_tot, cb_tot_new, cf
DOUBLE PRECISION :: cb_upper, machep, tol
DOUBLE PRECISION :: Jc, f, df, dHdcb
DOUBLE PRECISION :: DTHETAFDJ, dcbdthetaf, dsficdthetaf(ndi,ndi), dsigisodthetaf(ndi,ndi)
!     FULLY IMPLICIT FREE FRACTION F(theta_f) = H(theta_f; mu, J, cb_tot(theta_f)), RTSAFE
INTEGER :: it_fi
INTEGER, PARAMETER :: MAXIT_FI = 50
DOUBLE PRECISION, PARAMETER :: TOL_FI_F = 1.0d-10, TOL_FI_TH = 1.0d-12
DOUBLE PRECISION :: cb_try(ndir), lamn_try(ndir), theta0, theta1, thx, thl, thh, Ffi, Fpfi, dthfi, dthold
DOUBLE PRECISION :: Heff, Hth0
!
!     JAUMMAN RATE CONTRIBUTION (REQUIRED FOR ABAQUS UMAT)
DOUBLE PRECISION :: cjr(ndi,ndi,ndi,ndi)
!     CAUCHY STRESS AND ELASTICITY TENSOR
DOUBLE PRECISION :: sigma(ndi,ndi),ddsigdde(ndi,ndi,ndi,ndi),  &
    ddpkdde(ndi,ndi,ndi,ndi)
!     OTHER TANGENT AUX TENSORS
DOUBLE PRECISION :: dpk2dthetaf(ndi,ndi), dthetafdc(ndi,ndi), dPK2ficdcb(ndi,ndi), &
                    dcbdc(ndi,ndi), cvolchem(ndi,ndi,ndi,ndi), &
                    dsigmadcb(ndi,ndi)
DOUBLE PRECISION :: stest(ndi,ndi), ctest(ndi,ndi,ndi,ndi)

INTEGER :: I1, J1, K1, L1


!----------------------------------------------------------------------
!-------------------------- INITIALIZATIONS ---------------------------
!----------------------------------------------------------------------
!     IDENTITY AND PROJECTION TENSORS
unit2=zero
unit4=zero
unit4s=zero
proje=zero
projl=zero
!     KINEMATICS
distgr=zero
c=zero
b=zero
cbar=zero
bbar=zero
ubar=zero
vbar=zero
rot=zero
det=zero
cbari1=zero
cbari2=zero
!     VOLUMETRIC
pkvol=zero
svol=zero
cvol=zero
k=zero
pv=zero
ppv=zero
ppv_chem=zero
ssev=zero
!     ISOCHORIC
siso=zero
pkiso=zero
pk2=zero
ciso=zero
cfic=zero
sfic=zero
pkfic=zero
!     ISOTROPIC
c10=zero
c01=zero
sseiso=zero
diso=zero
pkmatfic=zero
smatfic=zero
sisomatfic=zero
cmisomatfic=zero
cisomatfic=zero
!     FILAMENTS NETWORK
snetfic=zero
cnetfic=zero
pknetfic=zero
pknetficaf=zero
pknetficnaf=zero
snetficaf=zero
snetficnaf=zero
cmnetfic=zero
cmnetficaf=zero
cmnetficnaf=zero
cnetficaf=zero
cnetficnaf=zero
!     JAUMANN RATE
cjr=zero
!     TOTAL CAUCHY STRESS AND ELASTICITY TENSORS
sigma=zero
ddsigdde=zero
cvolchem=zero
!     FLUID FLUX
jfluid=zero
!----------------------------------------------------------------------
!------------------------ IDENTITY TENSORS ----------------------------
!----------------------------------------------------------------------
CALL onem(unit2,unit4,unit4s,ndi)
!----------------------------------------------------------------------
!------------------------ RANDOM GENERATION ---------------------------
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!------------------- MATERIAL CONSTANTS AND DATA ----------------------
!----------------------------------------------------------------------
!     VOLUMETRIC
k        = props(1)
!     ISOCHORIC ISOTROPIC
c10      = props(2)
c01      = props(3)
phinet   = props(4)
!     ACTIN/CROSSLINKERS
a        = props(5)  ! Ratio between contour length and end-to-end distance
r0c      = props(6)
etac     = props(7)
mu0str   = props(8)
beta     = props(9)
Lp       = props(10) ! Persistence length
theta    = props(11) ! Absolute temperature
dx       = props(12) ! CL reactive distance / bond length
!     AFFINE NETWORK
bb       = props(13)
lambda0  = props(14)
cactin   = props(15)
R        = props(16) ! CL to actin ratio
Rfmax    = props(17) ! Maximum free CL to actin ratio
Rbmax   = props(18) ! Maximum bound CL to actin ratio
!     SOLVENT
CHI    = PROPS(19)
D      = PROPS(20)
MU0    = PROPS(21)
VMOL   = PROPS(22)
Koff0  = PROPS(23)
Keq    = PROPS(24)
Kcatch0 = PROPS(25)
dxc     = PROPS(26)
!     STRETCH-FREE REBINDING (as 6_CL_DIFF_CS if PROPS(27) = 0)
IF (NPROPS >= 28) THEN
  sfr  = PROPS(27) ! 1 = on (bonds rebind at the current stretch, reference stretch lamn), 0 = off (lamn = 1)
  xi_t = PROPS(28) ! fraction of the network with stretch-free rebinding (0..1; the rest is permanent)
ELSE
  sfr  = zero
  xi_t = zero
END IF

!Other parameters (Check which of these will be actually needed in the UMAT and not only in the AFFCL subroutine)
! (double-precision literals: single-precision ones carry ~1e-8 relative errors)
kb = 1.380649d-5
b0 = Lp * theta * kb
rgas = 8.314462618d0
Mactin = 42.0d-3       ! [MDa]
rhoactin = 16.0d0      ! [MDa/microm]
NA = 6.022d5           ! [1/amol]
Kon0 = (Koff0 + Kcatch0) * Keq
! Maximum allowable CL concentration
cfmax = Rfmax * cactin
cbmax = Rbmax * cactin

filprops = (/a, r0c, etac, mu0str, beta, Lp, theta, dx, kb, NA/)
affprops = (/bb, lambda0, cactin, Mactin, rhoactin, sfr, xi_t/)
chemprops = (/cbmax, cfmax, CHI, D, MU0, VMOL, Koff0, Keq, Kcatch0, dxc/)
eqprops = (/a, r0c, etac, mu0str, beta, b0, lambda0, Koff0, Kcatch0, dx/(kb*theta), dxc/(kb*theta)/)


!        STATE VARIABLES AND CHEMICAL PARAMETERS
! IF ((kinc <= 1).AND.(kstep == 1)) THEN
IF (STATEV(1) == 0.0d0) THEN
!     CL CONCENTRATION
  cabp = cactin*R 
! Initial bound and free CL concentrations
  cb_upper = MIN(cabp, cbmax)
  machep = 2.22d-16
  tol = 1.0d-12
  ! write(*,*) 'Calling pullchem at t=0'
  CALL pullchem(cb0, zero, cb_upper, machep, tol, cabp, cfmax, cbmax, CHI, Keq, eqprops)
  ! write(*,*) 'cb0 = ', cb0
  thetaf0 = (cabp - cb0) / cfmax
  ! write(*,*) 'thetaf0 = ', thetaf0
  CALL initialize(statev,thetaf0,vmol,cb0,cfmax)
END IF
!        READ STATEV
CALL sdvread(statev, thetaf_t, cb, cb_tot, lamn)
!----------------------------------------------------------------------
!---------------------------- KINEMATICS ------------------------------
!----------------------------------------------------------------------
!     DISTORTION GRADIENT
CALL fslip(dfgrd1,distgr,det,ndi)
!     INVERSE OF DEFORMATION GRADIENT
CALL matinv3d(dfgrd1,dfgrd1inv,ndi)
!     INVERSE OF DISTORTION GRADIENT
CALL matinv3d(distgr,distgrinv,ndi)
!     CAUCHY-GREEN DEFORMATION TENSORS
CALL deformation(dfgrd1,c,b,ndi)
CALL deformation(distgr,cbar,bbar,ndi)
!     INVERSE OF CAUCHY-GREEN DEFORMATION TENSORS
CALL matinv3d(c,cinv,ndi)
CALL matinv3d(cbar,cbarinv,ndi)
!     INVARIANTS OF DEVIATORIC DEFORMATION TENSORS
CALL invariants(cbar,cbari1,cbari2,ndi)
!     STRETCH TENSORS
CALL stretch(cbar,bbar,ubar,vbar,ndi)
!     ROTATION TENSORS
CALL rotation(distgr,rot,ubar,ndi)
!     DEVIATORIC PROJECTION TENSORS
CALL projeul(unit2,unit4s,proje,ndi)

CALL projlag(c,unit4,projl,ndi)
!----------------------------------------------------------------------
!---------------------- COUPLED DIFFUSION -----------------------------
!----------------------------------------------------------------------

!     1. Solve for current free crosslinker fraction (THETAF_TAU)
      ARGS(1) = MU_TAU
      ARGS(2) = MU0
      ARGS(3) = RGAS
      ARGS(4) = THETA
      ARGS(5) = CHI
      ARGS(6) = VMOL
      ARGS(7) = K
      ARGS(8) = DET
      ARGS(9) = CB_TOT
      ARGS(10) = CFMAX
      ARGS(11) = cactin * R      ! initial total crosslinker content c0 (reference state of J^c)
      ! write(*,*) 'ARGS = ', ARGS
      ! write(*,*) 'THETAF_T =', THETAF_T
      CALL SOLVETHETAF(THETAF_TAU, ARGS, NARGS, THETAF_T, PNEWDT)

      CALL erfi(efi,bb)
      IF ((phinet > zero) .AND. (PNEWDT >= one)) THEN
        theta0 = THETAF_TAU
        CALL fi_eval(theta0, Ffi, Fpfi)
        IF ((PNEWDT >= one) .AND. (ABS(Ffi) >= TOL_FI_F)) THEN
          CALL SOLVETHETAF(theta1, ARGS, NARGS, theta0, PNEWDT)   ! phi(theta0)
          Jc = one + VMOL * (ARGS(9) + cfmax * theta0 - cactin * R)
          dHdcb = k * vmol**two / (RGAS * THETA) * (one + DLOG(det/Jc)) / Jc**two
          IF (ABS(dHdcb * dcbdthetaf) * ABS(theta1 - theta0) < TOL_FI_F) THEN
            THETAF_TAU = theta1
          ELSE
            IF (dcbdthetaf <= zero) THEN                     ! SFR: [theta0, phi(theta0)] not guaranteed
              thl = 1.0d-8                                   ! F < 0
              thh = one - 1.0d-8
            ELSE IF (Ffi < zero) THEN
              thl = theta0
              thh = theta1
            ELSE
              thl = theta1
              thh = theta0
            END IF
            thx = theta0
            dthold = ABS(thh - thl)
            dthfi = dthold
            DO it_fi = 1, MAXIT_FI
              ! Bisect if Newton is out of the bracket or not decreasing fast enough
              IF ((((thx - thh) * Fpfi - Ffi) * ((thx - thl) * Fpfi - Ffi) >= zero) .OR. &
                  (ABS(two * Ffi) > ABS(dthold * Fpfi))) THEN
                dthold = dthfi
                dthfi = half * (thh - thl)
                thx = thl + dthfi
              ELSE
                dthold = dthfi
                dthfi = Ffi / Fpfi
                thx = thx - dthfi
              END IF
              IF (ABS(dthfi) < TOL_FI_TH) EXIT
              CALL fi_eval(thx, Ffi, Fpfi)
              IF ((PNEWDT < one) .OR. (ABS(Ffi) < TOL_FI_F)) EXIT
              IF (Ffi < zero) THEN
                thl = thx
              ELSE
                thh = thx
              END IF
            END DO
            IF (it_fi > MAXIT_FI) THEN
              PNEWDT = MIN(PNEWDT, half)
              write(*,*) 'WARNING: fully implicit thetaf not converged, noel =', noel, ' |F| =', ABS(Ffi), &
                         ' -> PNEWDT =', PNEWDT
            END IF
            THETAF_TAU = thx
          END IF
        END IF
      END IF

      ! Final network evaluation at the accepted theta_f: bound contents cb^tau
      dcbdc = zero
      dcbdthetaf = zero
      dsficdthetaf = zero
      cb_tot_new = cb_tot
      IF (phinet > zero) THEN
        CALL affclnetfic_discrete(snetficaf,cnetficaf,distgr,unit2,filprops,  &
             affprops,efi,noel,det,prefdir,ndi,cb,dtime,chemprops, &
             thetaf_tau, cb_tot_new, dPK2ficdcb, dcbdc, pnewdt, dcbdthetaf, dsficdthetaf, lamn, 2)
      END IF
      ARGS(9) = cb_tot_new

      cf = THETAF_TAU * cfmax
      ! Swelling relative to the initial content c0 = cactin*R (stress-free reference: Jc = 1 at t = 0);
      ! dJc/dc = VMOL as before, so the derivatives below are unchanged (same definition in thetafFunc)
      Jc = 1.0d0 + VMOL * (cb_tot_new + cf - cactin * R)

      ! Evaluate tangent at converged root
      CALL thetafFunc(THETAF_TAU, f, df, ARGS, NARGS)
      dHdcb = k * vmol**two / (RGAS * THETA) * (one + DLOG(det/Jc)) / Jc**two
      Heff = df + dHdcb * dcbdthetaf
      Hth0 = df - dHdcb * cfmax
      DTHETAFDMU = one / (RGAS * THETA * Heff)

      ! write(*,*) 'kinc = ', kinc
      ! write(*,*) 'kstep = ', kstep
      ! write(*,*) 'DTIME = ', DTIME

      IF (DTIME > 1.0d-12) THEN
        DTHETAFDT = (THETAF_TAU - THETAF_T) / DTIME            
      ELSE
        DTHETAFDT = 0.0d0
      END IF


      ! Fluid mobility (spatial, consistent with the spatial flux j = -m grad(mu)):
      !   m = D/(RT) * (cf/J) * (1 - thetaf), with cf/J the free CL concentration per unit
      !   deformed volume (cf is per unit reference volume; J = det = Jbar with F-bar).
      !   D is the Fickian diffusion coefficient; the 1/(RT) (Einstein relation) gives
      !   D_eff = m * dmu/d(cf/J) = D * [1 - 2*chi*thetaf*(1-thetaf)] ~ D
      !   (Chester, Di Leo & Anand 2015, Eqs. 2.11 and 4.1)
      MFLUID = D / (RGAS * THETA) * cf * (1.0d0 - THETAF_TAU) / det

      ! Mobility tangents
      DMDMU = D / (RGAS * THETA) * cfmax * (1.0d0 - 2.0d0 * THETAF_TAU) * DTHETAFDMU / det
      ! dm/dJ: through thetaf (Implicit Function Theorem on H(THETAF_TAU, mu, J) = 0,
      !   dthetaf/dJ = (k*Vmol) / (RT * Jc * det * df)) and through the 1/J of the concentration:
      !   dm/dJ = dm/dthetaf * dthetaf/dJ - m/J
      DTHETAFDJ = (k * VMOL) / (RGAS * THETA * Jc * det * Heff)
      DMDJ  = D / (RGAS * THETA) * cfmax * (1.0d0 - 2.0d0 * THETAF_TAU) * DTHETAFDJ / det &
            - MFLUID / det

      ! Fluid flux vector (for visualization/SVARS)
      jfluid = -MFLUID * DMUDX

      
!----------------------------------------------------------------------
!--------------------- CONSTITUTIVE RELATIONS  ------------------------
!----------------------------------------------------------------------
!---- VOLUMETRIC ------------------------------------------------------
!     STRAIN-ENERGY
CALL vol(ssev,pv,ppv,k,det,Jc)
! Add chemical contribution (Kuu Term 2)
ppv_chem = - (k**two * VMOL**two * (CFMAX + dcbdthetaf)) / (det * Jc**two) * DTHETAFDMU
ppv = ppv+ ppv_chem
! write(*,*) 'DTHETAFDMU = ', DTHETAFDMU
! write(*,*) 'ppv = ', ppv
! write(*,*) 'ppv_chem = ', ppv_chem
!---- ISOCHORIC ISOTROPIC ---------------------------------------------
IF (phinet < one) THEN
!     STRAIN-ENERGY
  CALL isomat(sseiso,diso,c10,c01,cbari1,cbari2)
!     PK2 'FICTICIOUS' STRESS TENSOR
  CALL pk2isomatfic(pkmatfic,diso,cbar,cbari1,unit2,ndi)
!     CAUCHY 'FICTICIOUS' STRESS TENSOR
  CALL sigisomatfic(sisomatfic,pkmatfic,distgr,det,ndi)
!     'FICTICIOUS' MATERIAL ELASTICITY TENSOR
  CALL cmatisomatfic(cmisomatfic,cbar,cbari1,cbari2, diso,unit2,unit4,det,ndi)
!     'FICTICIOUS' SPATIAL ELASTICITY TENSOR
  CALL csisomatfic(cisomatfic,cmisomatfic,distgr,det,ndi)
  
END IF
!---- FILAMENTS NETWORK -----------------------------------------------
!     Computed with the fully implicit free fraction: snetficaf, cnetficaf, cb^tau,
!     cb_tot_new, dcbdthetaf (G), dcbdc, dsficdthetaf

! Macroscopic reaction source (homogenized binding rate) and the source derivatives (Part C, item 1)
IF (DTIME > 1.0d-12) THEN
  RMACRO = (cb_tot_new - cb_tot) / DTIME
  DSRCDMU = (CFMAX + dcbdthetaf) * DTHETAFDMU / DTIME
  DSRCDJ  = (CFMAX + dcbdthetaf) * DTHETAFDJ / DTIME
  DSRCDGU = - dcbdc / DTIME * Hth0 / Heff
ELSE
  RMACRO = 0.0d0
  DSRCDMU = 0.0d0
  DSRCDJ  = 0.0d0
  DSRCDGU = 0.0d0
END IF
! write(*,*) 'cb_tot = ', cb_tot
! write(*,*) 'cb_tot_new = ', cb_tot_new

!      PKNETFIC=PKNETFICNAF+PKNETFICAF
snetfic=snetficnaf+snetficaf
!      CMNETFIC=CMNETFICNAF+CMNETFICAF
cnetfic=cnetficnaf+cnetficaf
!----------------------------------------------------------------------
!     STRAIN-ENERGY
SSE=SSEV+SSEISO
!     PK2 'FICTICIOUS' STRESS
pkfic=(one-phinet)*pkmatfic+pknetfic
!     CAUCHY 'FICTICIOUS' STRESS
sfic=(one-phinet)*sisomatfic+snetfic
!     MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
cmfic=(one-phinet)*cmisomatfic+cmnetfic
!     SPATIAL 'FICTICIOUS' ELASTICITY TENSOR
cfic=(one-phinet)*cisomatfic+cnetfic
!----------------------------------------------------------------------
!-------------------------- STRESS MEASURES ---------------------------
!----------------------------------------------------------------------
!---- VOLUMETRIC ------------------------------------------------------
!      PK2 STRESS
! CALL pk2vol(pkvol,pv,c,ndi)
CALL pk2vol(pkvol,pv,c,ndi,det)
!      CAUCHY STRESS
CALL sigvol(svol,pv,unit2,ndi)
!---- ISOCHORIC -------------------------------------------------------
!      PK2 STRESS
CALL pk2iso(pkiso,pkfic,projl,det,ndi)
!      CAUCHY STRESS
CALL sigiso(siso,sfic,proje,ndi)
!      ACTIVE CAUCHY STRESS
!      CALL SIGISO(SACTISO,SNETFICAF,PROJE,NDI)

!      CALL SPECTRAL(SACTISO,SACTVL,SACTVC)
!---- VOLUMETRIC + ISOCHORIC ------------------------------------------
!      PK2 STRESS
pk2 = pkvol + pkiso
!      CAUCHY STRESS
sigma = svol + siso

!----------------------------------------------------------------------
!-------------------- MATERIAL ELASTICITY TENSOR ----------------------
!----------------------------------------------------------------------

!---- VOLUMETRIC ------------------------------------------------------

!      CALL METVOL(CMVOL,C,PV,PPV,DET,NDI)

!---- ISOCHORIC -------------------------------------------------------

!      CALL METISO(CMISO,CMFIC,PROJL,PKISO,PKFIC,C,UNIT2,DET,NDI)

!----------------------------------------------------------------------

!      DDPKDDE=CMVOL+CMISO

!----------------------------------------------------------------------
!--------------------- SPATIAL ELASTICITY TENSOR ----------------------
!----------------------------------------------------------------------

!---- VOLUMETRIC ------------------------------------------------------

CALL setvol(cvol,pv,ppv,unit2,unit4s,ndi)

!---- ISOCHORIC -------------------------------------------------------

CALL setiso(ciso,cfic,proje,siso,sfic,unit2,ndi)

!-----JAUMMAN RATE ----------------------------------------------------

CALL setjr(cjr,sigma,unit2,ndi)

!----------------------------------------------------------------------

!     ELASTICITY TENSOR
ddsigdde=cvol+ciso ! +cvolchem ! +cjr
! if (npt==1) then
!   write(*,*) 'cvolchem / ddsigdde = ', cvolchem / ddsigdde
! end if
!----------------------------------------------------------------------
!------------------------- CROSS-COUPLINGS ----------------------------
!----------------------------------------------------------------------
!     DISPLACEMENT - CHEMICAL POTENTIAL MODULUS
! DO I1 = 1, NDI
!     DO J1 = 1, NDI
!       ! Derivative of Cauchy stress with respect to phi
!       SPUCMOD(I1,J1) = (K / (DETFE * PHI_TAU)) * UNIT2(I1,J1) * DPHIDMU
!     END DO
! END DO

!     CHEMICAL POTENTIAL - DISPLACEMENT MODULUS: assembled in the element from MFLUID, DMDJ,
!     DSRCDJ and DSRCDGU (u3d8.f90, Part C items 1 and 5)


!     CAUCHY STRESS - CHEMICAL POTENTIAL MODULUS (dS / dMu)
!     volumetric (free CL swelling) + isochoric network through cb_i(thetaf) (Part C, item 4):
!     dsigma/dthetaf = -K*Vmol*cfmax/(J*Jc) I + P : sum_i dsfic_i/dcb_i * S^theta_i
  CALL contraction42(dsigisodthetaf, proje, dsficdthetaf, ndi)
  DO I1 = 1, NDI
      DO J1 = 1, NDI
        DSIGDMU(I1,J1) = ( -((K * VMOL * (CFMAX + dcbdthetaf)) / (DET * Jc)) * UNIT2(I1,J1) &
                         + dsigisodthetaf(I1,J1) ) * DTHETAFDMU
      END DO
  END DO

!     NETWORK STRESS THROUGH THETAF(J) (Kuu term 4, Part C item 7): push-forward of
!     2 dSiso/dthetaf (x) dthetaf/dC, with dthetaf/dC = J/2 dthetaf/dJ C^-1:
!     c_ijkl += J * dthetaf/dJ * dsigiso/dthetaf_ij * delta_kl
!     (the volumetric counterpart is zero: Jc uses the lagged cb_tot, as for cvolchem)
  DO I1 = 1, NDI
    DO J1 = 1, NDI
      DO K1 = 1, NDI
        ddsigdde(I1,J1,K1,K1) = ddsigdde(I1,J1,K1,K1) + DET * DTHETAFDJ * dsigisodthetaf(I1,J1)
      END DO
    END DO
  END DO

!     FULLY IMPLICIT: stretch dependence through cb_tot^tau (dc_mech = -dcbdc:eps at fixed J and mu):
!     Term 5  - pressure through J^c:   dc = Hth0/Heff * dc_mech  ->  +K*Vmol/(J*Jc) * Hth0/Heff * I (x) dcbdc
!     Term 4b - network through theta_f: dthetaf = -dHdcb/Heff * dc_mech -> dsigiso/dthetaf (x) dHdcb/Heff * dcbdc
  DO I1 = 1, NDI
    DO J1 = 1, NDI
      DO K1 = 1, NDI
        DO L1 = 1, NDI
          ddsigdde(I1,J1,K1,L1) = ddsigdde(I1,J1,K1,L1) &
              + ((K * VMOL) / (DET * Jc) * Hth0 / Heff * UNIT2(I1,J1) + dHdcb / Heff * dsigisodthetaf(I1,J1)) &
              * dcbdc(K1,L1)
        END DO
      END DO
    END DO
  END DO


!----------------------------------------------------------------------
!------------------------- INDEX ALLOCATION ---------------------------
!----------------------------------------------------------------------
!     VOIGT NOTATION  - FULLY SIMMETRY IMPOSED
CALL indexx(stress,ddsdde,sigma,ddsigdde,ntens,ndi)

!----------------------------------------------------------------------
!--------------------------- STATE VARIABLES --------------------------
!----------------------------------------------------------------------
!     DO K1 = 1, NTENS
!      STATEV(1:27) = VISCOUS TENSORS
CALL sdvwrite(det,statev,stress,thetaf_tau,dmudx,Vmol,jfluid,cb,cb_tot_new,cfmax,lamn)
! CALL sdvwrite(det,etac_sdv,statev)
!     END DO
!----------------------------------------------------------------------
RETURN

!----------------------------------------------------------------------
!----------------------- AUXILIAR SUBROUTINES -------------------------
!----------------------------------------------------------------------

CONTAINS

  SUBROUTINE fi_eval(thv, Fv, Fpv)
    ! Kinetics of all directions at theta_f = thv from the state at t (on a copy of cb, kinetics only),
    ! then F = H(thv; mu, J, cb_tot(thv)) (thetafFunc) and F' = dH/dthetaf + dH/dcb_tot * G.
    ! Sets ARGS(9) = cb_tot(thv).
    DOUBLE PRECISION, INTENT(IN)  :: thv
    DOUBLE PRECISION, INTENT(OUT) :: Fv, Fpv
    DOUBLE PRECISION :: thc, Jcv
    cb_try = cb
    lamn_try = lamn
    dcbdthetaf = zero
    CALL affclnetfic_discrete(snetficaf,cnetficaf,distgr,unit2,filprops,  &
        affprops,efi,noel,det,prefdir,ndi,cb_try,dtime,chemprops, &
        thv, cb_tot_new, dPK2ficdcb, dcbdc, pnewdt, dcbdthetaf, dsficdthetaf, lamn_try, 1)
    ARGS(9) = cb_tot_new
    thc = thv
    CALL thetafFunc(thc, Fv, Fpv, ARGS, NARGS)
    Jcv = one + VMOL * (cb_tot_new + cfmax * thv - cactin * R)
    Fpv = Fpv + k * vmol**two / (RGAS * THETA) * (one + DLOG(det/Jcv)) / Jcv**two * dcbdthetaf
  END SUBROUTINE fi_eval

END SUBROUTINE material
!----------------------------------------------------------------------
!--------------------------- END OF UMAT ------------------------------
!----------------------------------------------------------------------

