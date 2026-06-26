!>********************************************************************
!> Record of revisions:                                              |
!>        Date        Programmer        Description of change        |
!>        ====        ==========        =====================        |
!>                                                                   |
!>--------------------------------------------------------------------
!>     Description:
!C>     UMAT: USER MATERIAL FOR THE FULL NETWORK MODEL.
!C>                 AFFINE DEFORMATIONS
!C>     UEXTERNALDB: READ FILAMENTS ORIENTATION AND PREFERED DIRECTION
!>--------------------------------------------------------------------
!>---------------------------------------------------------------------

! SUBROUTINE material(stress,statev,ddsdde,sse,spd,scd, rpl,ddsddt,drplde,drpldt,  &
!     stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,  &
!     ndi,nshr,ntens,nstatev,props,nprops,coords,drot,pnewdt,  &
!     celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

    SUBROUTINE MATERIAL(SIGMA,STATEV,DDSIGDDE,DFGRD0,DFGRD1,DET, &
    TIME,DTIME,PREDEF,NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS, &
    PNEWDT,NOEL,NPT,KSTEP,KINC,MU_TAU,THETAF_T,THETAF_TAU,DTHETAFDT, &
      DTHETAFDMU,RMACRO,MFLUID,DMDMU,DMUDX,DMDJ,VMOL,CFMAX,DSIGDMU,SPCUMODFAC)
!
use global  
IMPLICIT NONE
!----------------------------------------------------------------------
!--------------------------- DECLARATIONS -----------------------------
!----------------------------------------------------------------------
INTEGER :: NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT, &
            LAYER, KSPT, KSTEP, KINC

INTEGER, PARAMETER :: nargs = 10

REAL(KIND=8) :: STRESS(NTENS), STATEV(NSTATEV), &
                DDSDDE(NTENS,NTENS), DDSDDT(NTENS), DRPLDE(NTENS), &
                STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1), DPRED(1), &
                PROPS(NPROPS), COORDS(3), DROT(3,3), DFGRD0(3,3), DFGRD1(3,3), &
                FIBORI(NELEM,4), ARGS(NARGS)

REAL(8), INTENT(IN)      :: MU_TAU, THETAF_T, DMUDX(3,1)
! REAL(8), INTENT(OUT)     :: SPUCMOD(NDI,NDI), SPCUMODFAC(NDI,NDI)
REAL(8), INTENT(OUT)     :: DSIGDMU(NDI,NDI), SPCUMODFAC(NDI,NDI)
REAL(8), INTENT(OUT)     :: THETAF_TAU, DTHETAFDT, DTHETAFDMU, RMACRO! DPHIDMU, DPHIDOTDMU
REAL(8), INTENT(OUT)     :: MFLUID, DMDMU, DMDJ, VMOL, CFMAX

! cfmax can probably be defined at the element level

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
    ubar(ndi,ndi),vbar(ndi,ndi),rot(ndi,ndi), dfgrd1inv(ndi,ndi)
DOUBLE PRECISION :: det,detfe, detfs,cbari1,cbari2
!     VOLUMETRIC CONTRIBUTION
DOUBLE PRECISION :: pkvol(ndi,ndi),svol(ndi,ndi),  &
    cvol(ndi,ndi,ndi,ndi),cmvol(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: k,pv,ppv,ssev
!     ISOCHORIC CONTRIBUTION
DOUBLE PRECISION :: siso(ndi,ndi),pkiso(ndi,ndi),pk2(ndi,ndi),  &
    ciso(ndi,ndi,ndi,ndi),cmiso(ndi,ndi,ndi,ndi),  &
    sfic(ndi,ndi),cfic(ndi,ndi,ndi,ndi), pkfic(ndi,ndi),cmfic(ndi,ndi,ndi,ndi)
!     ISOCHORIC ISOTROPIC CONTRIBUTION
DOUBLE PRECISION :: c10,c01,sseiso,diso(5),pkmatfic(ndi,ndi),  &
    smatfic(ndi,ndi),sisomatfic(ndi,ndi), cmisomatfic(ndi,ndi,ndi,ndi),  &
    cisomatfic(ndi,ndi,ndi,ndi)
!     FILAMENTS NETWORK CONTRIBUTION
DOUBLE PRECISION :: filprops(10), affprops(5) ! affprops(6)
DOUBLE PRECISION :: cactin,cabp,ll,lambda0,mu0str,beta,nn,b0,bb
DOUBLE PRECISION :: phinet,r0,r0c,r0f,a,p,etac,na,mactin,rhoactin
DOUBLE PRECISION :: pknetfic(ndi,ndi),cmnetfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: snetfic(ndi,ndi),cnetfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: pknetficaf(ndi,ndi),pknetficnaf(ndi,ndi)
DOUBLE PRECISION :: snetficaf(ndi,ndi),snetficnaf(ndi,ndi)
DOUBLE PRECISION :: cmnetficaf(ndi,ndi,ndi,ndi), cmnetficnaf(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: cnetficaf(ndi,ndi,ndi,ndi), cnetficnaf(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: efi, kb, dx, Lp, theta
DOUBLE PRECISION :: R, Rfmax, Rbmax, Keq, Koff0, Kon0
DOUBLE PRECISION :: cb(ndir), cb0, cbmax, thetab, thetaf !, cfmax
DOUBLE PRECISION :: cb_tot, cb_tot_new, cf
DOUBLE PRECISION :: cb_upper, machep, tol
DOUBLE PRECISION :: Jc, f, df

! INTEGER :: nterm,factor 
!
!     JAUMMAN RATE CONTRIBUTION (REQUIRED FOR ABAQUS UMAT)
DOUBLE PRECISION :: cjr(ndi,ndi,ndi,ndi)
!     CAUCHY STRESS AND ELASTICITY TENSOR
DOUBLE PRECISION :: sigma(ndi,ndi),ddsigdde(ndi,ndi,ndi,ndi),  &
    ddpkdde(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: stest(ndi,ndi), ctest(ndi,ndi,ndi,ndi)

! DECLARATIONS FOR RANDOM GENERATION
INTEGER (kind=4) :: seed1, seed2
INTEGER (kind=4) :: test, test_num
INTEGER (kind=4) :: l, i, idx
CHARACTER(len=100) :: phrase
!REAL(kind=4) , allocatable :: etac_array(:), array(:)
DOUBLE PRECISION :: etac_sdv(nsdv-1)
!REAL(kind=4) :: l_bound, h_bound
REAL(kind=4) :: mean, sd

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

!Other parameters (Check which of these will be actually needed in the UMAT and not only in the AFFCL subroutine)
kb = 1.380649e-5      
b0 = Lp * theta * kb
rgas = 8.314462618
Mactin = 42.0e-3       ! [MDa]
rhoactin = 16.0        ! [MDa/microm]
NA = 6.022e5           ! [1/amol]
Kon0 = Koff0 * Keq

filprops = (/a, r0c, etac, mu0str, beta, Lp, theta, dx, kb, NA/)
affprops = (/bb, lambda0, cactin, Mactin, rhoactin/)
! affprops = (/bb, lambda0, cactin, R, Rfmax, Rbmax, kb, b0, rgas, Mactin, rhoactin, NA/)

! All of these will be needed (but not here)
! Check whether they should be at UEL/UMAT/AFFCL DIRECTION
!     CL CONCENTRATION
! cabp = cactin*R
! write(*,*) 'cabp = ', cabp
!     FILAMENT END-TO-END DISTANCE
! r0f = 1.6 * cabp**(-2.0/5.0) ! AFFCL DIRECTION
! write(*,*) 'r0f = ', r0f
!     FILAMENT CONTOUR LENGTH
! ll = a * r0f ! AFFCL DIRECTION
! write(*,*) 'll = ', ll
!     FILAMENT DENSITY
! na = 6.022e23
! mactin = 42.0          ! [kDa]
! rhoactin = 16.0        ! [MDa/microm]
! nn = cactin/ll * na * mactin / rhoactin * 1.0e-24 ! AFFCL DIRECTION
! write(*,*) 'nn = ', nn

!     CL CONCENTRATION
!!! THIS NEEDS TO BE CHANGED AFTER DIFFUSION IS IMPLEMENTED IN UEL
cabp = cactin*R  ! <-- Placeholder: Replace with true UEL cR later!
! Maximum allowable CL concentration
cfmax = Rfmax * cactin
cbmax = Rbmax * cactin

!        STATE VARIABLES AND CHEMICAL PARAMETERS
IF ((time(1) == zero).AND.(kstep == 1)) THEN
  ! Initial bound and free CL concentrations
  cb_upper = MIN(cabp, cbmax)
  machep = 2.22d-16
  tol = 1.0d-8
  CALL pullchem(cb0, zero, cb_upper, machep, tol, cabp, cfmax, cbmax, CHI, Keq)
  CALL initialize(statev,thetaf_t,vmol,cb0)
END IF
!        READ STATEV
CALL sdvread(statev, cb, cb_tot)
! --------------------------------------------
cf = cabp - cb_tot
thetaf = cf / cfmax
! Avoid numerical issues
thetaf = MIN(MAX(thetaf, 1.0d-6), 1.0d0 - 1.0d-6)
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


!     1. Solve for current free crosslinker fraction (THETAF)
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
      CALL SOLVETHETAF(THETAF_TAU, ARGS, NARGS, THETAF_T)

      thetaf = THETAF_TAU
      cf = thetaf * cfmax

      ! Evaluate tangent at converged root
      CALL thetafFunc(thetaf, f, df, ARGS, NARGS)
      DTHETAFDMU = one / df

      ! Rate of free crosslinker fraction
      DTHETAFDT = (THETAF_TAU - THETAF_T) / DTIME

      ! Fluid mobility and permeability
      MFLUID = D * cf * (1.0d0 - thetaf)

      ! Mobility tangents
      DMDMU = D * cfmax * (1.0d0 - 2.0d0 * thetaf) * DTHETAFDMU
      DMDJ  = 0.0d0   ! Mobility no longer depends on volume!

      ! Fluid flux vector (just visualization/SVARS)
      jfluid = -MFLUID * DMUDX
      
!       DETFE = DET * PHI_TAU
!       !     2. Time rate of swelling
!       DPDT = (PHI_TAU - PHI_T) / DTIME
      
!       !     3. Analytical derivatives of PHI
!       DPHIDMU = (ONE / (RGAS * THETA)) / &
!       ( (ONE / (PHI_TAU - ONE)) + ONE + TWO * CHI * PHI_TAU &
!       - ((VMOL * K) / (RGAS * THETA * PHI_TAU)) &
!       + ((VMOL * K) / (RGAS * THETA * PHI_TAU)) * DLOG(DETFE) )
      
!       DPHIDJ  = ( ((VMOL * K) / (RGAS * THETA * DET)) &
!       - ((VMOL * K) / (RGAS * THETA * DET)) * DLOG(DETFE) ) / &
!       ( (ONE / (PHI_TAU - ONE)) + ONE + TWO * CHI * PHI_TAU &
!       - ((VMOL * K) / (RGAS * THETA * PHI_TAU)) &
!       + ((VMOL * K) / (RGAS * THETA * PHI_TAU)) * DLOG(DETFE) )
      
!       !     4. Numerical Perturbation for D(PHIDOT)/DMU
!       IF (DABS(MU_TAU) > ONE) THEN
!         DELTAMU = DABS(MU_TAU) * 1.D-8
!       ELSE
!         DELTAMU = 1.D-8
!       END IF
      
!       ARGS(1) = MU_TAU + DELTAMU
!       CALL SOLVEPHI(PHI_PER, ARGS, NARGS, PHI_T)
!       DPDT_PER = (PHI_PER - PHI_T) / DTIME
      
!       ARGS(1) = MU_TAU - DELTAMU
!       CALL SOLVEPHI(PHI_M, ARGS, NARGS, PHI_T)
!       DPDT_M = (PHI_M - PHI_T) / DTIME
      
!       DPHIDOTDMU = (DPDT_PER - DPDT_M) / (TWO * DELTAMU)
      
!       !     5. Fluid mobility and permeability
!       MFLUID = (D * (ONE / PHI_TAU - ONE)) / (DET * VMOL * RGAS * THETA)
!       DMDPHI = -(D / (DET * VMOL * PHI_TAU * PHI_TAU * RGAS * THETA))
!       DMDMU  = DMDPHI * DPHIDMU
!       DMDJ   = DMDPHI * DPHIDJ

!       !    6. Fluid flux vector (just for plotting)
!       JFLUID = -MFLUID * DMUDX
      
!----------------------------------------------------------------------
!--------------------- CONSTITUTIVE RELATIONS  ------------------------
!----------------------------------------------------------------------
!---- VOLUMETRIC ------------------------------------------------------
!     STRAIN-ENERGY
! THIS NEEDS TO BE CHANGED!!!!!!!!!!!!!!
Jc = 1.0d0 + VMOL * cb_tot + VMOL * cfmax * thetaf
CALL vol(ssev,pv,ppv,k,det,Jc)

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
!     IMAGINARY ERROR FUNCTION BASED ON DISPERSION PARAMETER
! CALL erfi(efi,bb,nterm) ! (original)
CALL erfi(efi,bb)
!     'FICTICIOUS' PK2 STRESS AND MATERIAL ELASTICITY TENSORS
!------------ AFFINE NETWORK --------------
IF ((phinet > zero) .AND. (nn > zero)) THEN
  ! GET CL STIFFNESS DISTRIBUTION FOR CURRENT GP
  !CALL getprops_gp(noel, npt, etadir, etadir_array)
  CALL affclnetfic_discrete(snetficaf,cnetficaf,distgr,filprops,  &
      affprops,efi,noel,det,prefdir,ndi,cb,dtime,cabp,cfmax,cbmax,chi,Keq,Koff0, &
      thetaf, cb_tot_new)
END IF

! Macroscopic reaction source (homogenized binding rate)
RMACRO = (cb_tot_new - cb_tot) / DTIME

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
ddsigdde=cvol+ciso+cjr

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

!     CHEMICAL POTENTIAL - DISPLACEMENT MODULUS
DO I1 = 1, NDI
    DO J1 = 1, NDI
      SPCUMODFAC(I1,J1) = MFLUID * UNIT2(I1,J1)
    END DO
END DO


!     CAUCHY STRESS - CHEMICAL POTENTIAL MODULUS (dS / dMu)
  DO I1 = 1, NDI
      DO J1 = 1, NDI
        DSIGDMU(I1,J1) = -((K * VMOL * CFMAX) / (DET * Jc)) * UNIT2(I1,J1) * DTHETAFDMU
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
CALL sdvwrite(det,statev,stress,thetaf_tau,dmudx,Vmol,jfluid,cb,cb_tot_new)
! CALL sdvwrite(det,etac_sdv,statev)
!     END DO
!----------------------------------------------------------------------
RETURN
END SUBROUTINE material
!----------------------------------------------------------------------
!--------------------------- END OF UMAT ------------------------------
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!----------------------- AUXILIAR SUBROUTINES -------------------------
!----------------------------------------------------------------------
!                         INPUT FILES
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!                         KINEMATIC QUANTITIES
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!                         STRESS TENSORS
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!                   LINEARISED ELASTICITY TENSORS
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------- UTILITY SUBROUTINES --------------------------
!----------------------------------------------------------------------

