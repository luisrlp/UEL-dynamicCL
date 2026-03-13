!********************************************************************
! Record of revisions:                                              |
!        Date        Programmer        Description of change        |
!        ====        ==========        =====================        |
!     15/11/2017    Joao Ferreira      cont mech general eqs        |
!     01/11/2018    Joao Ferreira      comments added               |
!--------------------------------------------------------------------
!     Description:
!     UMAT: IMPLEMENTATION OF THE CONSTITUTIVE EQUATIONS BASED UPON 
!     A STRAIN-ENERGY FUNCTION (SEF).
!     THIS CODE, AS IS, EXPECTS A SEF BASED ON THE INVARIANTS OF THE 
!     CAUCHY-GREEN TENSORS. A VISCOELASTIC COMPONENT IS ALSO 
!     INCLUDED IF NEEDED. 
!     YOU CAN CHOOSE TO COMPUTE AT THE MATERIAL FRAME AND THEN 
!     PUSHFORWARD OR  COMPUTE AND THE SPATIAL FRAME DIRECTLY.
!--------------------------------------------------------------------
!     IF YOU WANT TO ADAPT THE CODE ACCORDING TO YOUR SEF:
!    ISOMAT - DERIVATIVES OF THE SEF IN ORDER TO THE INVARIANTS
!    ADD OTHER CONTRIBUTIONS: STRESS AND TANGENT MATRIX
!-------------------------------------------------------------------- 
!      STATE VARIABLES: CHECK ROUTINES - INITIALIZE, WRITESDV, READSDV.
!--------------------------------------------------------------------              
!     UEXTERNALDB: READ FILAMENTS ORIENTATION AND PREFERED DIRECTION
!--------------------------------------------------------------------
!---------------------------------------------------------------------
      SUBROUTINE MATERIAL(SIGMA,STATEV,DDSIGDDE,DFGRD0,DFGRD1,DET, &
      TIME,DTIME,PREDEF,NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS, &
      PNEWDT,NOEL,NPT,KSTEP,KINC)
!----------------------------------------------------------------------
!--------------------------- DECLARATIONS -----------------------------
!----------------------------------------------------------------------
      use global
      IMPLICIT NONE

!     ADD COMMON BLOCKS HERE IF NEEDED (and in uexternal)
!      COMMON /KBLOCK/KBLOCK

      CHARACTER(len=8) :: CMNAME

      INTEGER :: NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT, &
                 LAYER, KSPT, KSTEP, KINC

      REAL(KIND=8) :: STRESS(NTENS), STATEV(NSTATEV), &
                      DDSDDE(NTENS,NTENS), DDSDDT(NTENS), DRPLDE(NTENS), &
                      STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1), DPRED(1), &
                      PROPS(NPROPS), COORDS(3), DROT(3,3), DFGRD0(3,3), DFGRD1(3,3), &
                      FIBORI(NELEM,4)

      REAL(KIND=8) :: SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP, &
                      DTEMP, PNEWDT, CELENT

      INTEGER :: NTERM

!     FLAGS
!      INTEGER :: FLAG1

!     UTILITY TENSORS
      REAL(KIND=8) :: UNIT2(NDI,NDI), UNIT4(NDI,NDI,NDI,NDI), &
                      UNIT4S(NDI,NDI,NDI,NDI), &
                      PROJE(NDI,NDI,NDI,NDI), PROJL(NDI,NDI,NDI,NDI)

!     KINEMATICS
      REAL(KIND=8) :: DISTGR(NDI,NDI), C(NDI,NDI), B(NDI,NDI), &
                      CBAR(NDI,NDI), BBAR(NDI,NDI), DISTGRINV(NDI,NDI), &
                      UBAR(NDI,NDI), VBAR(NDI,NDI), ROT(NDI,NDI), &
                      DFGRD1INV(NDI,NDI)
      REAL(KIND=8) :: DET, CBARI1, CBARI2

!     VOLUMETRIC CONTRIBUTION
      REAL(KIND=8) :: PKVOL(NDI,NDI), SVOL(NDI,NDI), &
                      CVOL(NDI,NDI,NDI,NDI), CMVOL(NDI,NDI,NDI,NDI)
      REAL(KIND=8) :: KBULK, PV, PPV, SSEV

!     ISOCHORIC CONTRIBUTION
      REAL(KIND=8) :: SISO(NDI,NDI), PKISO(NDI,NDI), PK2(NDI,NDI), &
                      CISO(NDI,NDI,NDI,NDI), CMISO(NDI,NDI,NDI,NDI), &
                      SFIC(NDI,NDI), CFIC(NDI,NDI,NDI,NDI), &
                      PKFIC(NDI,NDI), CMFIC(NDI,NDI,NDI,NDI)

!     ISOCHORIC ISOTROPIC CONTRIBUTION
      REAL(KIND=8) :: C10, C01, SSEISO, DISO(5), PKMATFIC(NDI,NDI), &
                      SMATFIC(NDI,NDI), SISOMATFIC(NDI,NDI), &
                      CMISOMATFIC(NDI,NDI,NDI,NDI), &
                      CISOMATFIC(NDI,NDI,NDI,NDI)
      REAL(KIND=8) :: VORIF(3), VD(3), M0(3,3), MM(3,3), &
                      VORIF2(3), VD2(3), N0(3,3), NN(3,3)

!     LIST VARS OF OTHER CONTRIBUTIONS HERE

!     JAUMMAN RATE CONTRIBUTION (REQUIRED FOR ABAQUS UMAT)
      REAL(KIND=8) :: CJR(NDI,NDI,NDI,NDI)

!     CAUCHY STRESS AND ELASTICITY TENSOR
      REAL(KIND=8) :: SIGMA(NDI,NDI), DDSIGDDE(NDI,NDI,NDI,NDI), &
                      DDPKDDE(NDI,NDI,NDI,NDI)

!     TESTING/DEBUG VARS
      REAL(KIND=8) :: STEST(NDI,NDI), CTEST(NDI,NDI,NDI,NDI)
      INTEGER :: I1, J1, K1, L1
!----------------------------------------------------------------------
!-------------------------- INITIALIZATIONS ---------------------------
!----------------------------------------------------------------------
!     IDENTITY AND PROJECTION TENSORS
      UNIT2 = 0.0
      UNIT4 = 0.0
      UNIT4S = 0.0
      PROJE = 0.0
      PROJL = 0.0

!     KINEMATICS
      DISTGR = 0.0
      C = 0.0
      B = 0.0
      CBAR = 0.0
      BBAR = 0.0
      UBAR = 0.0
      VBAR = 0.0
      ROT = 0.0
      DET = 0.0
      CBARI1 = 0.0
      CBARI2 = 0.0

!     VOLUMETRIC
      PKVOL = 0.0
      SVOL = 0.0
      CVOL = 0.0
      KBULK = 0.0
      PV = 0.0
      PPV = 0.0
      SSEV = 0.0

!     ISOCHORIC
      SISO = 0.0
      PKISO = 0.0
      PK2 = 0.0
      CISO = 0.0
      CFIC = 0.0
      SFIC = 0.0
      PKFIC = 0.0

!     ISOTROPIC
      C10 = 0.0
      C01 = 0.0
      SSEISO = 0.0
      DISO = 0.0
      PKMATFIC = 0.0
      SMATFIC = 0.0
      SISOMATFIC = 0.0
      CMISOMATFIC = 0.0
      CISOMATFIC = 0.0

!     INITIALIZE OTHER CONT HERE

!     JAUMANN RATE
      CJR = 0.0

!     TOTAL CAUCHY STRESS AND ELASTICITY TENSORS
      SIGMA = 0.0
      DDSIGDDE = 0.0
!----------------------------------------------------------------------
!------------------------ IDENTITY TENSORS ----------------------------
!----------------------------------------------------------------------
      CALL ONEM(UNIT2, UNIT4, UNIT4S, NDI)
!----------------------------------------------------------------------
!------------------- MATERIAL CONSTANTS AND DATA ----------------------
!----------------------------------------------------------------------
!     VOLUMETRIC
      KBULK = PROPS(1)
!     ISOCHORIC ISOTROPIC NEO HOOKE
      C10 = PROPS(2)
!     NUMERICAL COMPUTATIONS
      NTERM = 60
!
!     STATE VARIABLES
!
      IF ((TIME(1) == 0.0) .AND. (KSTEP == 1)) THEN
         CALL INITIALIZE(STATEV)
      END IF
!        READ STATEV
      CALL SDVREAD(STATEV)
!      
!----------------------------------------------------------------------
!---------------------------- KINEMATICS ------------------------------
!----------------------------------------------------------------------
!     DISTORTION GRADIENT
      CALL FSLIP(DFGRD1, DISTGR, DET, NDI)
!     INVERSE OF DISTORTION GRADIENT
      CALL MATINV3D(DFGRD1, DFGRD1INV, NDI)
!     INVERSE OF DISTORTION GRADIENT
      CALL MATINV3D(DISTGR, DISTGRINV, NDI)
!     CAUCHY-GREEN DEFORMATION TENSORS
      CALL DEFORMATION(DFGRD1, C, B, NDI)
      CALL DEFORMATION(DISTGR, CBAR, BBAR, NDI)      
!     INVARIANTS OF DEVIATORIC DEFORMATION TENSORS
      CALL INVARIANTS(CBAR, CBARI1, CBARI2, NDI)
!     STRETCH TENSORS
      CALL STRETCH(CBAR, BBAR, UBAR, VBAR, NDI)
!     ROTATION TENSORS
      CALL ROTATION(DISTGR, ROT, UBAR, NDI)
!     DEVIATORIC PROJECTION TENSORS
      CALL PROJEUL(UNIT2, UNIT4S, PROJE, NDI)
!
      CALL PROJLAG(C, UNIT4, PROJL, NDI)
!----------------------------------------------------------------------
!--------------------- CONSTITUTIVE RELATIONS  ------------------------
!----------------------------------------------------------------------
!
!---- VOLUMETRIC ------------------------------------------------------
!     STRAIN-ENERGY AND DERIVATIVES (CHANGE HERE ACCORDING TO YOUR MODEL)
      CALL VOL(SSEV, PV, PPV, KBULK, DET)
      CALL ISOMAT(SSEISO, DISO, C10, CBARI1)
!
!---- ISOCHORIC ISOTROPIC ---------------------------------------------
!     PK2 'FICTICIOUS' STRESS TENSOR
      CALL PK2ISOMATFIC(PKMATFIC, DISO, CBAR, CBARI1, UNIT2, NDI)
!     CAUCHY 'FICTICIOUS' STRESS TENSOR
      CALL SIGISOMATFIC(SISOMATFIC, PKMATFIC, DISTGR, DET, NDI)
!     'FICTICIOUS' MATERIAL ELASTICITY TENSOR
      CALL CMATISOMATFIC(CMISOMATFIC, CBAR, CBARI1, CBARI2, DISO, UNIT2, UNIT4, DET, NDI)
!     'FICTICIOUS' SPATIAL ELASTICITY TENSOR
      CALL CSISOMATFIC(CISOMATFIC, CMISOMATFIC, DISTGR, DET, NDI)
!
!----------------------------------------------------------------------
!     SUM OF ALL ELASTIC CONTRIBUTIONS
!----------------------------------------------------------------------
!     STRAIN-ENERGY
      SSE = SSEV + SSEISO
!     PK2 'FICTICIOUS' STRESS
      PKFIC = PKMATFIC
!     CAUCHY 'FICTICIOUS' STRESS
      SFIC = SISOMATFIC
!     MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
      CMFIC = CMISOMATFIC
!     SPATIAL 'FICTICIOUS' ELASTICITY TENSOR
      CFIC = CISOMATFIC
!
!----------------------------------------------------------------------
!-------------------------- STRESS MEASURES ---------------------------
!----------------------------------------------------------------------
!
!---- VOLUMETRIC ------------------------------------------------------
!      PK2 STRESS
      CALL PK2VOL(PKVOL, PV, C, NDI)
!      CAUCHY STRESS
      CALL SIGVOL(SVOL, PV, UNIT2, NDI)
!
!---- ISOCHORIC -------------------------------------------------------
!      PK2 STRESS
      CALL PK2ISO(PKISO, PKFIC, PROJL, DET, NDI)
!      CAUCHY STRESS
      CALL SIGISO(SISO, SFIC, PROJE, NDI)
!
!---- VOLUMETRIC + ISOCHORIC ------------------------------------------
!      PK2 STRESS
      PK2 = PKVOL + PKISO
!      CAUCHY STRESS
      SIGMA = SVOL + SISO
!----------------------------------------------------------------------
!-------------------- MATERIAL ELASTICITY TENSOR ----------------------
!----------------------------------------------------------------------
!
!---- VOLUMETRIC ------------------------------------------------------
!
      CALL METVOL(CMVOL, C, PV, PPV, DET, NDI)
!
!---- ISOCHORIC -------------------------------------------------------
!
      CALL METISO(CMISO, CMFIC, PROJL, PKISO, PKFIC, C, UNIT2, DET, NDI)
!
!----------------------------------------------------------------------
!
      DDPKDDE = CMVOL + CMISO
!
!----------------------------------------------------------------------
!--------------------- SPATIAL ELASTICITY TENSOR ----------------------
!----------------------------------------------------------------------
!
!---- VOLUMETRIC ------------------------------------------------------
!
      CALL SETVOL(CVOL, PV, PPV, UNIT2, UNIT4S, NDI)
!
!---- ISOCHORIC -------------------------------------------------------
!
      CALL SETISO(CISO, CFIC, PROJE, SISO, SFIC, UNIT2, NDI)
!
!-----JAUMMAN RATE ----------------------------------------------------
!
      CALL SETJR(CJR, SIGMA, UNIT2, NDI)
!
!----------------------------------------------------------------------
!
!     ELASTICITY TENSOR
      !DDSIGDDE = CVOL + CISO + CJR
      DDSIGDDE = CVOL + CISO
!
!----------------------------------------------------------------------
!------------------------- INDEX ALLOCATION ---------------------------
!----------------------------------------------------------------------
!     VOIGT NOTATION  - FULLY SYMMETRY IMPOSED
      CALL INDEXX(STRESS, DDSDDE, SIGMA, DDSIGDDE, NTENS, NDI)
!
!----------------------------------------------------------------------
!--------------------------- STATE VARIABLES --------------------------
!----------------------------------------------------------------------
!     DO K1 = 1, NTENS
!      STATEV(1:27) = VISCOUS TENSORS
       CALL SDVWRITE(DET, STATEV)
!     END DO
!----------------------------------------------------------------------
      !write(*,*) 'F0'
      !write(*,*) DFGRD0
      !write(*,*) 'F1'
      !write(*,*) DFGRD1(1,1), DFGRD1(1,2), DFGRD1(1,3)
      !write(*,*) DFGRD1(2,1), DFGRD1(2,2), DFGRD1(2,3)
      !write(*,*) DFGRD1(3,1), DFGRD1(3,2), DFGRD1(3,3)
      !write(*,*) 'Cauchy stress'
      !write(*,*) STRESS(1), STRESS(2), STRESS(3), STRESS(4), STRESS(5), STRESS(6)
      RETURN
      END
!----------------------------------------------------------------------
!--------------------------- END OF UMAT ------------------------------
!----------------------------------------------------------------------
!
