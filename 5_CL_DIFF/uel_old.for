      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is used to transfer SDV's from the UEL
      !  onto the dummy mesh for viewing.  Note that an offset of
      !  ElemOffset is used between the real mesh and the dummy mesh.
      !  If your model has more than ElemOffset UEL elements, then
      !  this will need to be modified.
     
      use global
     
      include 'aba_param.inc'
      

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
      integer i1

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.

C      uvar(1) = globalSdv(noel-ElemOffset,npt,1)
c      for example
C      uvar(2) = globalSdv(noel-ElemOffset,npt,2)
      do i1=1,nsdv
      uvar(i1) = globalSdv(noel-ElemOffset,npt,i1)
      enddo
c      uvar(3) = globalSdv(noel-ElemOffset,npt,3)
c      uvar(4) = globalSdv(noel-ElemOffset,npt,4)

      return
      end subroutine uvarm
C>********************************************************************
C> Record of revisions:                                              |
C>        Date        Programmer        Description of change        |
C>        ====        ==========        =====================        |
C>     15/11/2017    Joao Ferreira      cont mech general eqs        |
C>     01/11/2018    Joao Ferreira      comments added               |
C>--------------------------------------------------------------------
C>     Description:
C>     UMAT: IMPLEMENTATION OF THE CONSTITUTIVE EQUATIONS BASED UPON 
C>     A STRAIN-ENERGY FUNCTION (SEF).
C>     THIS CODE, AS IS, EXPECTS A SEF BASED ON THE INVARIANTS OF THE 
C>     CAUCHY-GREEN TENSORS. A VISCOELASTIC COMPONENT IS ALSO 
C>     INCLUDED IF NEEDED. 
C>     YOU CAN CHOOSE TO COMPUTE AT THE MATERIAL FRAME AND THEN 
C>     PUSHFORWARD OR  COPUTE AND THE SPATIAL FRAME DIRECTLY.
C>--------------------------------------------------------------------
C>     IF YOU WANT TO ADAPT THE CODE ACCORDING TO YOUR SEF:
C>    ISOMAT - DERIVATIVES OF THE SEF IN ORDER TO THE INVARIANTS
C>    ADD OTHER CONTRIBUTIONS: STRESS AND TANGENT MATRIX
C>-------------------------------------------------------------------- 
C      STATE VARIABLES: CHECK ROUTINES - INITIALIZE, WRITESDV, READSDV.
C>--------------------------------------------------------------------              
C>     UEXTERNALDB: READ FILAMENTS ORIENTATION AND PREFERED DIRECTION
C>--------------------------------------------------------------------
C>---------------------------------------------------------------------
      SUBROUTINE MATERIAL(SIGMA,STATEV,DDSIGDDE,
     1 DFGRD0,DFGRD1,DET, 
     2 TIME,DTIME,PREDEF,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,PNEWDT,
     4 NOEL,NPT,KSTEP,KINC)
C
C----------------------------------------------------------------------
C--------------------------- DECLARATIONS -----------------------------
C----------------------------------------------------------------------
      use global
      IMPLICIT NONE
      
C     ADD COMMON BLOCKS HERE IF NEEDED (and in uexternal)
C      COMMON /KBLOCK/KBLOCK
C
      CHARACTER*8 CMNAME
C
      INTEGER NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT,
     1        LAYER, KSPT, KSTEP, KINC
C     
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 FIBORI(NELEM,4)
C
      DOUBLE PRECISION SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP,
     1                 DTEMP,PNEWDT,CELENT
C
      INTEGER NTERM
C
C     FLAGS
C      INTEGER FLAG1
C     UTILITY TENSORS
      DOUBLE PRECISION UNIT2(NDI,NDI),UNIT4(NDI,NDI,NDI,NDI),
     1                 UNIT4S(NDI,NDI,NDI,NDI),
     2                 PROJE(NDI,NDI,NDI,NDI),PROJL(NDI,NDI,NDI,NDI)
C     KINEMATICS
      DOUBLE PRECISION DISTGR(NDI,NDI),C(NDI,NDI),B(NDI,NDI),
     1                 CBAR(NDI,NDI),BBAR(NDI,NDI),DISTGRINV(NDI,NDI),
     2                 UBAR(NDI,NDI),VBAR(NDI,NDI),ROT(NDI,NDI),
     3                 DFGRD1INV(NDI,NDI)
      DOUBLE PRECISION DET,CBARI1,CBARI2
C     VOLUMETRIC CONTRIBUTION
      DOUBLE PRECISION PKVOL(NDI,NDI),SVOL(NDI,NDI),
     1                 CVOL(NDI,NDI,NDI,NDI),CMVOL(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION KBULK,PV,PPV,SSEV
C     ISOCHORIC CONTRIBUTION
      DOUBLE PRECISION SISO(NDI,NDI),PKISO(NDI,NDI),PK2(NDI,NDI),
     1                 CISO(NDI,NDI,NDI,NDI),CMISO(NDI,NDI,NDI,NDI),
     2                 SFIC(NDI,NDI),CFIC(NDI,NDI,NDI,NDI),
     3                 PKFIC(NDI,NDI),CMFIC(NDI,NDI,NDI,NDI)
C     ISOCHORIC ISOTROPIC CONTRIBUTION
      DOUBLE PRECISION C10,C01,SSEISO,DISO(5),PKMATFIC(NDI,NDI),
     1                 SMATFIC(NDI,NDI),SISOMATFIC(NDI,NDI),
     2                 CMISOMATFIC(NDI,NDI,NDI,NDI),
     3                 CISOMATFIC(NDI,NDI,NDI,NDI)   
      DOUBLE PRECISION VORIF(3),VD(3),M0(3,3),MM(3,3),
     1        VORIF2(3),VD2(3),N0(3,3),NN(3,3)
C     LIST VARS OF OTHER CONTRIBUTIONS HERE
C
C     JAUMMAN RATE CONTRIBUTION (REQUIRED FOR ABAQUS UMAT)
      DOUBLE PRECISION CJR(NDI,NDI,NDI,NDI)
C     CAUCHY STRESS AND ELASTICITY TENSOR
      DOUBLE PRECISION SIGMA(NDI,NDI),DDSIGDDE(NDI,NDI,NDI,NDI),
     1                                 DDPKDDE(NDI,NDI,NDI,NDI)
C     TESTING/DEBUG VARS
      DOUBLE PRECISION STEST(NDI,NDI), CTEST(NDI,NDI,NDI,NDI)
      INTEGER I1,J1,K1,L1
C----------------------------------------------------------------------
C-------------------------- INITIALIZATIONS ---------------------------
C----------------------------------------------------------------------
C     IDENTITY AND PROJECTION TENSORS
      UNIT2=ZERO
      UNIT4=ZERO
      UNIT4S=ZERO
      PROJE=ZERO
      PROJL=ZERO
C     KINEMATICS
      DISTGR=ZERO
      C=ZERO
      B=ZERO
      CBAR=ZERO
      BBAR=ZERO
      UBAR=ZERO
      VBAR=ZERO
      ROT=ZERO
      DET=ZERO
      CBARI1=ZERO
      CBARI2=ZERO
C     VOLUMETRIC
      PKVOL=ZERO
      SVOL=ZERO
      CVOL=ZERO
      KBULK=ZERO
      PV=ZERO
      PPV=ZERO
      SSEV=ZERO
C     ISOCHORIC
      SISO=ZERO
      PKISO=ZERO
      PK2=ZERO
      CISO=ZERO
      CFIC=ZERO
      SFIC=ZERO
      PKFIC=ZERO
C     ISOTROPIC
      C10=ZERO
      C01=ZERO
      SSEISO=ZERO
      DISO=ZERO
      PKMATFIC=ZERO
      SMATFIC=ZERO
      SISOMATFIC=ZERO
      CMISOMATFIC=ZERO
      CISOMATFIC=ZERO
C     INITIALIZE OTHER CONT HERE
C
C     JAUMANN RATE
      CJR=ZERO
C     TOTAL CAUCHY STRESS AND ELASTICITY TENSORS
      SIGMA=ZERO
      DDSIGDDE=ZERO
C
C----------------------------------------------------------------------
C------------------------ IDENTITY TENSORS ----------------------------
C----------------------------------------------------------------------
            CALL ONEM(UNIT2,UNIT4,UNIT4S,NDI)
C----------------------------------------------------------------------
C------------------- MATERIAL CONSTANTS AND DATA ----------------------
C----------------------------------------------------------------------
C     VOLUMETRIC
      KBULK    = PROPS(1)
C     ISOCHORIC ISOTROPIC NEO HOOKE
      C10      = PROPS(2)
C     NUMERICAL COMPUTATIONS
      NTERM    = 60
C
C     STATE VARIABLES
C
      IF ((TIME(1).EQ.ZERO).AND.(KSTEP.EQ.1)) THEN
      CALL INITIALIZE(STATEV)
      ENDIF
C        READ STATEV
      CALL SDVREAD(STATEV)
C      
C----------------------------------------------------------------------
C---------------------------- KINEMATICS ------------------------------
C----------------------------------------------------------------------
C     DISTORTION GRADIENT
      CALL FSLIP(DFGRD1,DISTGR,DET,NDI)
C     INVERSE OF DISTORTION GRADIENT
      CALL MATINV3D(DFGRD1,DFGRD1INV,NDI)
C     INVERSE OF DISTORTION GRADIENT
      CALL MATINV3D(DISTGR,DISTGRINV,NDI)
C     CAUCHY-GREEN DEFORMATION TENSORS
      CALL DEFORMATION(DFGRD1,C,B,NDI)
      CALL DEFORMATION(DISTGR,CBAR,BBAR,NDI)      
C     INVARIANTS OF DEVIATORIC DEFORMATION TENSORS
      CALL INVARIANTS(CBAR,CBARI1,CBARI2,NDI)
C     STRETCH TENSORS
      CALL STRETCH(CBAR,BBAR,UBAR,VBAR,NDI)
C     ROTATION TENSORS
      CALL ROTATION(DISTGR,ROT,UBAR,NDI)
C     DEVIATORIC PROJECTION TENSORS
      CALL PROJEUL(UNIT2,UNIT4S,PROJE,NDI)
C
      CALL PROJLAG(C,UNIT4,PROJL,NDI)
C----------------------------------------------------------------------
C--------------------- CONSTITUTIVE RELATIONS  ------------------------
C----------------------------------------------------------------------
C
C---- VOLUMETRIC ------------------------------------------------------
C     STRAIN-ENERGY AND DERIVATIVES (CHANGE HERE ACCORDING TO YOUR MODEL)
      CALL VOL(SSEV,PV,PPV,KBULK,DET)
      CALL ISOMAT(SSEISO,DISO,C10,CBARI1)
C
C---- ISOCHORIC ISOTROPIC ---------------------------------------------
C     PK2 'FICTICIOUS' STRESS TENSOR
      CALL PK2ISOMATFIC(PKMATFIC,DISO,CBAR,CBARI1,UNIT2,NDI)
C     CAUCHY 'FICTICIOUS' STRESS TENSOR
      CALL SIGISOMATFIC(SISOMATFIC,PKMATFIC,DISTGR,DET,NDI)
C     'FICTICIOUS' MATERIAL ELASTICITY TENSOR
      CALL CMATISOMATFIC(CMISOMATFIC,CBAR,CBARI1,CBARI2,
     1                          DISO,UNIT2,UNIT4,DET,NDI)
C     'FICTICIOUS' SPATIAL ELASTICITY TENSOR
      CALL CSISOMATFIC(CISOMATFIC,CMISOMATFIC,DISTGR,DET,NDI)
C
C----------------------------------------------------------------------
C     SUM OF ALL ELASTIC CONTRIBUTIONS
C----------------------------------------------------------------------
C     STRAIN-ENERGY
      SSE=SSEV+SSEISO
C     PK2 'FICTICIOUS' STRESS
      PKFIC=PKMATFIC
C     CAUCHY 'FICTICIOUS' STRESS
      SFIC=SISOMATFIC
C     MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
      CMFIC=CMISOMATFIC
C     SPATIAL 'FICTICIOUS' ELASTICITY TENSOR
      CFIC=CISOMATFIC
C
C----------------------------------------------------------------------
C-------------------------- STRESS MEASURES ---------------------------
C----------------------------------------------------------------------
C
C---- VOLUMETRIC ------------------------------------------------------
C      PK2 STRESS
      CALL PK2VOL(PKVOL,PV,C,NDI)
C      CAUCHY STRESS
      CALL SIGVOL(SVOL,PV,UNIT2,NDI)
C
C---- ISOCHORIC -------------------------------------------------------
C      PK2 STRESS
      CALL PK2ISO(PKISO,PKFIC,PROJL,DET,NDI)
C      CAUCHY STRESS
      CALL SIGISO(SISO,SFIC,PROJE,NDI)
C
C---- VOLUMETRIC + ISOCHORIC ------------------------------------------
C      PK2 STRESS
      PK2   = PKVOL + PKISO
C      CAUCHY STRESS
      SIGMA = SVOL  + SISO
C
C----------------------------------------------------------------------
C-------------------- MATERIAL ELASTICITY TENSOR ----------------------
C----------------------------------------------------------------------
C
C---- VOLUMETRIC ------------------------------------------------------
C
      CALL METVOL(CMVOL,C,PV,PPV,DET,NDI)
C
C---- ISOCHORIC -------------------------------------------------------
C
      CALL METISO(CMISO,CMFIC,PROJL,PKISO,PKFIC,C,UNIT2,DET,NDI)
C
C----------------------------------------------------------------------
C
      DDPKDDE=  CMVOL + CMISO
C
C----------------------------------------------------------------------
C--------------------- SPATIAL ELASTICITY TENSOR ----------------------
C----------------------------------------------------------------------
C
C---- VOLUMETRIC ------------------------------------------------------
C
      CALL SETVOL(CVOL,PV,PPV,UNIT2,UNIT4S,NDI)
C
C---- ISOCHORIC -------------------------------------------------------
C
      CALL SETISO(CISO,CFIC,PROJE,SISO,SFIC,UNIT2,NDI)
C
C-----JAUMMAN RATE ----------------------------------------------------
C
      CALL SETJR(CJR,SIGMA,UNIT2,NDI)
C
C----------------------------------------------------------------------
C
C     ELASTICITY TENSOR
      DDSIGDDE=CVOL+CISO+CJR
C
C----------------------------------------------------------------------
C------------------------- INDEX ALLOCATION ---------------------------
C----------------------------------------------------------------------
C     VOIGT NOTATION  - FULLY SIMMETRY IMPOSED
      CALL INDEXX(STRESS,DDSDDE,SIGMA,DDSIGDDE,NTENS,NDI)
C
C----------------------------------------------------------------------
C--------------------------- STATE VARIABLES --------------------------
C----------------------------------------------------------------------
C     DO K1 = 1, NTENS
C      STATEV(1:27) = VISCOUS TENSORS
       CALL SDVWRITE(DET,STATEV)
C     END DO
C----------------------------------------------------------------------
      RETURN
      END
C----------------------------------------------------------------------
C--------------------------- END OF UMAT ------------------------------
C----------------------------------------------------------------------
C

      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and 
      !   the dummy mesh needs to have the same number of integ
      !   points.  You must set that parameter value here.
      !
      !  ElemOffset
      !   Offset between element numbers on the real mesh and
      !    dummy mesh.  That is set in the input file, and 
      !    that value must be set here the same.

      integer numElem,ElemOffset,err
      INTEGER NWP,NELEM,NCH,NSDV
      DOUBLE PRECISION  ONE, TWO, THREE, FOUR, SIX, ZERO
      DOUBLE PRECISION HALF,THIRD
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of UEL elements used here
      parameter(numElem=1)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ElemOffset=1000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PARAMETER(NELEM=1, NSDV=1)
      PARAMETER(ZERO=0.D0, ONE=1.0D0,TWO=2.0D0)
      PARAMETER(THREE=3.0D0,FOUR=4.0D0,SIX=6.0D0)
      PARAMETER(HALF=0.5d0,THIRD=1.d0/3.d0)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable :: globalSdv(:,:,:)

      end module global
