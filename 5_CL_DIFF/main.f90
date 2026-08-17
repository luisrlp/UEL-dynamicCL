PROGRAM TEST
      use, intrinsic :: ISO_Fortran_env
      use global
      implicit none

      !VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS

      REAL(8) :: RHS,AMATRX,SVARS,ENERGY

      !VARIABLES PASSED INTO UEL 
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME, &
      DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC, &
      JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP

      ! 8 nodes x 4 DOFs/node (3 disp + 1 chem potential) = 32
      PARAMETER (NDOFEL=32, MLVARX=32, NRHS=1,NSVARS=8*NSDV,NPROPS=24)
      PARAMETER (NJPROP=2, MCRD=3,NNODE=8, JTYPE=3)
      PARAMETER (JELEM=1, NDLOAD=0,MDLOAD=0,NPREDF=1)

      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),PROPS(NPROPS),  &
      SVARS(NSVARS),ENERGY(8),COORDS(MCRD,NNODE),Uall(NDOFEL),    &
      DUall(MLVARX,1),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(1), &
      JDLTYP(MDLOAD,1),ADLMAG(MDLOAD,1),DDLMAG(MDLOAD,1),         &
      PREDEF(2,NPREDF,NNODE),LFLAGS(4),JPROPS(NJPROP)

      ! TEST CONTROL VARIABLES
      LOGICAL :: DO_RELAXATION = .TRUE.
      REAL(8) :: T_SHEAR = 1.0d0
      REAL(8) :: T_RELAX = 10.0d0
      REAL(8) :: TARGET_SHEAR = 0.3d0
      REAL(8) :: DT_INIT = 0.1d0
      REAL(8) :: CURRENT_TIME = 0.0d0
      REAL(8) :: TOTAL_TIME
      REAL(8) :: D_SHEAR

      ! Initialize LFLAGS
      LFLAGS(1) = 1
      LFLAGS(2) = 1
      LFLAGS(3) = 0
      LFLAGS(4) = 0      

! MATERIAL PROPERTIES (must match properties.inp and _umat_.f90 ordering)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VOLUMETRIC
      PROPS(1)  = 1000.d0    ! K: Bulk modulus (penalty)
! ISOCHORIC ISOTROPIC
      PROPS(2)  = 0.0d0      ! C10: Mooney-Rivlin
      PROPS(3)  = 0.0d0      ! C01: Mooney-Rivlin
      PROPS(4)  = 1.0d0      ! phinet: Network volume fraction
! ACTIN/CROSSLINKERS
      PROPS(5)  = 1.2d0      ! a: Contour/end-to-end ratio
      PROPS(6)  = 0.014d0    ! R0C: Reference CL distance
      PROPS(7)  = 0.5d0      ! etac: CL stiffness
      PROPS(8)  = 38600.d0   ! mu0str: CL stretch modulus
      PROPS(9)  = 0.5d0      ! beta: CL stretch parameter
      PROPS(10) = 16.0d0     ! Lp: Persistence length
      PROPS(11) = 298.0d0    ! theta: Absolute temperature (25+273)
      PROPS(12) = 0.001d0    ! dx: CL reactive distance
! AFFINE NETWORK
      PROPS(13) = 0.001d0    ! bb: von Mises concentration
      PROPS(14) = 1.0d0      ! lambda0: Reference stretch
      PROPS(15) = 0.0095d0   ! cactin: Actin concentration
      PROPS(16) = 0.1d0      ! R: CL to actin ratio
      PROPS(17) = 1.0d0      ! Rfmax: Max free CL ratio
      PROPS(18) = 0.25d0     ! Rbmax: Max bound CL ratio
! SOLVENT / DIFFUSION
      PROPS(19) = 0.1d0      ! CHI: Flory-Huggins chi
      PROPS(20) = 0.1d0      ! D: Diffusion coefficient
      PROPS(21) = 0.0d0      ! MU0: Reference chemical potential
      PROPS(22) = 0.15d0     ! VMOL: Molar volume
      PROPS(23) = 0.05d0     ! Koff0: Baseline off-rate
      PROPS(24) = 0.25d0     ! Keq: Equilibrium constant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! JPROPS: integer properties (must match U3D8 expectations)
      JPROPS(1) = NSDV       ! nlSdv: local SDVs per integration point
      JPROPS(2) = NSDV       ! ngSdv: global SDVs per integration point
      
!  8-node     8-----------7
!  brick     /|          /|       zeta
!           / |         / |       
!          5-----------6  |       |     eta
!          |  |        |  |       |   /
!          |  |        |  |       |  /
!          |  4--------|--3       | /
!          | /         | /        |/
!          |/          |/         O--------- xi
!          1-----------2        origin at cube center
!
!     Face numbering follows:
!       Face 1 = nodes 1,2,3,4
!       Face 2 = nodes 5,8,7,6
!       Face 3 = nodes 1,5,6,2
!       Face 4 = nodes 2,6,7,3
!       Face 5 = nodes 3,7,8,4
!       Face 6 = nodes 4,8,5,1
!              NODE 1      
      coords(1,1)=0.d0
      coords(2,1)=0.d0
      coords(3,1)=0.d0
      !        NODE 2      
      coords(1,2)=1.d0
      coords(2,2)=0.d0
      coords(3,2)=0.d0
      !        NODE 3      
      coords(1,3)=1.d0
      coords(2,3)=1.d0
      coords(3,3)=0.d0
      !        NODE 4      
      coords(1,4)=0.d0
      coords(2,4)=1.d0
      coords(3,4)=0.d0
      !        NODE 5      
      coords(1,5)=0.d0
      coords(2,5)=0.d0
      coords(3,5)=1.d0
      !        NODE 6      
      coords(1,6)=1.d0
      coords(2,6)=0.d0
      coords(3,6)=1.d0
      !        NODE 7      
      coords(1,7)=1.d0
      coords(2,7)=1.d0
      coords(3,7)=1.d0
      !        NODE 8      
      coords(1,8)=0.d0
      coords(2,8)=1.d0
      coords(3,8)=1.d0

      ! Open an output file to record results
      open(unit=10, file='test_results.dat', status='replace')
      write(10, *) 'Time, u_shear, Force_shear'

      ! Global Initialization
      Uall = 0.d0
      svars = 0.d0
      Vel    = 0.d0
      Accn   = 0.d0
      PREDEF = 0.d0
      PERIOD = 1.d0
      
      IF (DO_RELAXATION) THEN
          TOTAL_TIME = T_SHEAR + T_RELAX
      ELSE
          TOTAL_TIME = T_SHEAR
      END IF

      KSTEP = 1
      KINC = 1
      CURRENT_TIME = 0.d0

      ! Initialize common blocks (prefdir, etadir)
      ! Note: uexternaldb expects current time. At initialization it's 0.
      TIME(1) = 0.0d0
      TIME(2) = 0.0d0
      DTIME = DT_INIT
      CALL uexternaldb(0, LFLAGS(1), TIME, DTIME, KSTEP, KINC)

      print *, "Starting Element-Level Test..."
      print *, "Do Relaxation = ", DO_RELAXATION

      ! ---------------------------------------------------------
      ! TIME STEPPING LOOP
      ! ---------------------------------------------------------
      DO WHILE (CURRENT_TIME < TOTAL_TIME - 1.d-8)
          DTIME = DT_INIT
          
          ! Check if we need to switch from shear to relaxation step
          IF ((CURRENT_TIME >= T_SHEAR) .AND. (KSTEP == 1)) THEN
              KSTEP = 2
              KINC = 1
          END IF

          TIME(1) = CURRENT_TIME
          TIME(2) = CURRENT_TIME + DTIME
          
          DUall = 0.d0

          IF (KSTEP == 1) THEN
              ! Shear Step
              D_SHEAR = (TARGET_SHEAR / T_SHEAR) * DTIME
          ELSE
              ! Relaxation Step
              D_SHEAR = 0.0d0
          END IF

          ! Apply D_SHEAR to top nodes (shear in X)
          DUall(17,1) = D_SHEAR    ! Node 5, u_x
          DUall(21,1) = D_SHEAR    ! Node 6, u_x
          DUall(25,1) = D_SHEAR    ! Node 7, u_x
          DUall(29,1) = D_SHEAR    ! Node 8, u_x
          
          Uall = Uall + DUall(:,1)
          
          rhs    = 0.d0
          ENERGY = 0.d0
          AMATRX = 0.d0

          ! CALL UEL to advance state over DTIME
          CALL UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,       &
          PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,  &
          TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,   &
          PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,   &
          NJPROP,PERIOD)
         
          ! RHS holds the residual forces. In a displacement-controlled test,
          ! the residual force at the top nodes corresponds to the reaction force.
          ! Let's sum the X-reaction forces on the top face (Nodes 5,6,7,8 -> DOFs 17,21,25,29).
          ! Note: RHS in UEL is the internal force vector (or residual). 
          ! The reaction force is typically the negative of the internal force for prescribed DOFs.
          
          write(10, '(3E15.6)') TIME(2), Uall(17), -(RHS(17,1)+RHS(21,1)+RHS(25,1)+RHS(29,1))
          
          CURRENT_TIME = CURRENT_TIME + DTIME
          KINC = KINC + 1
      END DO

      print *, "Test finished successfully. Results written to test_results.dat."
      close(10)

      END PROGRAM
