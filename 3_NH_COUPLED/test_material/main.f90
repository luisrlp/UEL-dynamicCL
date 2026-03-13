PROGRAM TEST_GEL_MATERIAL
    IMPLICIT NONE

    ! Array dimensions
    INTEGER, PARAMETER :: NPROPS = 7, NDI = 3, NSHR = 3
    INTEGER, PARAMETER :: NTENS = 6, NSTATEV = 10
    INTEGER :: I, J, K, L

    ! Dummy UMAT variables required for the interface
    INTEGER :: NOEL, NPT, KSTEP, KINC, LAYER, KSPT
    REAL(8) :: STRESS(NTENS), STATEV(NSTATEV), DDSDDE(NTENS,NTENS), DDSDDT(NTENS)
    REAL(8) :: DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS), TIME(2)
    REAL(8) :: PREDEF(1), DPRED(1), COORDS(3), DROT(3,3), DFGRD0(3,3)
    REAL(8) :: PNEWDT, DET, CELENT

    ! Input Variables
    REAL(8) :: PROPS(NPROPS), DFGRD1(3,3)
    REAL(8) :: DTIME, MU_TAU, PHI_T, THETA, RT
    
    ! Output Variables
    REAL(8) :: SIGMA(3,3), DDSIGDDE(3,3,3,3)
    REAL(8) :: SPUCMOD(3,3), SPCUMODFAC(3,3)
    REAL(8) :: PHI_TAU, DPDT, DPHIDMU, DPHIDOTDMU
    REAL(8) :: MFLUID, DMDMU, DMDJ, VMOL_OUT

    !--------------------------------------------------------------
    ! 1. INITIALIZE ARRAYS TO ZERO
    !--------------------------------------------------------------
    SIGMA = 0.0D0
    DDSIGDDE = 0.0D0
    SPUCMOD = 0.0D0
    SPCUMODFAC = 0.0D0
    DFGRD1 = 0.0D0
    
    ! Initialize dummy UMAT variables
    STRESS = 0.0D0
    STATEV = 0.0D0
    DDSDDE = 0.0D0
    DFGRD0 = 0.0D0
    TIME = 0.0D0
    NOEL = 1
    NPT = 1
    KSTEP = 1
    KINC = 1

    !--------------------------------------------------------------
    ! 2. DEFINE MATERIAL PROPERTIES (GEL IN SOLVENT)
    !--------------------------------------------------------------
    PROPS(1) = 1.0D3   ! Kbulk
    PROPS(2) = 1.D0   ! C10 (GSHEAR/2)
    PROPS(3) = 0.1D0   ! chi (Flory-Huggins interaction parameter)
    PROPS(4) = 5.0D-9  ! D (Diffusion coefficient)
    PROPS(5) = 0.0D0   ! mu0 (Reference chemical potential)
    PROPS(6) = 1.D-4   ! Vmol (Molar volume of fluid, e.g., water)
    PROPS(7) = 8.314D0 ! Rgas (Universal Gas Constant)

    !--------------------------------------------------------------
    ! 3. DEFINE THERMODYNAMIC STATE AND KINEMATICS
    !--------------------------------------------------------------
    DTIME  = 1.0D0     ! Time increment
    PHI_T  = 0.999999D0   ! Previous polymer volume fraction
    THETA  = 298.D0    ! Temperature in Kelvin

    ! Compute thermal energy (R * T)
    RT = PROPS(7) * THETA

    ! Dynamically compute the initial chemical potential (MU_TAU)
    ! assuming an initially undeformed state (J = 1.0)
    MU_TAU = PROPS(5) + RT * (LOG(1.0D0 - PHI_T) + PHI_T + PROPS(3) * PHI_T**2) &
             - (PROPS(2) * PROPS(6)) * LOG(1.0D0) &
             + (PROPS(2) * PROPS(6) / (2.0D0 * RT)) * (LOG(1.0D0)**2)

    WRITE(*,*) 'Calculated Initial MU_TAU =', MU_TAU

    ! Tensile monotonic load test with shear component
    ! DFGRD1(1,1) = 1.3D0
    DFGRD1(1,1) = 1.1D0
    DFGRD1(1,2) = 0.3D0
    DFGRD1(1,3) = 0.0D0
    DFGRD1(2,1) = 0.0D0
    ! DFGRD1(2,2) = 1.0D0 / SQRT(DFGRD1(1,1))
    DFGRD1(2,2) = 1.0D0
    DFGRD1(2,3) = 0.0D0
    DFGRD1(3,1) = 0.0D0
    DFGRD1(3,2) = 0.0D0
    DFGRD1(3,3) = 1.0D0
    ! DFGRD1(3,3) = 1.0D0 / SQRT(DFGRD1(1,1)) 

    !--------------------------------------------------------------
    ! 4. CALL THE COUPLED MATERIAL SUBROUTINE
    !--------------------------------------------------------------
    WRITE(*,*) '===================================================='
    WRITE(*,*) '  TESTING COUPLED GEL MATERIAL SUBROUTINE           '
    WRITE(*,*) '===================================================='
    WRITE(*,*) 'PROPS: ', PROPS(1:NPROPS)
    WRITE(*,*) 'Running MATERIAL subroutine...'
    
    CALL MATERIAL(SIGMA, STATEV, DDSIGDDE, DFGRD0, DFGRD1, DET, &
                  TIME, DTIME, PREDEF, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, COORDS, &
                  PNEWDT, NOEL, NPT, KSTEP, KINC, MU_TAU, PHI_T, THETA, PHI_TAU, DPDT, &
                  DPHIDMU, DPHIDOTDMU, MFLUID, DMDMU, DMDJ, VMOL_OUT, SPUCMOD, SPCUMODFAC)
               
    WRITE(*,*) 'Done.'
    WRITE(*,*) '----------------------------------------------------'

    !--------------------------------------------------------------
    ! 5. PRINT RESULTS
    !--------------------------------------------------------------
    WRITE(*,*) '-> CAUCHY STRESS (SIGMA):'
    WRITE(*,'(3(ES12.4,1X))') SIGMA(1,1), SIGMA(1,2), SIGMA(1,3)
    WRITE(*,'(3(ES12.4,1X))') SIGMA(2,1), SIGMA(2,2), SIGMA(2,3)
    WRITE(*,'(3(ES12.4,1X))') SIGMA(3,1), SIGMA(3,2), SIGMA(3,3)
    WRITE(*,*) ''
    
    WRITE(*,*) '-> TANGENT MODULUS (DDSIGDDE):'
    WRITE(*,*) DDSIGDDE
    ! WRITE(*,'(3(3(ES12.4,1X))))') DDSIGDDE(1,1,:,:), DDSIGDDE(1,2,:,:), DDSIGDDE(1,3,:,:)
    ! WRITE(*,'(3(3(ES12.4,1X))))') DDSIGDDE(2,1,:,:), DDSIGDDE(2,2,:,:), DDSIGDDE(2,3,:,:)
    ! WRITE(*,'(3(3(ES12.4,1X))))') DDSIGDDE(3,1,:,:), DDSIGDDE(3,2,:,:), DDSIGDDE(3,3,:,:)
    ! WRITE(*,*) ''

    WRITE(*,*) '-> COUPLED CHEMISTRY VARIABLES:'
    WRITE(*,'(A,ES12.4)') ' New Polymer Vol Fraction (PHI_TAU) = ', PHI_TAU
    WRITE(*,'(A,ES12.4)') ' Time Rate of Swelling (DPDT)       = ', DPDT
    WRITE(*,'(A,ES12.4)') ' Scalar Fluid Mobility (MFLUID)     = ', MFLUID
    WRITE(*,*) ''

    WRITE(*,*) '-> DISPLACEMENT - CHEMICAL POTENTIAL MODULUS (SPUCMOD):'
    WRITE(*,'(3(ES12.4,1X))') SPUCMOD(1,1), SPUCMOD(1,2), SPUCMOD(1,3)
    WRITE(*,'(3(ES12.4,1X))') SPUCMOD(2,1), SPUCMOD(2,2), SPUCMOD(2,3)
    WRITE(*,'(3(ES12.4,1X))') SPUCMOD(3,1), SPUCMOD(3,2), SPUCMOD(3,3)
    WRITE(*,*) ''

    WRITE(*,*) '-> CHEM POTENTIAL - DISPLACEMENT MODULUS (SPCUMODFAC):'
    WRITE(*,'(3(ES12.4,1X))') SPCUMODFAC(1,1), SPCUMODFAC(1,2), SPCUMODFAC(1,3)
    WRITE(*,'(3(ES12.4,1X))') SPCUMODFAC(2,1), SPCUMODFAC(2,2), SPCUMODFAC(2,3)
    WRITE(*,'(3(ES12.4,1X))') SPCUMODFAC(3,1), SPCUMODFAC(3,2), SPCUMODFAC(3,3)
    WRITE(*,*) '===================================================='

END PROGRAM TEST_GEL_MATERIAL