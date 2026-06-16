PROGRAM TEST_GENERAL_UMAT
    use global
    use,intrinsic :: ISO_Fortran_env
    
    PARAMETER(NTENS = 6, NPROPS = 30, NDI=3, NSHR=3)
    PARAMETER(NOEL = 1, NPT = 1)
    
    CHARACTER*8 CMNAME
    ! ARRAYS FOR UMAT SIGNATURE
    REAL(8) :: STRESS(NDI,NDI),STATEV(NSDV),DDSDDE(NDI,NDI,NDI,NDI)
    REAL(8) :: TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3)
    REAL(8) :: DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
    REAL(8) :: DET, DTIME, PNEWDT
    INTEGER :: KSTEP, KINC
    
    ! CHEMISTRY VARIABLES
    REAL(8) :: MU_TAU, PHI_T, THETA, PHI_TAU, DPDT, DPHIDMU, DPHIDOTDMU
    REAL(8) :: MFLUID, DMDMU, DMDJ, VMOL
    REAL(8) :: DMUDX(3,1), SPUCMOD(3,3), SPCUMODFAC(3,3)
    
    ! INITIALIZE ALL ARRAYS TO ZERO
    DDSDDE=0.0D0
    STRESS=0.0D0
    STATEV=0.0D0
    DFGRD0=0.0D0
    DFGRD1=0.0D0
    DMUDX = 0.0D0
    SPUCMOD = 0.0D0
    SPCUMODFAC = 0.0D0
    
    ! IDENTITY DEFORMATION
    DFGRD0(1,1)=1.0D0; DFGRD0(2,2)=1.0D0; DFGRD0(3,3)=1.0D0
    DFGRD1(1,1)=1.0D0; DFGRD1(2,2)=1.0D0; DFGRD1(3,3)=1.0D0
    DET = 1.0D0
    
    TIME(1)=0.d0
    TIME(2)=0.d0
    DTIME = 0.1d0
    KSTEP = 1
    KINC = 1
    
    ! CHEMICAL INPUTS (Dummy values for testing)
    MU_TAU = -10.0D0
    PHI_T  = 0.5D0
    THETA  = 293.0D0  ! Kelvin
    PHI_TAU = 0.5D0
    VMOL = 1.0d5
    
    ! MATERIAL PROPERTIES
    PROPS=0.0D0
    PROPS(1) = 1000.d0   ! K
    PROPS(2) = 1.00d0    ! C10
    PROPS(3) = 0.0D0     ! C01
    PROPS(4) = 0.5D0     ! phinet
    PROPS(5) = 1.0D0     ! a
    PROPS(6) = 1.0D0     ! r0c
    PROPS(7) = 0.5D0     ! etac
    PROPS(8) = 1.0D0     ! mu0str
    PROPS(9) = 1.0D0     ! beta
    PROPS(10) = 1.0D0    ! Lp
    PROPS(11) = 293.0D0  ! theta (temperature)
    PROPS(12) = 0.1D0    ! dx
    PROPS(13) = 1.0D0    ! bb
    PROPS(14) = 1.0D0    ! lambda0
    PROPS(15) = 100.0D0  ! cactin
    PROPS(16) = 0.1D0    ! R
    PROPS(17) = 0.9D0    ! Rfmax
    PROPS(18) = 0.25D0   ! Rbbmax
    PROPS(19) = 0.1D0    ! CHI
    PROPS(20) = 0.01D0   ! D
    PROPS(21) = -10.0D0  ! MU0
    PROPS(22) = 1.0d5    ! VMOL
    PROPS(23) = 1.0D0    ! Koff0
    PROPS(24) = 2.0D0    ! Keq
    
    ! APPLY DEFORMATION
    DFGRD1(1,1)=  1.3D0
    DFGRD1(1,2)=  0.0D0
    DFGRD1(1,3)=  0.0d0
    DFGRD1(2,1)=  0.0d0
    DFGRD1(2,2)=  1/sqrt(DFGRD1(1,1))
    DFGRD1(2,3)=  0.0d0
    DFGRD1(3,1)=  0.0d0
    DFGRD1(3,2)=  0.0d0
    DFGRD1(3,3)=  1/sqrt(DFGRD1(1,1))
    
    CALL MATERIAL(STRESS,STATEV,DDSDDE,DFGRD0,DFGRD1,DET,TIME,DTIME,PREDEF,NDI, &
         NSHR,NTENS,NSDV,PROPS,NPROPS,COORDS,PNEWDT,NOEL,NPT,KSTEP,KINC,MU_TAU,PHI_T, &
         PHI_TAU,DPDT,DPHIDMU,DPHIDOTDMU,MFLUID,DMDMU,DMUDX,DMDJ,VMOL,SPUCMOD,SPCUMODFAC)
    
    ! Write output to file
    OPEN(UNIT=10, FILE='results.txt', STATUS='REPLACE', ACTION='WRITE')
    
    write(10,*) '=== STRESS TENSOR ==='
    write(10,*) STRESS(1,1), STRESS(1,2), STRESS(1,3)
    write(10,*) STRESS(2,1), STRESS(2,2), STRESS(2,3)
    write(10,*) STRESS(3,1), STRESS(3,2), STRESS(3,3)
    write(10,*)
    write(10,*) '=== CHEMICAL OUTPUTS ==='
    write(10,*) 'PHI_TAU    = ', PHI_TAU
    write(10,*) 'DPDT       = ', DPDT
    write(10,*) 'DPHIDMU    = ', DPHIDMU
    write(10,*) 'DPHIDOTDMU = ', DPHIDOTDMU
    write(10,*) 'MFLUID     = ', MFLUID
    write(10,*) 'DMDMU      = ', DMDMU
    write(10,*) 'DMDJ       = ', DMDJ
    write(10,*) 'VMOL       = ', VMOL
    write(10,*)
    write(10,*) '=== STATE VARIABLES (first 15) ==='
    DO i = 1, 15
        write(10,'(A,I3,A,E15.6)') '  STATEV(', i, ') = ', STATEV(i)
    END DO
    write(10,*)
    write(10,*) 'DONE.'
    
    CLOSE(10)
    write(*,*) 'Output written to results.txt'
    
END PROGRAM