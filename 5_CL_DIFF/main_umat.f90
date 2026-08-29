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
        REAL(8) :: MU_TAU, THETAF_T, THETA, THETAF_TAU, DTHETAFDT, DTHETAFDMU, RMACRO
        REAL(8) :: MFLUID, DMDMU, DMDJ, VMOL, CFMAX             
        REAL(8) :: DMUDX(3,1), DSIGDMU(3,3), SPCUMODFAC(3,3)    
        INTEGER :: i, NINC, lop, lrestart
    REAL(8) :: TFINAL, DEF_MIN, DEF_MAX, DEF
    REAL(4) :: time4(2)
    ! INITIALIZE ALL ARRAYS TO ZERO
    DDSDDE=0.0D0
    STRESS=0.0D0
    STATEV=0.0D0
    DFGRD0=0.0D0
    DFGRD1=0.0D0
    DMUDX = 0.0D0
    DSIGDMU = 0.0D0
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
    MU_TAU = -5885.596868D0
    THETAF_T  = 0.078734D0
    THETAF_TAU = THETAF_T
    
    ! MATERIAL PROPERTIES
    PROPS=0.0D0
    PROPS(1) = 1000.d0   ! K
    PROPS(2) = 1.00d0    ! C10
    PROPS(3) = 1.0D0     ! C01
    PROPS(4) = 1.D0     ! phinet
    PROPS(5) = 1.2D0     ! a
    PROPS(6) = 0.014D0     ! r0c
    PROPS(7) = two / three     ! etac
    PROPS(8) = 1.0D12     ! mu0str
    PROPS(9) = 0.5D0     ! beta
    PROPS(10) = 16.0D0    ! Lp
    PROPS(11) = 25.0D0 + 273.0D0  ! theta (temperature)
    PROPS(12) = 0.1D0    ! dx
    PROPS(13) = 0.000001D0    ! bb
    PROPS(14) = 1.0D0    ! lambda0
    PROPS(15) = 0.0095D0  ! cactin
    PROPS(16) = 0.1D0    ! R
    PROPS(17) = 1.0D0    ! Rfmax
    PROPS(18) = 0.25D0   ! Rbmax
    PROPS(19) = 0.1D0    ! CHI
    PROPS(20) = 0.1D0   ! D
    PROPS(21) = 0.0D0  ! MU0
    PROPS(22) = 0.15d0    ! VMOL
    PROPS(23) = 0.05D0    ! Koff0
    PROPS(24) = 0.25D0    ! Keq
    
    VMOL = PROPS(22)
    THETA = PROPS(11)
! =================================================
! SIMULATION LOOP
! =================================================
    TIME(1) = 0.0d0
    DTIME = 0.05d0
    TFINAL = 60.0d0
    NINC = INT(TFINAL/DTIME)
    DEF_MIN = 0.0d0
    DEF_MAX = 0.3d0

    OPEN(UNIT=10, FILE='results.txt', STATUS='REPLACE', ACTION='WRITE')
    write(10,*) 'INC, TIME, STRETCH, S11, S12, THETAF, CB_TOT'
    
    lop = 0
    lrestart = 0
    time4(1) = real(time(1), 4)
    time4(2) = real(time(2), 4)
    CALL uexternaldb(lop,lrestart,time4,dtime,kstep,kinc)
    
    DO KINC = 1, NINC

        DEF = DEF_MIN + (DEF_MAX - DEF_MIN) * REAL(KINC) / REAL(NINC)
        
        DFGRD0 = DFGRD1
        
        ! APPLY DEFORMATION
        DFGRD1(1,1)=  1.0d0
        DFGRD1(1,2)=  DEF
        DFGRD1(1,3)=  0.0d0
        DFGRD1(2,1)=  0.0d0
        DFGRD1(2,2)=  1.0d0 !1.0d0/sqrt(DFGRD1(1,1))
        DFGRD1(2,3)=  0.0d0
        DFGRD1(3,1)=  0.0d0
        DFGRD1(3,2)=  0.0d0
        DFGRD1(3,3)=  1.0d0 !/DFGRD1(1,1)**2

        TIME(2) = TIME(1) + DTIME

        write(*,*) 'Simulating Increment: ', KINC
        
        CALL MATERIAL(STRESS,STATEV,DDSDDE,DFGRD0,DFGRD1,DET,TIME,DTIME,PREDEF,NDI, &
            NSHR,NTENS,NSDV,PROPS,NPROPS,COORDS,PNEWDT,NOEL,NPT,KSTEP,KINC,MU_TAU,THETAF_T, &
            THETAF_TAU,DTHETAFDT,DTHETAFDMU,RMACRO,MFLUID,DMDMU,DMUDX,DMDJ,VMOL,CFMAX,DSIGDMU,SPCUMODFAC)
        
        ! Write output to file
        write(10,'(I4, 7(A,E15.6))') KINC, ', ', TIME(2), ',', DEF, ', ', STRESS(1,1) - STRESS(3,3), &
        ', ', STRESS(1,2), ', ', THETAF_TAU, ', ', STATEV(4)
    
        ! write(10,*) '=== STRESS TENSOR ==='
        ! write(10,*) STRESS(1,1), STRESS(1,2), STRESS(1,3)
        ! write(10,*) STRESS(2,1), STRESS(2,2), STRESS(2,3)
        ! write(10,*) STRESS(3,1), STRESS(3,2), STRESS(3,3)
        ! write(10,*)
        ! write(10,*) '=== CHEMICAL OUTPUTS ==='
        ! write(10,*) 'THETAF_TAU    = ', THETAF_TAU
        ! write(10,*) 'DTHETAFDT   = ', DTHETAFDT
        ! write(10,*) 'DTHETAFDMU  = ', DTHETAFDMU
        ! write(10,*) 'RMACRO     = ', RMACRO
        ! write(10,*) 'MFLUID     = ', MFLUID
        ! write(10,*) 'DMDMU      = ', DMDMU
        ! write(10,*) 'DMDJ       = ', DMDJ
        ! write(10,*) 'VMOL       = ', VMOL
        ! write(10,*)
        ! write(10,*) '=== STATE VARIABLES (first 15) ==='
        ! DO i = 1, 15
        !     write(10,'(A,I3,A,E15.6)') '  STATEV(', i, ') = ', STATEV(i)
        ! END DO
        ! write(10,*)
        ! write(10,*) 'DONE.'

        TIME(1) = TIME(2)
        THETAF_T = THETAF_TAU

    END DO

    CLOSE(10)
    write(*,*) 'Simulation complete. Output written to results.txt'
    
END PROGRAM