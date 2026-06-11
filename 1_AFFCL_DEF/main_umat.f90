PROGRAM TEST_GENERAL_UMAT
    use,intrinsic :: ISO_Fortran_env
    INCLUDE 'aba_param.inc'
    INCLUDE 'param_umat.inc'
    
    !C     ADD COMMON BLOCKS HERE IF NEEDED ()
    !C      COMMON /KBLOCK/KBLOCK
    
    PARAMETER(NTENS = 6, NSTATEV = NSDV, NPROPS = 14, NDI=3, NSHR=3)
    PARAMETER(NOEL = 1, NPT = 1)
    !
    CHARACTER*8 CMNAME
    DIMENSION STRESS(NDI, NDI),STATEV(NSTATEV),DDSDDE(NDI,NDI,NDI,NDI),DDSDDT(NTENS),      &
    DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),            &
    PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
    !
    INTEGER :: lop, lrestart
    REAL(4) :: time4(2)
    !
    i=1.0d0
    j=1.0d0
    DDSDDE=0.0D0
    DO i=1,NDI
        DO j=1,NDI
            STRESS(i,j)=0.0D0
        ENDDO
    ENDDO
    
    !
    ! DEFORMATION GRADIENT
    DFGRD1(1,1)= 1.1D0
    DFGRD1(1,2)= 0.0D0
    DFGRD1(1,3)= 0.0D0
    DFGRD1(2,1)= 0.0D0
    DFGRD1(2,2)= 1.0D0/DFGRD1(1,1)
    DFGRD1(2,3)= 0.0D0
    DFGRD1(3,1)= 0.0D0
    DFGRD1(3,2)= 0.0D0
    DFGRD1(3,3)= 1.0D0/DFGRD1(1,1)
    !
    time(1)=0.d0
    time(2)=0.d0
    dtime = 0.1d0
    kstep = 1
    !
    ! MATERIAL PROPERTIES
    !
    ! k PENALTY PARAMETER
PROPS(1)=1000.000d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ISOTROPIC MATRIX PARAMS
! C10=
PROPS(2)=1.00d0
! C01
PROPS(3)=1.00d0
!PHI....
PROPS(4)=1.0000d0 ! 1.0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SINGLE FILAMENT PARAMS
!L (default: 1.96)
PROPS(5)= 1.96d0
!CACTIN
! PROPS(5)=9.5d0
!R0F
PROPS(6)= 1.63d0
!R
! PROPS(6)=0.1d0
!R0C
PROPS(7) = 0.014d0
!ETAC
PROPS(8)= 0.5d0
!mu0
!PROPS(7)=38600.0d0
PROPS(9)= 38600.0d0
!beta
PROPS(10)=0.5d0
!PROPS(8)=0.5d0
!B0 = tk*lp*k0
! PROPS(11)=294.d0*16.d0*1.38d-5
PROPS(11) = 16.0
!lambda0.
PROPS(12)=1.00d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!AFFINE NETWORK PARAMS
!n - isotropic filaments per unit volume
PROPS(13)=7.66D0 ! 7.6627
!A
! PROPS(13)=1.2d0
!B....
PROPS(14)=0.001d0
    ! !
    STATEV=0.D0
    !
    erf=0.d0
    RHO=0.D0
    !
    !
     DFGRD1(1,1)=  1.0D0
     DFGRD1(1,2)=  0.0D0
     DFGRD1(1,3)=  0.0d0
     DFGRD1(2,1)=  0.0d0
     DFGRD1(2,2)=  1.0D0
     DFGRD1(2,3)=  0.0d0
     DFGRD1(3,1)=  0.0d0
     DFGRD1(3,2)=  0.0d0
     DFGRD1(3,3)=  1.0D0
    !
    !################################################################################################!
    !!     TENSILE MONOTONIC LOAD TEST
    !  DFGRD1(1,1)=  1.3D0
     DFGRD1(1,2)=  0.2D0
     DFGRD1(1,3)=  0.0d0
     DFGRD1(2,1)=  0.0d0
    !  DFGRD1(2,2)=  1/sqrt(DFGRD1(1,1))
     DFGRD1(2,3)=  0.0d0
     DFGRD1(3,1)=  0.0d0
     DFGRD1(3,2)=  0.0d0
    !  DFGRD1(3,3)=  1/sqrt(DFGRD1(1,1))
    !
DFGRD1(1,1)=  0.999409865502331
DFGRD1(1,2)=  6.840823995929941E-001
DFGRD1(1,3)=  1.415375547882356E-002
DFGRD1(2,1)=  -2.710760077512183E-01
DFGRD1(2,2)=  1.00134406533680
DFGRD1(2,3)=  -2.710760077512183E-02
DFGRD1(3,1)=  -3.978368599508040E-01
DFGRD1(3,2)=  -8.550969990486904E-05
DFGRD1(3,3)=  0.999249478911563

    
     lop = 0
     lrestart = 0
     time4(1) = real(time(1), 4)
     time4(2) = real(time(2), 4)
     CALL uexternaldb(lop,lrestart,time4,dtime,kstep,kinc)
     
     CALL MATERIAL(STRESS,STATEV,DDSDDE,DFGRD0,DFGRD1,DET,TIME,DTIME,PREDEF,NDI, &
     NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,PNEWDT,NOEL,NPT,KSTEP,KINC)
    !  CALL MATERIAL(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
    ! DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
    ! NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
    !
    
     write(*,*) 'STRESS'
     write(*,*) STRESS
    !  write(*,*) STRESS(1,1), STRESS(1,2), STRESS(1,3)
    !  write(*,*) STRESS(2,1), STRESS(2,2), STRESS(2,3)
    !  write(*,*) STRESS(3,1), STRESS(3,2), STRESS(3,3)
     write(*,*)
     write(*,*) 'DDSDDE'
     write(*,*) DDSDDE

    !################################################################################################!
    !
    END PROGRAM
    