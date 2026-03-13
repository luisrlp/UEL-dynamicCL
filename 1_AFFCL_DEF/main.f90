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

      PARAMETER (NDOFEL=24, MLVARX=24, NRHS=1,NSVARS=8*NSDV,NPROPS=3)
      PARAMETER (NJPROP=2, MCRD=3,NNODE=8, JTYPE=3,KSTEP=1,KINC=1)
      PARAMETER (JELEM=1, NDLOAD=0,MDLOAD=0,NPREDF=1)

!     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD   

      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),PROPS(NPROPS),  &
      SVARS(NSVARS),ENERGY(8),COORDS(MCRD,NNODE),Uall(NDOFEL),    &
      DUall(MLVARX,1),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(1), &
      JDLTYP(MDLOAD,1),ADLMAG(MDLOAD,1),DDLMAG(MDLOAD,1),         &
      PREDEF(2,NPREDF,NNODE),LFLAGS(4),JPROPS(NJPROP)

      ! Initialize LFLAGS
      !LFLAGS(1) = 64
      LFLAGS(1) = 1
      LFLAGS(2) = 1
      LFLAGS(3) = 0
      LFLAGS(4) = 0      
! MATERIAL PROPERTIES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! k PENALTY PARAMETER
            PROPS(1)=1000.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ISOTROPIC MATRIX PARAMS
! C10=
            PROPS(2)=1.00d0
! C01=
            PROPS(3)=0.00d0
! JPROPS is used to pass the number of properties       
            JPROPS(1)=NSDV        
            JPROPS(2)=NSVARS
      
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
!
      Uall=0.d0
      !        NODE 1
      DUall(1,1)=0.d0
      DUall(2,1)=0.d0
      DUall(3,1)=0.d0
      !        NODE 2      !
      DUall(4,1)=0.d0
      DUall(5,1)=0.d0
      DUall(6,1)=0.d0
      !        NODE 3      !
      DUall(7,1)=0.d0
      DUall(8,1)=0.d0
      DUall(9,1)=0.d0
      !        NODE 4      !
      DUall(10,1)=0.d0
      DUall(11,1)=0.d0
      DUall(12,1)=0.d0
      !        NODE 5      !
      DUall(13,1)=0.3d0
      DUall(14,1)=0.d0
      DUall(15,1)=0.D0
      !        NODE 6      !
      DUall(16,1)=0.3d0
      DUall(17,1)=0.d0
      DUall(18,1)=0.D0
      !        NODE 7      !
      DUall(19,1)=0.3d0
      DUall(20,1)=0.d0
      DUall(21,1)=0.D0
      !        NODE 8      !
      DUall(22,1)=0.3d0
      DUall(23,1)=0.D0
      DUall(24,1)=0.D0
      
      Uall=Uall+DUall(:,1)
      !  
      rhs=0.d0
      svars=0.d0    
      Vel=0.d0
      Accn=0.d0
      TIME(1)=0.2d0
      DTIME=0.1d0
      TIME(2)=TIME(1) + DTIME
      !
      CALL UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,       &
      PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,  &
      TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,   &
      PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,   &
      NJPROP,PERIOD)
     
      write(*,*) AMATRX
      write(*,*) rhs

      END PROGRAM
