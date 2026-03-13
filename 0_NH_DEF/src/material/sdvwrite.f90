SUBROUTINE sdvwrite(det,statev)
! VISCOUS DISSIPATION: WRITE STATE VARS
    use global
    IMPLICIT NONE

    DOUBLE PRECISION STATEV(NSDV),DET
    !write your sdvs here. they should be allocated 
    !after the viscous terms (check hvwrite)
     STATEV(1)=DET
    RETURN

END SUBROUTINE sdvwrite
