SUBROUTINE pfdlambdadc(sfic,rho,lambda,unit2,m,rw,det,ndi)



!>    SINGLE FILAMENT:  PUSH FORWARD OF THE DERIVATIVE OF THE FILAMENT STRETCH WRT THE RIGHT CAUCHY-GREEN TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi !number of dimensions
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi) !ficticious cauchy stress 
DOUBLE PRECISION, INTENT(IN)             :: rho  !angular density at m
DOUBLE PRECISION, INTENT(IN)             :: lambda !filament stretch
DOUBLE PRECISION, INTENT(IN)             :: unit2(ndi,ndi) !inverse of modified right cauchy-green
DOUBLE PRECISION, INTENT(IN)             :: m(ndi) !direction vector
DOUBLE PRECISION, INTENT(IN)             :: rw ! integration weights
DOUBLE PRECISION, INTENT(IN)             :: det !determinant of the deformation gradient



INTEGER :: i1,j1

DOUBLE PRECISION :: aux

aux=rho * rw
DO i1=1,ndi
  DO j1=1,ndi
    sfic(i1,j1)=aux * (m(i1)*m(j1) - one/three * lambda**two * unit2(i1,j1))
  END DO
END DO

RETURN
END SUBROUTINE pfdlambdadc
