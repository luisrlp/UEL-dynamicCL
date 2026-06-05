SUBROUTINE getprops_gp(noel, npt, etadir, etadir_array)

use global
IMPLICIT NONE

INTEGER, INTENT(IN)              :: noel, npt
DOUBLE PRECISION, INTENT(IN)     :: etadir(nelem*ngp, ndir+2)
DOUBLE PRECISION, INTENT(OUT)    :: etadir_array(ndir)
INTEGER                          :: i, l, idx

l = (noel-1)*NGP + npt
IF ((etadir(l,1)==noel).AND.(etadir(l,2)==npt)) THEN
  DO i = 1, ndir
    etadir_array(i) = etadir(l,i+2)
  END DO
ELSE
  DO i = 1, NELEM*NGP
    IF ((etadir(i,1)==noel).AND.(etadir(i,2)==npt)) THEN
      idx = i
      exit
  END IF
END DO
END IF

RETURN
END SUBROUTINE getprops_gp
