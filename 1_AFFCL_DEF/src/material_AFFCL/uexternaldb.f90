SUBROUTINE uexternaldb(lop,lrestart,time,dtime,kstep,kinc)



!>    READ FILAMENTS ORIENTATION AND PREFERED DIRECTIONS
use global
INCLUDE 'aba_param.inc'
!       this subroutine get the directions and weights for
!      the numerical integration

!     UEXTERNAL just called once; work in parallel computing

INTEGER, INTENT(IN OUT)                  :: lop
INTEGER, INTENT(IN OUT)                  :: lrestart
REAL, INTENT(IN OUT)                     :: time(2)
real(8), INTENT(IN OUT)                  :: dtime
INTEGER, INTENT(IN OUT)                  :: kstep
INTEGER, INTENT(IN OUT)                  :: kinc

COMMON /kfilp/prefdir
COMMON /kfile/etadir

DOUBLE PRECISION :: prefdir(nelem,4)
DOUBLE PRECISION :: etadir(nelem*8, 2+ndir)
CHARACTER (LEN=256) ::  filename, jobdir, etafile
INTEGER :: lenjobdir,i,j,k

!     LOP=0 --> START OF THE ANALYSIS
IF(lop == 0.OR.lop == 4) THEN
  
  CALL getoutdir(jobdir,lenjobdir)
  write(*,*) 'nelem = ', nelem
  write(*,*) 'jobdir = ', jobdir
  write(*,*) 'lenjobdir = ', lenjobdir
  
  !preferential direction
  filename=jobdir(:lenjobdir)//'/'//dir2
  OPEN(16,FILE=filename,STATUS='OLD')
  DO i=1,nelem
    READ(16,*) (prefdir(i,j),j=1,4)
  END DO
  CLOSE(16)

  !random CL stiffness eta
  !etafile = jobdir(:lenjobdir)//'/'//dir3
  !OPEN(17,FILE=etafile)
  !DO i=1,nelem*ngp
  !  READ(17,*) (etadir(i,j),j=1,ndir+2)
  !END DO
  !CLOSE(17)


END IF

RETURN

END SUBROUTINE uexternaldb
