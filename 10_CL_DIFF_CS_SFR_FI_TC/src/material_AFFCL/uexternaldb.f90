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

! DOUBLE PRECISION :: prefdir(nelem,4)
DOUBLE PRECISION :: prefdir(1,4)
DOUBLE PRECISION :: etadir(nelem*8, 2+ndir)
CHARACTER (LEN=256) ::  filename, jobdir, etafile
INTEGER :: lenjobdir,i,j,k
INTEGER :: alloc_stat

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

  ! globalSdv (UEL integration-point state, read by UVARM for the dummy mesh) allocated once here, at the start of
  ! the analysis and before any element call: with element operations in threads, the first UEL calls would
  ! otherwise all find it unallocated and allocate it at the same time (race; lost copies, NaN in the first frame)
  IF (.NOT. ALLOCATED(globalSdv)) THEN
    ALLOCATE(globalSdv(numElem, ngp, nsdv), STAT=alloc_stat)
    IF (alloc_stat /= 0) THEN
      write(*,*) 'UEXTERNALDB: error when allocating globalSdv, stat =', alloc_stat, ', numElem =', numElem
      call exit
    END IF
    globalSdv = 0.0d0
    write(*,*) 'UEXTERNALDB: globalSdv allocated, (numElem, ngp, nsdv) =', numElem, ngp, nsdv
  END IF

END IF

RETURN

END SUBROUTINE uexternaldb
