SUBROUTINE indexx(stress,ddsdde,sig,tng,ntens,ndi)



!>    INDEXATION: FULL SIMMETRY  IN STRESSES AND ELASTICITY TENSORS
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
INTEGER, INTENT(IN)                      :: ntens
DOUBLE PRECISION, INTENT(OUT)            :: stress(ntens)
DOUBLE PRECISION, INTENT(OUT)            :: ddsdde(ntens,ntens)
DOUBLE PRECISION, INTENT(IN)             :: sig(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: tng(ndi,ndi,ndi,ndi)



INTEGER :: ii1(6),ii2(6), i1,j1


DOUBLE PRECISION :: pp1,pp2

ii1(1)=1
ii1(2)=2
ii1(3)=3
ii1(4)=1
ii1(5)=1
ii1(6)=2

ii2(1)=1
ii2(2)=2
ii2(3)=3
ii2(4)=2
ii2(5)=3
ii2(6)=3

DO i1=1,ntens
!       STRESS VECTOR
  stress(i1)=sig(ii1(i1),ii2(i1))
  DO j1=1,ntens
!       DDSDDE - FULLY SIMMETRY IMPOSED
    pp1=tng(ii1(i1),ii2(i1),ii1(j1),ii2(j1))
    pp2=tng(ii1(i1),ii2(i1),ii2(j1),ii1(j1))
    ddsdde(i1,j1)=(one/two)*(pp1+pp2)
  END DO
END DO

RETURN

END SUBROUTINE indexx
SUBROUTINE pull2(pk,sig,finv,det,ndi)



!>       PULL-BACK TIMES DET OF A 2ND ORDER TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pk(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: sig(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: finv(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: det



INTEGER :: i1,j1,ii1,jj1


DOUBLE PRECISION :: aux


DO i1=1,ndi
  DO j1=1,ndi
    aux=zero
    DO ii1=1,ndi
      DO jj1=1,ndi
        aux=aux+det*finv(i1,ii1)*finv(j1,jj1)*sig(ii1,jj1)
      END DO
    END DO
    pk(i1,j1)=aux
  END DO
END DO

RETURN
END SUBROUTINE pull2
SUBROUTINE setjr(cjr,sigma,unit2,ndi)


use global
IMPLICIT NONE
!>    JAUMAN RATE CONTRIBUTION FOR THE SPATIAL ELASTICITY TENSOR

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cjr(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: sigma(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit2(ndi,ndi)



INTEGER :: i1,j1,k1,l1


DO i1 = 1, ndi
  DO j1 = 1, ndi
    DO k1 = 1, ndi
      DO l1 = 1, ndi
        
        cjr(i1,j1,k1,l1)= (one/two)*(unit2(i1,k1)*sigma(j1,l1)  &
            +sigma(i1,k1)*unit2(j1,l1)+unit2(i1,l1)*sigma(j1,k1)  &
            +sigma(i1,l1)*unit2(j1,k1))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE setjr
SUBROUTINE spectral(a,d,v)



!>    EIGENVALUES AND EIGENVECTOR OF A 3X3 MATRIX
!     THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
!     A SYMMETRIC 3X3 MATRIX A.

!     THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
!     EIGENVALUES IN ASCENDING ORDER, AND A MATRIX V WHOSE
!     COLUMNS CONTAIN THE CORRESPONDING EIGENVECTORS.

use global

DOUBLE PRECISION, INTENT(IN OUT)         :: a(3,3)
DOUBLE PRECISION                         :: e(3,3)
DOUBLE PRECISION, INTENT(IN OUT)         :: d(3)
DOUBLE PRECISION, INTENT(IN OUT)         :: v(3,3)

INTEGER :: nrot
INTEGER :: np=3



e = a

CALL jacobi(e,3,np,d,v,nrot)
CALL eigsrt(d,v,3,np)

RETURN
END SUBROUTINE spectral

!***********************************************************************

SUBROUTINE jacobi(a,n,np,d,v,nrot)

! COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
!  MATRIX A, WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL
!  NP BY NP ARRAY.  ON OUTPUT, ELEMENTS OF A ABOVE THE DIAGONAL
!  ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
!  AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
!  VECTOR D RETURNS THE EIGENVALUES OF A IN ITS FIRST N ELEMENTS.
!  V IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
!  A WHOSE COLUMNS CONTAIN, UPON OUTPUT, THE NORMALIZED
!  EIGENVECTORS OF A.  NROT RETURNS THE NUMBER OF JACOBI ROTATION
!  WHICH WERE REQUIRED.

! THIS SUBROUTINE IS TAKEN FROM 'NUMERICAL RECIPES.'
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: np
DOUBLE PRECISION, INTENT(IN OUT)         :: a(np,np)
INTEGER, INTENT(IN)                      :: n
DOUBLE PRECISION, INTENT(OUT)            :: d(np)
DOUBLE PRECISION, INTENT(OUT)            :: v(np,np)
INTEGER, INTENT(OUT)                     :: nrot

INTEGER :: ip,iq,  i,j
INTEGER, PARAMETER :: nmax=100

DOUBLE PRECISION :: b(nmax),z(nmax), sm,tresh,g,t,h,theta,s,c,tau


! INITIALIZE V TO THE IDENTITY MATRIX
DO i=1,3
  v(i,i)=one
  DO j=1,3
    IF (i /= j)THEN
      v(i,j)=zero
    END IF
  END DO
END DO
! INITIALIZE B AND D TO THE DIAGONAL OF A, AND Z TO ZERO.
!  THE VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
!  IN EQUATION (11.1.14)

DO ip = 1,n
  b(ip) = a(ip,ip)
  d(ip) = b(ip)
  z(ip) = 0.d0
END DO


! BEGIN ITERATION

nrot = 0
DO i=1,50
  
!         SUM OFF-DIAGONAL ELEMENTS
  
  sm = 0.d0
  DO ip=1,n-1
    DO iq=ip+1,n
      sm = sm + DABS(a(ip,iq))
    END DO
  END DO
  
!          IF SM = 0., THEN RETURN.  THIS IS THE NORMAL RETURN,
!          WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE
!          UNDERFLOW.
  
  IF (sm == 0.d0) RETURN
  
!          IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
!           |A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE,
!           SEE EQUATION (11.1.25).  THEREAFTER TRESH = 0.
  
  IF (i < 4) THEN
    tresh = 0.2D0*sm/n**2
  ELSE
    tresh = 0.d0
  END IF
  
  DO ip=1,n-1
    DO iq=ip+1,n
      g = 100.d0*DABS(a(ip,iq))
      
!              AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE
!               OFF-DIAGONAL ELEMENT IS SMALL.
      
      IF ((i > 4).AND.(DABS(d(ip))+g == DABS(d(ip)))  &
            .AND.(DABS(d(iq))+g == DABS(d(iq)))) THEN
        a(ip,iq) = 0.d0
      ELSE IF (DABS(a(ip,iq)) > tresh) THEN
        h = d(iq) - d(ip)
        IF (DABS(h)+g == DABS(h)) THEN
          
!                  T = 1./(2.*THETA), EQUATION (11.1.10)
          
          t =a(ip,iq)/h
        ELSE
          theta = 0.5D0*h/a(ip,iq)
          t =1.d0/(DABS(theta)+DSQRT(1.d0+theta**2.d0))
          IF (theta < 0.d0) t = -t
        END IF
        c = 1.d0/DSQRT(1.d0 + t**2.d0)
        s = t*c
        tau = s/(1.d0 + c)
        h = t*a(ip,iq)
        z(ip) = z(ip) - h
        z(iq) = z(iq) + h
        d(ip) = d(ip) - h
        d(iq) = d(iq) + h
        a(ip,iq) = 0.d0
        
!               CASE OF ROTATIONS 1 <= J < P
        
        DO j=1,ip-1
          g = a(j,ip)
          h = a(j,iq)
          a(j,ip) = g - s*(h + g*tau)
          a(j,iq) = h + s*(g - h*tau)
        END DO
        
!                CASE OF ROTATIONS P < J < Q
        
        DO j=ip+1,iq-1
          g = a(ip,j)
          h = a(j,iq)
          a(ip,j) = g - s*(h + g*tau)
          a(j,iq) = h + s*(g - h*tau)
        END DO
        
!                 CASE OF ROTATIONS Q < J <= N
        
        DO j=iq+1,n
          g = a(ip,j)
          h = a(iq,j)
          a(ip,j) = g - s*(h + g*tau)
          a(iq,j) = h + s*(g - h*tau)
        END DO
        DO j = 1,n
          g = v(j,ip)
          h = v(j,iq)
          v(j,ip) = g - s*(h + g*tau)
          v(j,iq) = h + s*(g - h*tau)
        END DO
        nrot = nrot + 1
      END IF
    END DO
  END DO
  
!          UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z
  
  DO ip=1,n
    b(ip) = b(ip) + z(ip)
    d(ip) = b(ip)
    z(ip) = 0.d0
  END DO
END DO

! IF THE ALGORITHM HAS REACHED THIS STAGE, THEN THERE
!  ARE TOO MANY SWEEPS.  PRINT A DIAGNOSTIC AND CUT THE
!  TIME INCREMENT.

WRITE (*,'(/1X,A/)') '50 ITERATIONS IN JACOBI SHOULD NEVER HAPPEN'

RETURN
END SUBROUTINE jacobi

!**********************************************************************

SUBROUTINE eigsrt(d,v,n,np)

!     GIVEN THE EIGENVALUES D AND EIGENVECTORS V AS OUTPUT FROM
!     JACOBI, THIS SUBROUTINE SORTS THE EIGENVALUES INTO ASCENDING
!     ORDER AND REARRANGES THE COLMNS OF V ACCORDINGLY.

!     THE SUBROUTINE WAS TAKEN FROM 'NUMERICAL RECIPES.'
use global

DOUBLE PRECISION, INTENT(IN OUT)         :: d(np)
DOUBLE PRECISION, INTENT(IN OUT)         :: v(np,np)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: np


INTEGER :: i,j,k

DOUBLE PRECISION :: p

DO i=1,n-1
  k = i
  p = d(i)
  DO j=i+1,n
    IF (d(j) >= p) THEN
      k = j
      p = d(j)
    END IF
  END DO
  IF (k /= i) THEN
    d(k) = d(i)
    d(i) = p
    DO j=1,n
      p = v(j,i)
      v(j,i) = v(j,k)
      v(j,k) = p
    END DO
  END IF
END DO

RETURN
END SUBROUTINE eigsrt
subroutine pk2anisomatfic(afic, daniso, cbar, inv4, st0, ndi)
      !> ANISOTROPIC MATRIX: 2PK 'FICTICIOUS' STRESS TENSOR
      !! INPUT:
      !! DANISO - ANISOTROPIC STRAIN-ENERGY DERIVATIVES
      !! CBAR - DEVIATORIC LEFT CAUCHY-GREEN TENSOR
      !! INV1, INV4 - CBAR INVARIANTS
      !! UNIT2 - 2ND ORDER IDENTITY TENSOR
      !! OUTPUT:
      !! AFIC - 2ND PIOLA KIRCHOOF 'FICTICIOUS' STRESS TENSOR
      use global
      implicit none

      integer, intent(in) :: ndi
      real(8), intent(in) :: daniso(3), cbar(3, 3), inv4
      real(8), intent(in) :: st0(3, 3)
      real(8), intent(out) :: afic(ndi, ndi)
      real(8) :: dudi4, di4dc(3, 3)

      ! FIRST DERIVATIVE OF SSEANISO IN ORDER TO I4
      dudi4 = daniso(1)

      di4dc = st0

      afic = 2.0d0 * (dudi4 * di4dc)

      return
end subroutine pk2anisomatfic
SUBROUTINE invariants(a,inv1,inv2,ndi)



!>    1ST AND 2ND INVARIANTS OF A TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: inv1
DOUBLE PRECISION, INTENT(OUT)            :: inv2



INTEGER :: i1
DOUBLE PRECISION :: aa(ndi,ndi)
DOUBLE PRECISION :: inv1aa

inv1=zero
inv1aa=zero
aa=matmul(a,a)
DO i1=1,ndi
  inv1=inv1+a(i1,i1)
  inv1aa=inv1aa+aa(i1,i1)
END DO
inv2=(one/two)*(inv1*inv1-inv1aa)

RETURN
END SUBROUTINE invariants
SUBROUTINE pk2isomatfic(fic,diso,cbar,cbari1,unit2,ndi)



!>     ISOTROPIC MATRIX: 2PK 'FICTICIOUS' STRESS TENSOR
!      INPUT:
!       DISO - STRAIN-ENERGY DERIVATIVES
!       CBAR - DEVIATORIC LEFT CAUCHY-GREEN TENSOR
!       CBARI1,CBARI2 - CBAR INVARIANTS
!       UNIT2 - 2ND ORDER IDENTITY TENSOR
!      OUTPUT:
!       FIC - 2ND PIOLA KIRCHOOF 'FICTICIOUS' STRESS TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: fic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: diso(5)
DOUBLE PRECISION, INTENT(IN)             :: cbar(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: cbari1
DOUBLE PRECISION, INTENT(IN)             :: unit2(ndi,ndi)



INTEGER :: i1,j1

DOUBLE PRECISION :: dudi1,dudi2
DOUBLE PRECISION :: aux1,aux2

dudi1=diso(1)
dudi2=diso(2)

aux1=two*(dudi1+cbari1*dudi2)
aux2=-two*dudi2

DO i1=1,ndi
  DO j1=1,ndi
    fic(i1,j1)=aux1*unit2(i1,j1)+aux2*cbar(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE pk2isomatfic
SUBROUTINE fslip(f,fbar,det,ndi)



!>     DISTORTION GRADIENT
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: fbar(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: det


INTEGER :: i1,j1

DOUBLE PRECISION :: scale1

!     JACOBIAN DETERMINANT/VOLUME RATIO (J = det(F))
det = f(1,1) * f(2,2) * f(3,3) - f(1,2) * f(2,1) * f(3,3)

IF (ndi == 3) THEN
  det = det + f(1,2) * f(2,3) * f(3,1) + f(1,3) * f(3,2) * f(2,1)  &
      - f(1,3) * f(3,1) * f(2,2) - f(2,3) * f(3,2) * f(1,1)
END IF

scale1=det**(-one /three)

DO i1=1,ndi
  DO j1=1,ndi
    fbar(i1,j1)=scale1*f(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE fslip
SUBROUTINE sdvread(statev)
use global
implicit none
!>    VISCOUS DISSIPATION: READ STATE VARS
DOUBLE PRECISION, INTENT(IN)             :: statev(nsdv)




RETURN

END SUBROUTINE sdvread
SUBROUTINE isomat(sseiso,diso,c10,cbari1)



!>     ISOTROPIC MATRIX : ISOCHORIC SEF AND DERIVATIVES
use global
IMPLICIT NONE


DOUBLE PRECISION, INTENT(OUT)            :: sseiso
DOUBLE PRECISION, INTENT(OUT)            :: diso(5)
DOUBLE PRECISION, INTENT(IN)             :: c10
!DOUBLE PRECISION, INTENT(IN)             :: c01
DOUBLE PRECISION, INTENT(IN OUT)         :: cbari1
!DOUBLE PRECISION, INTENT(IN OUT)         :: cbari2


SSEISO=C10*(CBARI1-THREE)

!FIRST DERIVATIVE OF SSEISO IN ORDER TO I1
DISO(1)=C10
!FIRST DERIVATIVE OF SSEISO IN ORDER TO I2
DISO(2)=ZERO
!SECOND DERIVATIVE OF SSEISO IN ORDER TO I1
DISO(3)=ZERO
!SECOND DERIVATIVE OF SSEISO IN ORDER TO I2
DISO(4)=ZERO
!SECOND DERIVATIVE OF SSEISO IN ORDER TO I1 AND I2
DISO(5)=ZERO

RETURN
END SUBROUTINE isomat
SUBROUTINE setvol(cvol,pv,ppv,unit2,unit4s,ndi)



!>    VOLUMETRIC SPATIAL ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cvol(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pv
DOUBLE PRECISION, INTENT(IN OUT)         :: ppv
DOUBLE PRECISION, INTENT(IN OUT)         :: unit2(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit4s(ndi,ndi,ndi,ndi)


INTEGER :: i1,j1,k1,l1



DO i1 = 1, ndi
  DO j1 = 1, ndi
    DO k1 = 1, ndi
      DO l1 = 1, ndi
        cvol(i1,j1,k1,l1)= ppv*unit2(i1,j1)*unit2(k1,l1)  &
            -two*pv*unit4s(i1,j1,k1,l1)
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE setvol
SUBROUTINE getoutdir(outdir, lenoutdir)



!>     GET CURRENT WORKING DIRECTORY
INCLUDE 'aba_param.inc'


CHARACTER (LEN=256), INTENT(IN OUT)      :: outdir
INTEGER, INTENT(OUT)                     :: lenoutdir



CALL getcwd(outdir)
!        OUTDIR=OUTDIR(1:SCAN(OUTDIR,'\',BACK=.TRUE.)-1)
lenoutdir=len_trim(outdir)

RETURN
END SUBROUTINE getoutdir
subroutine anisomat(sseaniso, daniso, diso, k1, k2, kdisp, i4, i1)
      use global
      implicit none

      ! Arguments
      double precision, intent(out) :: sseaniso
      double precision, intent(out) :: daniso(4)
      double precision, intent(inout) :: diso(5)
      double precision, intent(in) :: k1, k2, kdisp, i4, i1

      ! Local variables
      double precision :: dudi1, d2ud2i1
      double precision :: e1, ee2, ee3, dudi4, d2ud2i4, d2dudi1di4, d2dudi2di4

      dudi1 = diso(1)
      d2ud2i1 = diso(3)

      e1 = i4 * (one - three * kdisp) + i1 * kdisp - one
      sseaniso = (k1 / k2) * (dexp(k1 * e1 * e1) - one)

      if (e1 > zero) then
            ee2 = dexp(k2 * e1 * e1)
            ee3 = (one + two * k2 * e1 * e1)

            dudi1 = dudi1 + k1 * kdisp * e1 * ee2
            d2ud2i1 = d2ud2i1 + k1 * kdisp * kdisp * ee3 * ee2

            dudi4 = k1 * (one - three * kdisp) * e1 * ee2
            d2ud2i4 = k1 * ((one - three * kdisp)**two) * ee3 * ee2
            d2dudi1di4 = k1 * (one - three * kdisp) * kdisp * ee3 * ee2
            d2dudi2di4 = zero
      else
            dudi4 = zero
            d2ud2i4 = zero
            d2dudi1di4 = zero
            d2dudi2di4 = zero
            d2ud2i1 = zero
      end if

      ! First derivative of sseaniso with respect to i1
      daniso(1) = dudi4
      ! First derivative of sseaniso with respect to i2
      daniso(2) = d2ud2i4
      ! Second derivative of sseaniso with respect to i1
      daniso(3) = d2dudi1di4
      ! Second derivative of sseaniso with respect to i2
      daniso(4) = d2dudi2di4

      diso(1) = dudi1
      diso(3) = d2ud2i1

      return
end subroutine anisomat
SUBROUTINE contraction42(s,LT,rt,ndi)



!>       DOUBLE CONTRACTION BETWEEN 4TH ORDER AND 2ND ORDER  TENSOR
!>      INPUT:
!>       LT - left 4TH ORDER TENSOR
!>       RT - right  2ND ODER TENSOR
!>      OUTPUT:
!>       S - DOUBLE CONTRACTED TENSOR (2ND ORDER)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: s(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: LT(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rt(ndi,ndi)


INTEGER :: i1,j1,k1,l1


DOUBLE PRECISION :: aux



DO i1=1,ndi
  DO j1=1,ndi
    aux=zero
    DO k1=1,ndi
      DO l1=1,ndi
        aux=aux+LT(i1,j1,k1,l1)*rt(k1,l1)
      END DO
    END DO
    s(i1,j1)=aux
  END DO
END DO
RETURN
END SUBROUTINE contraction42
SUBROUTINE onem(a,aa,aas,ndi)



!>      THIS SUBROUTINE GIVES:
!>          2ND ORDER IDENTITY TENSORS - A
!>          4TH ORDER IDENTITY TENSOR - AA
!>          4TH ORDER SYMMETRIC IDENTITY TENSOR - AAS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: aa(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: aas(ndi,ndi,ndi,ndi)



INTEGER :: i,j,k,l

a = zero
aa = zero
aas = zero

DO i = 1, ndi
  a(i,i) = one
END DO

DO i=1,ndi
  DO j=1,ndi
    DO k=1,ndi
      DO l=1,ndi
        IF (i == k .and. j == l) then
          aa(i,j,k,l) = one
        END IF
        aas(i,j,k,l) = (one/two)*(a(i,k)*a(j,l)+a(i,l)*a(j,k))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE onem
SUBROUTINE csisomatfic(cisomatfic,cmisomatfic,distgr,det,ndi)



!>    ISOTROPIC MATRIX: SPATIAL 'FICTICIOUS' ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: cisomatfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: cmisomatfic(ndi,ndi,ndi,ndi) !(ndi,ndi) in previous version
DOUBLE PRECISION, INTENT(IN OUT)         :: distgr(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det



call push4(cisomatfic,cmisomatfic,distgr,det,ndi)

RETURN
END SUBROUTINE csisomatfic
subroutine tensorprod2(a, b, c, ndi)
       use global
       implicit none
       integer, intent(in) :: ndi
       real(8), intent(in) :: a(ndi, ndi), b(ndi, ndi)
       real(8), intent(out) :: c(ndi, ndi, ndi, ndi)
       integer :: i, j, k, l

       do i = 1, ndi
              do j = 1, ndi
                     do k = 1, ndi
                            do l = 1, ndi
                                   c(i, j, k, l) = a(i, j) * b(k, l)
                            end do
                     end do
              end do
       end do

end subroutine tensorprod2
SUBROUTINE pull4(mat,spatial,finv,det,ndi)



!>        PULL-BACK TIMES DET OF 4TH ORDER TENSOR

use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: mat(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: spatial(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: finv(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: det



INTEGER :: i1,j1,k1,l1,ii1,jj1,kk1,ll1


DOUBLE PRECISION :: aux


DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        aux=zero
        DO ii1=1,ndi
          DO jj1=1,ndi
            DO kk1=1,ndi
              DO ll1=1,ndi
                aux=aux+det* finv(i1,ii1)*finv(j1,jj1)*  &
                    finv(k1,kk1)*finv(l1,ll1)*spatial(ii1,jj1,kk1,ll1)
              END DO
            END DO
          END DO
        END DO
        mat(i1,j1,k1,l1)=aux
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE pull4
SUBROUTINE metiso(cmiso,cmfic,pl,pkiso,pkfic,c,unit2,det,ndi)



!>    ISOCHORIC MATERIAL ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cmiso(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: cmfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pl(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pkiso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pkfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit2(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: det



INTEGER :: i1,j1,k1,l1
DOUBLE PRECISION :: cisoaux(ndi,ndi,ndi,ndi), cisoaux1(ndi,ndi,ndi,ndi),  &
    plt(ndi,ndi,ndi,ndi),cinv(ndi,ndi), pll(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: trfic,xx,yy,zz, aux,aux1

CALL matinv3d(c,cinv,ndi)
cisoaux1=zero
cisoaux=zero
CALL contraction44(cisoaux1,pl,cmfic,ndi)
DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        plt(i1,j1,k1,l1)=pl(k1,l1,i1,j1)
      END DO
    END DO
  END DO
END DO

CALL contraction44(cisoaux,cisoaux1,plt,ndi)

trfic=zero
aux=det**(-two/three)
aux1=aux**two
CALL contraction22(trfic,aux*pkfic,c,ndi)

DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        xx=aux1*cisoaux(i1,j1,k1,l1)
        pll(i1,j1,k1,l1)=(one/two)*(cinv(i1,k1)*cinv(j1,l1)+  &
            cinv(i1,l1)*cinv(j1,k1))- (one/three)*cinv(i1,j1)*cinv(k1,l1)
        yy=trfic*pll(i1,j1,k1,l1)
        zz=pkiso(i1,j1)*cinv(k1,l1)+cinv(i1,j1)*pkiso(k1,l1)
        
        cmiso(i1,j1,k1,l1)=xx+(two/three)*yy-(two/three)*zz
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE metiso
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
SUBROUTINE matinv3d(a,a_inv,ndi)
!>    INVERSE OF A 3X3 MATRIX
!     RETURN THE INVERSE OF A(3,3) - A_INV
use global

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN)             :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: a_inv(ndi,ndi)

DOUBLE PRECISION :: det_a,det_a_inv

det_a = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) -  &
    a(2,1)*(a(1,2)*a(3,3) - a(3,2)*a(1,3)) +  &
    a(3,1)*(a(1,2)*a(2,3) - a(2,2)*a(1,3))

IF (det_a <= 0.d0) THEN
  WRITE(*,*) 'WARNING: SUBROUTINE MATINV3D:'
  WRITE(*,*) 'WARNING: DET OF MAT=',det_a
  RETURN
END IF

det_a_inv = 1.d0/det_a

a_inv(1,1) = det_a_inv*(a(2,2)*a(3,3)-a(3,2)*a(2,3))
a_inv(1,2) = det_a_inv*(a(3,2)*a(1,3)-a(1,2)*a(3,3))
a_inv(1,3) = det_a_inv*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
a_inv(2,1) = det_a_inv*(a(3,1)*a(2,3)-a(2,1)*a(3,3))
a_inv(2,2) = det_a_inv*(a(1,1)*a(3,3)-a(3,1)*a(1,3))
a_inv(2,3) = det_a_inv*(a(2,1)*a(1,3)-a(1,1)*a(2,3))
a_inv(3,1) = det_a_inv*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
a_inv(3,2) = det_a_inv*(a(3,1)*a(1,2)-a(1,1)*a(3,2))
a_inv(3,3) = det_a_inv*(a(1,1)*a(2,2)-a(2,1)*a(1,2))

RETURN
END SUBROUTINE matinv3d
SUBROUTINE projeul(a,aa,pe,ndi)



!>    EULERIAN PROJECTION TENSOR
!      INPUTS:
!          IDENTITY TENSORS - A, AA
!      OUTPUTS:
!          4TH ORDER SYMMETRIC EULERIAN PROJECTION TENSOR - PE
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: aa(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: pe(ndi,ndi,ndi,ndi)



INTEGER :: i,j,k,l



DO i=1,ndi
  DO j=1,ndi
    DO k=1,ndi
      DO l=1,ndi
        pe(i,j,k,l)=aa(i,j,k,l)-(one/three)*(a(i,j)*a(k,l))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE projeul
SUBROUTINE sigiso(siso,sfic,pe,ndi)



!>    ISOCHORIC CAUCHY STRESS
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: siso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pe(ndi,ndi,ndi,ndi)


CALL contraction42(siso,pe,sfic,ndi)

RETURN
END SUBROUTINE sigiso
SUBROUTINE projlag(c,aa,pl,ndi)



!>    LAGRANGIAN PROJECTION TENSOR
!      INPUTS:
!          IDENTITY TENSORS - A, AA
!          ISOCHORIC LEFT CAUCHY GREEN TENSOR - C
!          INVERSE OF C - CINV
!      OUTPUTS:
!          4TH ORDER SYMMETRIC LAGRANGIAN PROJECTION TENSOR - PL
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: aa(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: pl(ndi,ndi,ndi,ndi)



INTEGER :: i,j,k,l

DOUBLE PRECISION :: cinv(ndi,ndi)

CALL matinv3d(c,cinv,ndi)

DO i=1,ndi
  DO j=1,ndi
    DO k=1,ndi
      DO l=1,ndi
        pl(i,j,k,l)=aa(i,j,k,l)-(one/three)*(cinv(i,j)*c(k,l))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE projlag
SUBROUTINE vol(ssev,pv,ppv,k,det)

! Code converted using TO_F90 by Alan Miller
! Date: 2020-12-12  Time: 12:08:12

!>     VOLUMETRIC CONTRIBUTION :STRAIN ENERGY FUNCTION AND DERIVATIVES
use global
implicit none


DOUBLE PRECISION :: g, aux
DOUBLE PRECISION, INTENT(OUT)            :: ssev
DOUBLE PRECISION, INTENT(OUT)            :: pv
DOUBLE PRECISION, INTENT(OUT)            :: ppv
DOUBLE PRECISION, INTENT(IN)             :: k
DOUBLE PRECISION, INTENT(IN)             :: det


g=(one/four)*(det*det-one-two*LOG(det))

ssev=k*g

pv=k*(one/two)*(det-one/det)
aux=k*(one/two)*(one+one/(det*det))
ppv=pv+det*aux

RETURN
END SUBROUTINE vol
SUBROUTINE stretch(c,b,u,v,ndi)



!>    STRETCH TENSORS

use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi

DOUBLE PRECISION, INTENT(IN OUT)         :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: b(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: u(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: v(ndi,ndi)




DOUBLE PRECISION :: eigval(ndi),omega(ndi),eigvec(ndi,ndi)

CALL spectral(c,omega,eigvec)

eigval(1) = DSQRT(omega(1))
eigval(2) = DSQRT(omega(2))
eigval(3) = DSQRT(omega(3))

u(1,1) = eigval(1)
u(2,2) = eigval(2)
u(3,3) = eigval(3)

u = matmul(matmul(eigvec,u),transpose(eigvec))

CALL spectral(b,omega,eigvec)

eigval(1) = DSQRT(omega(1))
eigval(2) = DSQRT(omega(2))
eigval(3) = DSQRT(omega(3))
!      write(*,*) eigvec(1,1),eigvec(2,1),eigvec(3,1)

v(1,1) = eigval(1)
v(2,2) = eigval(2)
v(3,3) = eigval(3)

v = matmul(matmul(eigvec,v),transpose(eigvec))
RETURN
END SUBROUTINE stretch
subroutine CMATANISOMATFIC(CMANISOMATFIC, M0, DANISO, UNIT2, DET, NDI)
      ! ANISOTROPIC MATRIX: MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
      use global
      implicit none

      integer :: NDI, I, J, K, L
      real(8) :: CMANISOMATFIC(NDI,NDI,NDI,NDI), UNIT2(NDI,NDI)
      real(8) :: M0(NDI,NDI), DANISO(3), DET
      real(8) :: CINV4(NDI,NDI,NDI,NDI), CINV14(NDI,NDI,NDI,NDI)
      real(8) :: D2UDI4, D2UDI1DI4
      real(8) :: IMM(NDI,NDI,NDI,NDI), MMI(NDI,NDI,NDI,NDI)
      real(8) :: MM0(NDI,NDI,NDI,NDI)

      ! 2ND DERIVATIVE OF SSEANISO IN ORDER TO I4
      D2UDI4 = DANISO(2)
      ! 2ND DERIVATIVE OF SSEANISO IN ORDER TO I1 AND I4
      D2UDI1DI4 = DANISO(3)

      call TENSORPROD2(M0, M0, MM0, NDI)
      call TENSORPROD2(UNIT2, M0, IMM, NDI)
      call TENSORPROD2(M0, UNIT2, MMI, NDI)

      do I = 1, NDI
            do J = 1, NDI
                  do K = 1, NDI
                        do L = 1, NDI
                              CINV4(I,J,K,L) = D2UDI4 * MM0(I,J,K,L)
                              CINV14(I,J,K,L) = D2UDI1DI4 * (IMM(I,J,K,L) + MMI(I,J,K,L))
                              CMANISOMATFIC(I,J,K,L) = FOUR * (CINV4(I,J,K,L) + CINV14(I,J,K,L))
                        end do
                  end do
            end do
      end do

      return
end subroutine CMATANISOMATFIC
SUBROUTINE initialize(statev)
use global
IMPLICIT NONE

!      DOUBLE PRECISION TIME(2),KSTEP
INTEGER :: pos1, i
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)


pos1=0
!       DETERMINANT
statev(pos1+1)=one
!       CL RELATIVE STIFFNESS
DO i = pos1+2, nsdv
    statev(i)=zero
END DO
!        CONTRACTION VARIANCE
!statev(pos1+2)=zero

RETURN

END SUBROUTINE initialize
SUBROUTINE sigisomatfic(sfic,pkfic,f,det,ndi)



!>    ISOTROPIC MATRIX:  ISOCHORIC CAUCHY STRESS
use global
IMPLICIT NONE


INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pkfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det





CALL push2(sfic,pkfic,f,det,ndi)

RETURN
END SUBROUTINE sigisomatfic
SUBROUTINE sigvol(svol,pv,unit2,ndi)



!>    VOLUMETRIC CAUCHY STRESS

use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: svol(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pv
DOUBLE PRECISION, INTENT(IN)             :: unit2(ndi,ndi)



INTEGER :: i1,j1



DO i1=1,ndi
  DO j1=1,ndi
    svol(i1,j1)=pv*unit2(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE sigvol
SUBROUTINE deformation(f,c,b,ndi)



!>     RIGHT AND LEFT CAUCHY-GREEN DEFORMATION TENSORS
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: b(ndi,ndi)


!     RIGHT CAUCHY-GREEN DEFORMATION TENSOR
c=matmul(transpose(f),f)
!     LEFT CAUCHY-GREEN DEFORMATION TENSOR
b=matmul(f,transpose(f))
RETURN
END SUBROUTINE deformation
!********************************************************************
! Record of revisions:                                              |
!        Date        Programmer        Description of change        |
!        ====        ==========        =====================        |
!     15/11/2017    Joao Ferreira      cont mech general eqs        |
!     01/11/2018    Joao Ferreira      comments added               |
!--------------------------------------------------------------------
!     Description:
!     UMAT: IMPLEMENTATION OF THE CONSTITUTIVE EQUATIONS BASED UPON 
!     A STRAIN-ENERGY FUNCTION (SEF).
!     THIS CODE, AS IS, EXPECTS A SEF BASED ON THE INVARIANTS OF THE 
!     CAUCHY-GREEN TENSORS. A VISCOELASTIC COMPONENT IS ALSO 
!     INCLUDED IF NEEDED. 
!     YOU CAN CHOOSE TO COMPUTE AT THE MATERIAL FRAME AND THEN 
!     PUSHFORWARD OR  COMPUTE AND THE SPATIAL FRAME DIRECTLY.
!--------------------------------------------------------------------
!     IF YOU WANT TO ADAPT THE CODE ACCORDING TO YOUR SEF:
!    ISOMAT - DERIVATIVES OF THE SEF IN ORDER TO THE INVARIANTS
!    ADD OTHER CONTRIBUTIONS: STRESS AND TANGENT MATRIX
!-------------------------------------------------------------------- 
!      STATE VARIABLES: CHECK ROUTINES - INITIALIZE, WRITESDV, READSDV.
!--------------------------------------------------------------------              
!     UEXTERNALDB: READ FILAMENTS ORIENTATION AND PREFERED DIRECTION
!--------------------------------------------------------------------
!---------------------------------------------------------------------
      SUBROUTINE MATERIAL(SIGMA,STATEV,DDSIGDDE,DFGRD0,DFGRD1,DET, &
      TIME,DTIME,PREDEF,NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS, &
      PNEWDT,NOEL,NPT,KSTEP,KINC)
!----------------------------------------------------------------------
!--------------------------- DECLARATIONS -----------------------------
!----------------------------------------------------------------------
      use global
      IMPLICIT NONE

!     ADD COMMON BLOCKS HERE IF NEEDED (and in uexternal)
!      COMMON /KBLOCK/KBLOCK

      CHARACTER(len=8) :: CMNAME

      INTEGER :: NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT, &
                 LAYER, KSPT, KSTEP, KINC

      REAL(KIND=8) :: STRESS(NTENS), STATEV(NSTATEV), &
                      DDSDDE(NTENS,NTENS), DDSDDT(NTENS), DRPLDE(NTENS), &
                      STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1), DPRED(1), &
                      PROPS(NPROPS), COORDS(3), DROT(3,3), DFGRD0(3,3), DFGRD1(3,3), &
                      FIBORI(NELEM,4)

      REAL(KIND=8) :: SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP, &
                      DTEMP, PNEWDT, CELENT

      INTEGER :: NTERM

!     FLAGS
!      INTEGER :: FLAG1

!     UTILITY TENSORS
      REAL(KIND=8) :: UNIT2(NDI,NDI), UNIT4(NDI,NDI,NDI,NDI), &
                      UNIT4S(NDI,NDI,NDI,NDI), &
                      PROJE(NDI,NDI,NDI,NDI), PROJL(NDI,NDI,NDI,NDI)

!     KINEMATICS
      REAL(KIND=8) :: DISTGR(NDI,NDI), C(NDI,NDI), B(NDI,NDI), &
                      CBAR(NDI,NDI), BBAR(NDI,NDI), DISTGRINV(NDI,NDI), &
                      UBAR(NDI,NDI), VBAR(NDI,NDI), ROT(NDI,NDI), &
                      DFGRD1INV(NDI,NDI)
      REAL(KIND=8) :: DET, CBARI1, CBARI2

!     VOLUMETRIC CONTRIBUTION
      REAL(KIND=8) :: PKVOL(NDI,NDI), SVOL(NDI,NDI), &
                      CVOL(NDI,NDI,NDI,NDI), CMVOL(NDI,NDI,NDI,NDI)
      REAL(KIND=8) :: KBULK, PV, PPV, SSEV

!     ISOCHORIC CONTRIBUTION
      REAL(KIND=8) :: SISO(NDI,NDI), PKISO(NDI,NDI), PK2(NDI,NDI), &
                      CISO(NDI,NDI,NDI,NDI), CMISO(NDI,NDI,NDI,NDI), &
                      SFIC(NDI,NDI), CFIC(NDI,NDI,NDI,NDI), &
                      PKFIC(NDI,NDI), CMFIC(NDI,NDI,NDI,NDI)

!     ISOCHORIC ISOTROPIC CONTRIBUTION
      REAL(KIND=8) :: C10, C01, SSEISO, DISO(5), PKMATFIC(NDI,NDI), &
                      SMATFIC(NDI,NDI), SISOMATFIC(NDI,NDI), &
                      CMISOMATFIC(NDI,NDI,NDI,NDI), &
                      CISOMATFIC(NDI,NDI,NDI,NDI)
      REAL(KIND=8) :: VORIF(3), VD(3), M0(3,3), MM(3,3), &
                      VORIF2(3), VD2(3), N0(3,3), NN(3,3)

!     LIST VARS OF OTHER CONTRIBUTIONS HERE

!     JAUMMAN RATE CONTRIBUTION (REQUIRED FOR ABAQUS UMAT)
      REAL(KIND=8) :: CJR(NDI,NDI,NDI,NDI)

!     CAUCHY STRESS AND ELASTICITY TENSOR
      REAL(KIND=8) :: SIGMA(NDI,NDI), DDSIGDDE(NDI,NDI,NDI,NDI), &
                      DDPKDDE(NDI,NDI,NDI,NDI)

!     TESTING/DEBUG VARS
      REAL(KIND=8) :: STEST(NDI,NDI), CTEST(NDI,NDI,NDI,NDI)
      INTEGER :: I1, J1, K1, L1
!----------------------------------------------------------------------
!-------------------------- INITIALIZATIONS ---------------------------
!----------------------------------------------------------------------
!     IDENTITY AND PROJECTION TENSORS
      UNIT2 = 0.0
      UNIT4 = 0.0
      UNIT4S = 0.0
      PROJE = 0.0
      PROJL = 0.0

!     KINEMATICS
      DISTGR = 0.0
      C = 0.0
      B = 0.0
      CBAR = 0.0
      BBAR = 0.0
      UBAR = 0.0
      VBAR = 0.0
      ROT = 0.0
      DET = 0.0
      CBARI1 = 0.0
      CBARI2 = 0.0

!     VOLUMETRIC
      PKVOL = 0.0
      SVOL = 0.0
      CVOL = 0.0
      KBULK = 0.0
      PV = 0.0
      PPV = 0.0
      SSEV = 0.0

!     ISOCHORIC
      SISO = 0.0
      PKISO = 0.0
      PK2 = 0.0
      CISO = 0.0
      CFIC = 0.0
      SFIC = 0.0
      PKFIC = 0.0

!     ISOTROPIC
      C10 = 0.0
      C01 = 0.0
      SSEISO = 0.0
      DISO = 0.0
      PKMATFIC = 0.0
      SMATFIC = 0.0
      SISOMATFIC = 0.0
      CMISOMATFIC = 0.0
      CISOMATFIC = 0.0

!     INITIALIZE OTHER CONT HERE

!     JAUMANN RATE
      CJR = 0.0

!     TOTAL CAUCHY STRESS AND ELASTICITY TENSORS
      SIGMA = 0.0
      DDSIGDDE = 0.0
!----------------------------------------------------------------------
!------------------------ IDENTITY TENSORS ----------------------------
!----------------------------------------------------------------------
      CALL ONEM(UNIT2, UNIT4, UNIT4S, NDI)
!----------------------------------------------------------------------
!------------------- MATERIAL CONSTANTS AND DATA ----------------------
!----------------------------------------------------------------------
!     VOLUMETRIC
      KBULK = PROPS(1)
!     ISOCHORIC ISOTROPIC NEO HOOKE
      C10 = PROPS(2)
!     NUMERICAL COMPUTATIONS
      NTERM = 60
!
!     STATE VARIABLES
!
      IF ((TIME(1) == 0.0) .AND. (KSTEP == 1)) THEN
         CALL INITIALIZE(STATEV)
      END IF
!        READ STATEV
      CALL SDVREAD(STATEV)
!      
!----------------------------------------------------------------------
!---------------------------- KINEMATICS ------------------------------
!----------------------------------------------------------------------
!     DISTORTION GRADIENT
      CALL FSLIP(DFGRD1, DISTGR, DET, NDI)
!     INVERSE OF DISTORTION GRADIENT
      CALL MATINV3D(DFGRD1, DFGRD1INV, NDI)
!     INVERSE OF DISTORTION GRADIENT
      CALL MATINV3D(DISTGR, DISTGRINV, NDI)
!     CAUCHY-GREEN DEFORMATION TENSORS
      CALL DEFORMATION(DFGRD1, C, B, NDI)
      CALL DEFORMATION(DISTGR, CBAR, BBAR, NDI)      
!     INVARIANTS OF DEVIATORIC DEFORMATION TENSORS
      CALL INVARIANTS(CBAR, CBARI1, CBARI2, NDI)
!     STRETCH TENSORS
      CALL STRETCH(CBAR, BBAR, UBAR, VBAR, NDI)
!     ROTATION TENSORS
      CALL ROTATION(DISTGR, ROT, UBAR, NDI)
!     DEVIATORIC PROJECTION TENSORS
      CALL PROJEUL(UNIT2, UNIT4S, PROJE, NDI)
!
      CALL PROJLAG(C, UNIT4, PROJL, NDI)
!----------------------------------------------------------------------
!--------------------- CONSTITUTIVE RELATIONS  ------------------------
!----------------------------------------------------------------------
!
!---- VOLUMETRIC ------------------------------------------------------
!     STRAIN-ENERGY AND DERIVATIVES (CHANGE HERE ACCORDING TO YOUR MODEL)
      CALL VOL(SSEV, PV, PPV, KBULK, DET)
      CALL ISOMAT(SSEISO, DISO, C10, CBARI1)
!
!---- ISOCHORIC ISOTROPIC ---------------------------------------------
!     PK2 'FICTICIOUS' STRESS TENSOR
      CALL PK2ISOMATFIC(PKMATFIC, DISO, CBAR, CBARI1, UNIT2, NDI)
!     CAUCHY 'FICTICIOUS' STRESS TENSOR
      CALL SIGISOMATFIC(SISOMATFIC, PKMATFIC, DISTGR, DET, NDI)
!     'FICTICIOUS' MATERIAL ELASTICITY TENSOR
      CALL CMATISOMATFIC(CMISOMATFIC, CBAR, CBARI1, CBARI2, DISO, UNIT2, UNIT4, DET, NDI)
!     'FICTICIOUS' SPATIAL ELASTICITY TENSOR
      CALL CSISOMATFIC(CISOMATFIC, CMISOMATFIC, DISTGR, DET, NDI)
!
!----------------------------------------------------------------------
!     SUM OF ALL ELASTIC CONTRIBUTIONS
!----------------------------------------------------------------------
!     STRAIN-ENERGY
      SSE = SSEV + SSEISO
!     PK2 'FICTICIOUS' STRESS
      PKFIC = PKMATFIC
!     CAUCHY 'FICTICIOUS' STRESS
      SFIC = SISOMATFIC
!     MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
      CMFIC = CMISOMATFIC
!     SPATIAL 'FICTICIOUS' ELASTICITY TENSOR
      CFIC = CISOMATFIC
!
!----------------------------------------------------------------------
!-------------------------- STRESS MEASURES ---------------------------
!----------------------------------------------------------------------
!
!---- VOLUMETRIC ------------------------------------------------------
!      PK2 STRESS
      CALL PK2VOL(PKVOL, PV, C, NDI)
!      CAUCHY STRESS
      CALL SIGVOL(SVOL, PV, UNIT2, NDI)
!
!---- ISOCHORIC -------------------------------------------------------
!      PK2 STRESS
      CALL PK2ISO(PKISO, PKFIC, PROJL, DET, NDI)
!      CAUCHY STRESS
      CALL SIGISO(SISO, SFIC, PROJE, NDI)
!
!---- VOLUMETRIC + ISOCHORIC ------------------------------------------
!      PK2 STRESS
      PK2 = PKVOL + PKISO
!      CAUCHY STRESS
      SIGMA = SVOL + SISO
!----------------------------------------------------------------------
!-------------------- MATERIAL ELASTICITY TENSOR ----------------------
!----------------------------------------------------------------------
!
!---- VOLUMETRIC ------------------------------------------------------
!
      CALL METVOL(CMVOL, C, PV, PPV, DET, NDI)
!
!---- ISOCHORIC -------------------------------------------------------
!
      CALL METISO(CMISO, CMFIC, PROJL, PKISO, PKFIC, C, UNIT2, DET, NDI)
!
!----------------------------------------------------------------------
!
      DDPKDDE = CMVOL + CMISO
!
!----------------------------------------------------------------------
!--------------------- SPATIAL ELASTICITY TENSOR ----------------------
!----------------------------------------------------------------------
!
!---- VOLUMETRIC ------------------------------------------------------
!
      CALL SETVOL(CVOL, PV, PPV, UNIT2, UNIT4S, NDI)
!
!---- ISOCHORIC -------------------------------------------------------
!
      CALL SETISO(CISO, CFIC, PROJE, SISO, SFIC, UNIT2, NDI)
!
!-----JAUMMAN RATE ----------------------------------------------------
!
      CALL SETJR(CJR, SIGMA, UNIT2, NDI)
!
!----------------------------------------------------------------------
!
!     ELASTICITY TENSOR
      DDSIGDDE = CVOL + CISO + CJR
      !DDSIGDDE = CVOL + CISO
!
!----------------------------------------------------------------------
!------------------------- INDEX ALLOCATION ---------------------------
!----------------------------------------------------------------------
!     VOIGT NOTATION  - FULLY SYMMETRY IMPOSED
      CALL INDEXX(STRESS, DDSDDE, SIGMA, DDSIGDDE, NTENS, NDI)
!
!----------------------------------------------------------------------
!--------------------------- STATE VARIABLES --------------------------
!----------------------------------------------------------------------
!     DO K1 = 1, NTENS
!      STATEV(1:27) = VISCOUS TENSORS
       CALL SDVWRITE(DET, STATEV)
!     END DO
!----------------------------------------------------------------------
      !write(*,*) 'F0'
      !write(*,*) DFGRD0
      !write(*,*) 'F1'
      !write(*,*) DFGRD1(1,1), DFGRD1(1,2), DFGRD1(1,3)
      !write(*,*) DFGRD1(2,1), DFGRD1(2,2), DFGRD1(2,3)
      !write(*,*) DFGRD1(3,1), DFGRD1(3,2), DFGRD1(3,3)
      !write(*,*) 'Cauchy stress'
      !write(*,*) STRESS(1), STRESS(2), STRESS(3), STRESS(4), STRESS(5), STRESS(6)
      RETURN
      END
!----------------------------------------------------------------------
!--------------------------- END OF UMAT ------------------------------
!----------------------------------------------------------------------
!
SUBROUTINE contraction44(s,LT,rt,ndi)



!>       DOUBLE CONTRACTION BETWEEN 4TH ORDER TENSORS
!>      INPUT:
!>       LT - RIGHT 4TH ORDER TENSOR
!>       RT - LEFT  4TH ORDER TENSOR
!>      OUTPUT:
!>       S - DOUBLE CONTRACTED TENSOR (4TH ORDER)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: s(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: LT(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rt(ndi,ndi,ndi,ndi)



INTEGER :: i1,j1,k1,l1,m1,n1


DOUBLE PRECISION :: aux



DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        aux=zero
        DO m1=1,ndi
          DO n1=1,ndi
            aux=aux+LT(i1,j1,m1,n1)*rt(m1,n1,k1,l1)
          END DO
        END DO
        s(i1,j1,k1,l1)=aux
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE contraction44
SUBROUTINE cmatisomatfic(cmisomatfic,cbar,cbari1,cbari2,  &
        diso,unit2,unit4,det,ndi)



!>    ISOTROPIC MATRIX: MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: cmisomatfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: cbar(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: cbari1
DOUBLE PRECISION, INTENT(IN OUT)         :: cbari2
DOUBLE PRECISION, INTENT(IN)             :: diso(5)
DOUBLE PRECISION, INTENT(IN)             :: unit2(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit4(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1
    

DOUBLE PRECISION :: dudi1,dudi2,d2ud2i1,d2ud2i2,d2udi1i2
DOUBLE PRECISION :: aux,aux1,aux2,aux3,aux4
DOUBLE PRECISION :: uij,ukl,cij,ckl

dudi1=diso(1)
dudi2=diso(2)
d2ud2i1=diso(3)
d2ud2i2=diso(4)
d2udi1i2=diso(5)

aux1=four*(d2ud2i1+two*cbari1*d2udi1i2+ dudi2+cbari1*cbari1*d2ud2i2)
aux2=-four*(d2udi1i2+cbari1*d2ud2i2)
aux3=four*d2ud2i2
aux4=-four*dudi2

DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        uij=unit2(i1,j1)
        ukl=unit2(k1,l1)
        cij=cbar(i1,j1)
        ckl=cbar(k1,l1)
        aux=aux1*uij*ukl+ aux2*(uij*ckl+cij*ukl)+aux3*cij*ckl+  &
            aux4*unit4(i1,j1,k1,l1)
        cmisomatfic(i1,j1,k1,l1)=aux * det**(-four/three)
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE cmatisomatfic
subroutine resetdfgrd(dfgrd, ndi)
      use global
      implicit none

      integer, intent(in) :: ndi
      real(8), intent(out) :: dfgrd(ndi, ndi)

      dfgrd = 0.0
      dfgrd(1,1) = 1.0
      dfgrd(2,2) = 1.0
      dfgrd(3,3) = 1.0

end subroutine resetdfgrd
SUBROUTINE metvol(cvol,c,pv,ppv,det,ndi)



!>    VOLUMETRIC MATERIAL ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cvol(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pv
DOUBLE PRECISION, INTENT(IN OUT)         :: ppv
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1
DOUBLE PRECISION :: cinv(ndi,ndi)


CALL matinv3d(c,cinv,ndi)

DO i1 = 1, ndi
  DO j1 = 1, ndi
    DO k1 = 1, ndi
      DO l1 = 1, ndi
        cvol(i1,j1,k1,l1)= det*ppv*cinv(i1,j1)*cinv(k1,l1)  &
            -det*pv*(cinv(i1,k1)*cinv(j1,l1) +cinv(i1,l1)*cinv(j1,k1))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE metvol
SUBROUTINE pk2vol(pkvol,pv,c,ndi)
!SUBROUTINE pk2vol(pkvol,pv,c,ndi, det)

! Shouldn't det be included??

!>    VOLUMETRIC PK2 STRESS
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pkvol(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pv
DOUBLE PRECISION, INTENT(IN OUT)         :: c(ndi,ndi)
!DOUBLE PRECISION, INTENT(IN)             :: det

INTEGER :: i1,j1
DOUBLE PRECISION :: cinv(ndi,ndi)


CALL matinv3d(c,cinv,ndi)

DO i1=1,ndi
  DO j1=1,ndi
    !pkvol(i1,j1)=det*pv*cinv(i1,j1)
    pkvol(i1,j1)=pv*cinv(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE pk2vol
SUBROUTINE rotation(f,r,u,ndi)



!>    COMPUTES ROTATION TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: r(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: u(ndi,ndi)




DOUBLE PRECISION :: uinv(ndi,ndi)

CALL matinv3d(u,uinv,ndi)

r = matmul(f,uinv)
RETURN
END SUBROUTINE rotation
SUBROUTINE pk2iso(pkiso,pkfic,pl,det,ndi)



!>    ISOCHORIC PK2 STRESS TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pkiso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pkfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pl(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: det



INTEGER :: i1,j1

DOUBLE PRECISION :: scale2

CALL contraction42(pkiso,pl,pkfic,ndi)

scale2=det**(-two/three)
DO i1=1,ndi
  DO j1=1,ndi
    pkiso(i1,j1)=scale2*pkiso(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE pk2iso
SUBROUTINE setiso(ciso,cfic,pe,siso,sfic,unit2,ndi)


use global
IMPLICIT NONE

!>    ISOCHORIC SPATIAL ELASTICITY TENSOR

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: ciso(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pe(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: siso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: unit2(ndi,ndi)



INTEGER :: i1,j1,k1,l1
DOUBLE PRECISION :: cisoaux(ndi,ndi,ndi,ndi), cisoaux1(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: trfic,xx,yy,zz

cisoaux1=zero
cisoaux=zero

CALL contraction44(cisoaux1,pe,cfic,ndi)
CALL contraction44(cisoaux,cisoaux1,pe,ndi)

trfic=zero
DO i1=1,ndi
  trfic=trfic+sfic(i1,i1)
END DO

DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        xx=cisoaux(i1,j1,k1,l1)
        yy=trfic*pe(i1,j1,k1,l1)
        zz=siso(i1,j1)*unit2(k1,l1)+unit2(i1,j1)*siso(k1,l1)
        
        ciso(i1,j1,k1,l1)=xx+(two/three)*yy-(two/three)*zz
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE setiso
SUBROUTINE push4(spatial,mat,f,det,ndi)



!>        PIOLA TRANSFORMATION
!>      INPUT:
!>       MAT - MATERIAL ELASTICITY TENSOR
!>       F - DEFORMATION GRADIENT
!>       DET - DEFORMATION DETERMINANT
!>      OUTPUT:
!>       SPATIAL - SPATIAL ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: spatial(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: mat(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det


INTEGER :: i1,j1,k1,l1,ii1,jj1,kk1,ll1


DOUBLE PRECISION :: aux


DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        aux=zero
        DO ii1=1,ndi
          DO jj1=1,ndi
            DO kk1=1,ndi
              DO ll1=1,ndi
                aux=aux+(det**(-one))* f(i1,ii1)*f(j1,jj1)*  &
                    f(k1,kk1)*f(l1,ll1)*mat(ii1,jj1,kk1,ll1)
              END DO
            END DO
          END DO
        END DO
        spatial(i1,j1,k1,l1)=aux
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE push4
SUBROUTINE contraction22(aux,lt,rt,ndi)
!>       DOUBLE CONTRACTION BETWEEN 2nd ORDER AND 2ND ORDER  TENSOR
!>      INPUT:
!>       LT - RIGHT 2ND ORDER TENSOR
!>       RT - LEFT  2nd ODER TENSOR
!>      OUTPUT:
!>       aux - DOUBLE CONTRACTED TENSOR (scalar)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: lt(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rt(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: aux
INTEGER :: i1,j1


    aux=zero
    DO i1=1,ndi
      DO j1=1,ndi
        aux=aux+lt(i1,j1)*rt(j1,i1)
      END DO
    END DO
RETURN
END SUBROUTINE contraction22
SUBROUTINE push2(sig,pk,f,det,ndi)



!>        PIOLA TRANSFORMATION
!>      INPUT:
!>       PK - 2ND PIOLA KIRCHOOF STRESS TENSOR
!>       F - DEFORMATION GRADIENT
!>       DET - DEFORMATION DETERMINANT
!>      OUTPUT:
!>       SIG - CAUCHY STRESS TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sig(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pk(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det


INTEGER :: i1,j1,ii1,jj1


DOUBLE PRECISION :: aux

DO i1=1,ndi
  DO j1=1,ndi
    aux=zero
    DO ii1=1,ndi
      DO jj1=1,ndi
        aux=aux+(det**(-one))*f(i1,ii1)*f(j1,jj1)*pk(ii1,jj1)
      END DO
    END DO
    sig(i1,j1)=aux
  END DO
END DO

RETURN
END SUBROUTINE push2
SUBROUTINE contraction24(s,LT,rt,ndi)



!>       DOUBLE CONTRACTION BETWEEN 4TH ORDER AND 2ND ORDER  TENSOR
!>      INPUT:
!>       LT - RIGHT 2ND ORDER TENSOR
!>       RT - LEFT  4TH ODER TENSOR
!>      OUTPUT:
!>       S - DOUBLE CONTRACTED TENSOR (2ND ORDER)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: s(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: lt(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rt(ndi,ndi,ndi,ndi)



INTEGER :: i1,j1,k1,l1


DOUBLE PRECISION :: aux



DO k1=1,ndi
  DO l1=1,ndi
    aux=zero
    DO i1=1,ndi
      DO j1=1,ndi
        aux=aux+lt(k1,l1)*rt(i1,j1,k1,l1)
      END DO
    END DO
    s(k1,l1)=aux
  END DO
END DO
RETURN
END SUBROUTINE contraction24
