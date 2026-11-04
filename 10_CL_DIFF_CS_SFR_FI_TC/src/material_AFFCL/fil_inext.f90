SUBROUTINE fil_inext(f,dw,ddw, &
            lambdai,lambdaf,lambda0,lambda0f,&
            ll,r0,r0f,beta,b0,etac,&
            cb,DfDcb)



!>    SINGLE FILAMENT: STRAIN ENERGY DERIVATIVES
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: f
DOUBLE PRECISION, INTENT(OUT)            :: dw
DOUBLE PRECISION, INTENT(OUT)            :: ddw
DOUBLE PRECISION, INTENT(OUT)            :: DfDcb
DOUBLE PRECISION, INTENT(IN OUT)         :: lambdai
DOUBLE PRECISION, INTENT(IN OUT)         :: lambdaf
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda0
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda0f
DOUBLE PRECISION, INTENT(IN OUT)         :: ll
DOUBLE PRECISION, INTENT(IN OUT)         :: r0
DOUBLE PRECISION, INTENT(IN OUT)         :: r0f
DOUBLE PRECISION, INTENT(IN OUT)         :: beta
DOUBLE PRECISION, INTENT(IN OUT)         :: b0
DOUBLE PRECISION, INTENT(IN OUT)         :: etac
DOUBLE PRECISION, INTENT(IN OUT)         :: cb
DOUBLE PRECISION :: aratio, r0c
DOUBLE PRECISION :: pi, num_ddw, den_ddw
DOUBLE PRECISION :: aux00,aux01,aux02,aux03,aux04,aux05

f=zero
pi=four*ATAN(one)

aratio=ll/r0f
f = pi*pi*b0/(ll*ll) * (((aratio-1)/(aratio-lambdaf))**(1/beta) - 1)

! Strain energy derivatives
dw=lambda0*(r0)*f
num_ddw = pi*pi*b0*etac*lambda0*r0*r0*((aratio-1)/(aratio-lambdaf))**(1/beta)
den_ddw = ll*ll*beta*r0f*(aratio-lambdaf)
ddw = num_ddw / den_ddw

! Force derivative wrt cb
r0c = r0 - r0f
aux00 = two / 5.d0 * cb**(-one)
aux01 = 2 * f
aux02 = etac * r0c / r0f * (lambdai - 1)
aux03 = b0 * pi * pi / (r0f * aratio)**2
aux04 = (aratio - lambdaf) * beta
aux05 = (f + aux03) / aux04
DfDcb = aux00 * (aux01 + aux02 * aux05)

RETURN
END SUBROUTINE fil_inext
