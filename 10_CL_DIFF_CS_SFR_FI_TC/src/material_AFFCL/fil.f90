SUBROUTINE fil(f,ff,dw,ddw, &
            lambdai,lambdaf,lambda0,lambda0f,&
            ll,r0,r0f,mu0,beta,b0,etac,&
            cb,DfDcb)



!>    SINGLE FILAMENT: STRAIN ENERGY DERIVATIVES
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: f
DOUBLE PRECISION, INTENT(OUT)            :: ff
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
DOUBLE PRECISION, INTENT(IN OUT)         :: mu0
DOUBLE PRECISION, INTENT(IN OUT)         :: beta
DOUBLE PRECISION, INTENT(IN OUT)         :: b0
DOUBLE PRECISION, INTENT(IN OUT)         :: etac
DOUBLE PRECISION, INTENT(IN OUT)         :: cb
! DOUBLE PRECISION :: aratio, r0c
DOUBLE PRECISION :: a,b,machep,t
DOUBLE PRECISION :: aux, pi,alpha
DOUBLE PRECISION :: aux0,aux1,aux2,aux3,aux4,aux5,aux6,y
DOUBLE PRECISION :: tt,dGdf,dGdcb,dlfdcb
! DOUBLE PRECISION :: aux00,aux01,aux02,aux03,aux04,aux05

a=zero
b=1.0E09
machep=2.2204E-16
t=1.0E-14 ! 1.0E-6
f=zero

CALL pullforce(f, a, b, machep, t, lambdaf,lambda0f,ll,r0f,mu0,beta,b0)
! write(*,*) 'lambdaf', lambdaf
! write(*,*) 'f', f
pi=four*ATAN(one)
! ff=f*ll*(pi*pi*b0)**(-one)
ff=f*ll*ll*(pi*pi*b0)**(-one)
! ff = 100000000000000000.0

alpha=pi*pi*b0*(ll*ll*mu0)**(-one)

aux0=beta/alpha
aux=alpha*ff
aux1=one+ff+aux*ff
aux2=one+two*aux
aux3=one+aux
! aux4=lambda0*r0*r0*mu0*(ll**(-one))
! OLD
! aux4=lambda0*lambda0*r0*r0*mu0*(ll**(-one))
! aux4=lambda0*lambda0f*r0*r0*mu0*(ll**(-one))
! NEW
aux4=etac*lambda0f*lambda0*r0*r0*mu0*(ll**(-one))
aux5=((one+aux)*(aux1**(-one)))**beta
aux6=one-r0f*((ll)**(-one))

y=aux0*(aux2*aux2*(aux1**(-one)))-beta*(aux2*(aux3**(-one)))-two

! Strain energy derivatives
dw=lambda0*(r0)*f
! dw = pi*pi*r0*b0/(ll*ll)*(((ll/r0-1)/(ll/r0-lambda))**TWO - one)
ddw=aux4*((one+y*aux5*aux6)**(-one))

! Force derivative wrt cb: implicit function theorem on G(f, cb) = 0 (Eq. 80), DfDcb = -dG/dcb / dG/df
!   G = lambdaf*lambda0f*r0f/ll - (1 + f/mu0 - aux6*tt),   tt = aux2*aux5 = (1+2f/mu0)(1+f/mu0)^beta (pi^2/D)^beta
!   D = pi^2*aux1,  L = a*r0f  =>  r0f/ll and aux6 do not depend on cb,  dL/dcb = -2/5 L/cb
tt = aux2*aux5
dGdf = -one/mu0 + aux6*tt*(beta/(mu0 + f) + two/(mu0 + two*f) - beta*aux2*ll*ll/(pi*pi*b0*aux1))
dlfdcb = two/5.d0*etac*(r0 - r0f)/r0f/cb
dGdcb = r0f/ll*dlfdcb*(lambda0f*(lambdai - one) + lambdaf*(lambda0 - one)) &
      + four/5.d0*aux6*beta*tt*(aux1 - one)/(aux1*cb)
DfDcb = -dGdcb/dGdf

RETURN
END SUBROUTINE fil
