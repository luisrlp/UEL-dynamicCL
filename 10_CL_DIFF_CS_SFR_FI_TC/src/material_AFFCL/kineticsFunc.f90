subroutine kineticsFunc(cbtau, f, df, dHdlambda, dHdth, dHds, args, nargs)

! Outputs: f = H (Backward-Euler residual), df = dH/dcb, dHdlambda = dH/dlambdar (at fixed cb), with
! lambdar = args(1) the stretch seen by the kinetics (lambda_i/lambda_n with stretch-free rebinding);
! the caller converts it to dcb/dlambda_i with dlrdl (consistent tangent)

    use global
    implicit none

    integer, intent(in)              :: nargs
    DOUBLE PRECISION, intent(in out) :: cbtau
    DOUBLE PRECISION, intent(out)    :: f, df, dHdlambda, dHdth, dHds
    DOUBLE PRECISION, intent(in)     :: args(nargs)                                                              
                                                                                                                
    DOUBLE PRECISION                 :: r0f, etac, r0, r0c, fi, ffi, dwi, ddwi, l, mu0str, beta, b0 
    DOUBLE PRECISION                 :: cfmax, cbmax, dx_kT, dt, kon, koff0, koff, thetab, Ri
    DOUBLE PRECISION                 :: kcatch0, dxc_kT, dkoffdf
    DOUBLE PRECISION                 :: lambdar, lambdaif, lambda0, lambda0f, thetaf, cbt, det                                                          
    DOUBLE PRECISION                 :: DfDcb,DRiDcb, aratio
    DOUBLE PRECISION                 :: s, dbdth, bterm
    
    Ri = zero
                                                                                                                
    ! Obtain relevant quantities
    lambdar = args(1)
    lambda0 = args(2)
    aratio  = args(3)
    etac    = args(4)
    mu0str  = args(5)
    beta    = args(6)
    b0      = args(7)
    r0c     = args(8)                                                                                         
    cbmax   = args(9)
    cfmax   = args(10)
    dx_kT   = args(11)
    dt      = args(12)
    kon     = args(13)
    koff0   = args(14)
    thetaf  = args(15)
    cbt     = args(16)
    det     = args(17)
    kcatch0 = args(18)
    dxc_kT  = args(19)
    s = args(20)
    dbdth = args(21)

    r0f = 1.6d0 * (cbtau*1.d3)**(- two / 5.d0)
    l = aratio * r0f
    r0 = r0f + r0c
    
    IF((etac > zero).AND.(etac .LE. one))THEN
        lambdaif=etac*(r0/r0f)*(lambdar-one)+one
        lambda0f=etac*(r0/r0f)*(lambda0-one)+one
    ELSE
        lambdaif=lambdar 
        lambda0f=lambda0
    END IF

    ! Slack filament (total filament stretch, prestretch included, <= 1): no force.
    ! With LAMBDA0 = 1 this is lambdar <= 1; with LAMBDA0 > 1 the filament is loaded at rest
    IF (lambdaif*lambda0f .LE. one) then
        fi = 0.0
        DfDcb = 0.0
        ddwi = zero
    ELSE
        CALL fil(fi,ffi,dwi,ddwi,&
                lambdar,lambdaif,lambda0,lambda0f,&
                l,r0,r0f,mu0str,beta,b0,etac,&
                cbtau,DfDcb)
        !CALL fil_inext(fi,dwi,ddwi,lambdar,lambdaif,lambda0,lambda0f,l,r0,r0f,beta,b0,etac,cbtau,DfDcb)
            ! CALL filpce(lambdar, fi, dwi, ddwi)
    END IF

    ! Unbinding rate (catch-slip) and its force derivative
    CALL koffcs(koff, dkoffdf, fi, koff0, kcatch0, dx_kT, dxc_kT)
    thetab = cbtau / cbmax

    ! Reaction rate and residual                                                                         
    bterm = kon * cfmax * thetaf / (one - thetaf) * s
    Ri = koff * (bterm - cbmax * thetab / (one - thetab))
    f = cbtau - cbt - Ri * dt
    ! Check if the residual is NaN or Inf. Large finite values are legitimate with the
    ! extensible filament and the catch-slip law (koff is capped by koffcs, not bounded by locking)
    if (.not. (abs(f) < huge(f))) then
        write(*,*) 'Error: kinetics residual is NaN or Inf in kineticsFunc'
        write(*,*) 'cbtau =', cbtau
        write(*,*) 'cbt =', cbt
        write(*,*) 'Ri =', Ri
        write(*,*) 'koff =', koff
        write(*,*) 'fi=', fi
        write(*,*) 'lambdaif =', lambdaif
        write(*,*) 'r0f =', r0f
        write(*,*) 'dx_kT =', dx_kT
        write(*,*) 'det =', det
        write(*,*) 'dt =', dt
    end if
                                                                                                                
    ! Residual derivative
    dRiDcb = dkoffdf * DfDcb * (bterm - cbtau / (one - thetab)) - koff / (one - thetab)**2
    df = one - dRiDcb * dt

    ! Residual derivative wrt lambda_i at fixed cb: df/dlambda_i = ddwi / (lambda0 * r0), since dwi = lambda0 * r0 * fi
    dHdlambda = - dt * dkoffdf * ddwi / (lambda0 * r0) * (bterm - cbmax * thetab / (one - thetab))
    dHdth = - dt * koff * s * dbdth
    dHds  = - dt * koff * kon * cfmax * thetaf / (one - thetaf)
end subroutine kineticsFunc

