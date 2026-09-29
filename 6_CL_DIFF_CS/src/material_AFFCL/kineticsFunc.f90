subroutine kineticsFunc(cbtau, f, df, dHdl, args, nargs)
    ! This subroutine serves as the function we would like to solve for
    ! the bound crosslinker volume fraction (cbtau = cf/cfmax)
    ! by finding cbtau such that f = 0
    ! Outputs: f = H (Backward-Euler residual), df = dH/dcb, dHdl = dH/dlambda_i (at fixed cb),
    ! the latter used for the consistent tangent (sub-step sensitivity dcb/dlambda_i)
    use global
    implicit none

    integer, intent(in)              :: nargs
    DOUBLE PRECISION, intent(in out) :: cbtau
    DOUBLE PRECISION, intent(out)    :: f, df, dHdl
    DOUBLE PRECISION, intent(in)     :: args(nargs)                                                              
                                                                                                                
    DOUBLE PRECISION                 :: r0f, etac, r0, r0c, fi, ffi, dwi, ddwi, l, mu0str, beta, b0 
    DOUBLE PRECISION                 :: cfmax, cbmax, dx_kT, dt, kon, koff0, koff, thetab, Ri
    DOUBLE PRECISION                 :: kcatch0, dxc_kT, dkoffdf
    DOUBLE PRECISION                 :: lambdai, lambdaif, lambda0, lambda0f, lambdaic, thetaf, cbt, det                                                          
    DOUBLE PRECISION                 :: DfDcb,DRiDcb, aratio
    
    Ri = zero
                                                                                                                
    ! Obtain relevant quantities
    lambdai = args(1)
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

    r0f = 1.6 * (cbtau*1.d3)**(- two / 5.d0)
    l = aratio * r0f
    r0 = r0f + r0c
    
    IF (lambdai.LE.one) then
        fi = 0.0
        DfDcb = 0.0
        ddwi = zero
    ELSE
        IF((etac > zero).AND.(etac .LE. one))THEN
            lambdaif=etac*(r0/r0f)*(lambdai-one)+one
            lambda0f=etac*(r0/r0f)*(lambda0-one)+one
            lambdaic=(lambdai*r0-lambdaif*r0f)/r0c
        ELSE
            lambdaif=lambdai ! False for a filament attached to a stiff crosslinker (etac = 1), only valid for etac = 0 (???)
            lambdaic=zero ! False for a stiff crosslinker (etac = 1), only valid for etac = 0 (???)
        END IF
        CALL fil(fi,ffi,dwi,ddwi,&
                lambdai,lambdaif,lambda0,lambda0f,&
                l,r0,r0f,mu0str,beta,b0,etac,&
                cbtau,DfDcb)
        !CALL fil_inext(fi,dwi,ddwi,lambdai,lambdaif,lambda0,lambda0f,l,r0,r0f,beta,b0,etac,cbtau,DfDcb)
            ! CALL filpce(lambdai, fi, dwi, ddwi)
    END IF

    ! Unbinding rate (catch-slip) and its force derivative
    CALL koffcs(koff, dkoffdf, fi, koff0, kcatch0, dx_kT, dxc_kT)
    thetab = cbtau / cbmax

    ! Reaction rate and residual                                                                         
    Ri = kon * cfmax * thetaf / (1 - thetaf) - koff * cbmax * thetab / (1 - thetab)
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
    dRiDcb = - (dkoffdf * cbtau / (1 - thetab) * DfDcb + koff / (1 - thetab)**2)
    df = one - dRiDcb * dt

    ! Residual derivative wrt lambda_i at fixed cb: df/dlambda_i = ddwi / (lambda0 * r0), since dwi = lambda0 * r0 * fi
    dHdl = dt * cbmax * thetab / (1 - thetab) * dkoffdf * ddwi / (lambda0 * r0)

end subroutine kineticsFunc

