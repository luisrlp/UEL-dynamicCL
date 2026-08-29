subroutine kineticsFunc(cbtau, f, df, args, nargs)
    ! This subroutine serves as the function we would like to solve for                                         
    ! the bound crosslinker volume fraction (cbtau = cf/cfmax)                                                  
    ! by finding cbtau such that f = 0                                                                         
    use global                                                                                                        
    implicit none                                                                                               
                                                                                                                
    integer, intent(in)              :: nargs                                                                   
    DOUBLE PRECISION, intent(in out) :: cbtau                                                                   
    DOUBLE PRECISION, intent(out)    :: f, df                                                                    
    DOUBLE PRECISION, intent(in)     :: args(nargs)                                                              
                                                                                                                
    DOUBLE PRECISION                 :: r0f, etac, r0, r0c, fi, ffi, dwi, ddwi, l, mu0str, beta, b0 
    DOUBLE PRECISION                 :: cfmax, cbmax, dx_kT, dt, kon, koff0, koff, thetab, Ri
    DOUBLE PRECISION                 :: lambdai, lambdaif, lambda0, lambda0f, lambdaic, thetaf, cbt                                                           
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

    r0f = 1.6 * (cbtau*1.d3)**(- two / 5.d0)
    l = aratio * r0f
    r0 = r0f + r0c
    
    IF (lambdai.LE.one) then
        fi = 0.0
        DfDcb = 0.0
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
            ! CALL filpce(lambdai, fi, dwi, ddwi)
    END IF

    ! Unbinding rate
    koff = koff0 * exp(dx_kT * fi)
    thetab = cbtau / cbmax

    ! Reaction rate and residual                                                                         
    Ri = kon * cfmax * thetaf / (1 - thetaf) - koff * cbmax * thetab / (1 - thetab)
    f = cbtau - cbt - Ri * dt                                                      
                                                                                                                
    ! Residual derivative
    dRiDcb = - koff * (cbtau / (1 - thetab) * dx_kT * DfDcb + 1 / (1 - thetab)**2) 
    df = one - dRiDcb * dt
    
end subroutine kineticsFunc

