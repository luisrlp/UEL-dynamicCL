SUBROUTINE xmech(x, fle, dxdle, dfdle, cb, le, lamn, aratio, r0c, etac, mu0str, beta, b0, lambda0, kactin, rt)

    !> C4: x = mu_b^mech/RT = d[n (w + F0 ln lamn)]/dcb / RT at fixed le = lambda/lamn and lamn (material_TC._x_mech).
    !> w(le) = (lambda0 r0/aux) [f lf - int_0^f lam_f(f') df'] (exact extensible beta-model, 16-point Gauss-Legendre in
    !> s = ln(1 + aux1 f')). Also returns the force fle at le, dxdle = dx/dle = d(n w')/dcb / RT (w' = lambda0 r0 f) and
    !> dfdle = df/dle at fixed cb (for the tangent).

    use global
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(OUT) :: x, fle, dxdle, dfdle
    DOUBLE PRECISION, INTENT(IN)  :: cb, le, lamn, aratio, r0c, etac, mu0str, beta, b0, lambda0, kactin, rt
    INTEGER :: k
    DOUBLE PRECISION :: xg(16), wg(16), pi
    DOUBLE PRECISION :: r0f, dr0f, r0, ll, n, dndcb, aux, daux, l0f, dl0f, aux1, daux1, aux2, daux2, lf, kk
    DOUBLE PRECISION :: smax, s, fp, jac, u, v, dd, rhs, lam, drhs, dlam, qi, qj, w, dw, mu
    DOUBLE PRECISION :: fi, ffi, dwi, ddwi, dfdcb, f0, ff0, dw0, ddw0, df0
    DOUBLE PRECISION :: c_le, c_lf, c_l0, c_l0f, c_ll, c_r0, c_r0f, c_mu0, c_beta, c_b0, c_eta, c_cb
    xg = (/ -0.9894009349916499d0, -0.9445750230732326d0, -0.8656312023878318d0, -0.7554044083550030d0, &
            -0.6178762444026438d0, -0.4580167776572274d0, -0.2816035507792589d0, -0.0950125098376374d0, &
            0.0950125098376374d0, 0.2816035507792589d0, 0.4580167776572274d0, 0.6178762444026438d0, &
            0.7554044083550030d0, 0.8656312023878318d0, 0.9445750230732326d0, 0.9894009349916499d0 /)
    wg = (/ 0.0271524594117542d0, 0.0622535239386475d0, 0.0951585116824926d0, 0.1246289712555341d0, &
            0.1495959888165767d0, 0.1691565193950026d0, 0.1826034150449236d0, 0.1894506104550686d0, &
            0.1894506104550686d0, 0.1826034150449236d0, 0.1691565193950026d0, 0.1495959888165767d0, &
            0.1246289712555341d0, 0.0951585116824926d0, 0.0622535239386475d0, 0.0271524594117542d0 /)
    
    pi   = four * ATAN(one)
    r0f  = 1.6d0 * (cb * 1.0d3)**(-two / 5.0d0)
    dr0f = -0.4d0 * r0f / cb
    r0   = r0f + r0c
    ll   = aratio * r0f
    n    = kactin / ll
    dndcb   = 0.4d0 * n / cb
    aux  = etac * r0 / r0f
    daux = -etac * r0c * dr0f / r0f**2
    l0f  = aux * (lambda0 - one) + one
    dl0f = daux * (lambda0 - one)
    aux1    = ll * ll / (pi * pi * b0)
    daux1   = -0.8d0 * aux1 / cb
    aux2   = lambda0 * r0 / aux
    daux2  = lambda0 * (dr0f * aux - r0 * daux) / aux**2
    lf   = aux * (le - one) + one
    kk   = one - one / aratio
    w = zero; dw = zero; fle = zero; dxdle = zero; dfdle = zero

    IF (lf * l0f > one) THEN
        c_le = le; c_lf = lf; c_l0 = lambda0; c_l0f = l0f; c_ll = ll; c_r0 = r0; c_r0f = r0f
        c_mu0 = mu0str; c_beta = beta; c_b0 = b0; c_eta = etac; c_cb = cb
        CALL fil(fi, ffi, dwi, ddwi, c_le, c_lf, c_l0, c_l0f, c_ll, c_r0, c_r0f, c_mu0, c_beta, c_b0, c_eta, c_cb, dfdcb)
        fle   = fi
        dfdle = ddwi / (lambda0 * r0)                        ! dw' = lambda0 r0 df
        smax = LOG(one + aux1 * fi)
        qi = zero; qj = zero
        DO k = 1, 16
            s    = half * smax * (xg(k) + one)
            fp   = (EXP(s) - one) / aux1                          ! f' at the node
            jac  = EXP(s) / aux1                                  ! df'/ds
            u    = one + fp / mu0str
            v    = one + two * fp / mu0str
            dd   = one + fp * aux1 * u
            rhs  = u - kk * u**beta * v * dd**(-beta)
            lam  = aratio * rhs / l0f                          ! lam_f(f')
            drhs = beta * kk * u**beta * v * fp * u * dd**(-beta - one) * daux1
            dlam = aratio * drhs / l0f - lam * dl0f / l0f      ! d lam_f(f')/dcb at fixed f'
            qi = qi + wg(k) * lam * jac
            qj = qj + wg(k) * dlam * jac
        END DO
        qi = half * smax * qi
        qj = half * smax * qj
        w  = aux2 * (fi * lf - qi)
        dw = daux2 * (fi * lf - qi) + aux2 * (fi * daux * (le - one) - qj)
        dxdle = n * (0.4d0 * dwi / cb + lambda0 * (dr0f * fi + r0 * dfdcb)) / rt   ! d(n lambda0 r0 f)/dcb / RT
        END IF
        mu = dndcb * w + n * dw
        IF (lamn /= one) THEN                                  ! prestress reservoir: + n F0 ln(lamn)
        f0 = zero; df0 = zero
        IF (l0f > one) THEN
            c_le = one; c_lf = one; c_l0 = lambda0; c_l0f = l0f; c_ll = ll; c_r0 = r0; c_r0f = r0f
            c_mu0 = mu0str; c_beta = beta; c_b0 = b0; c_eta = etac; c_cb = cb
            CALL fil(f0, ff0, dw0, ddw0, c_le, c_lf, c_l0, c_l0f, c_ll, c_r0, c_r0f, c_mu0, c_beta, c_b0, c_eta, c_cb, df0)
        END IF
        mu = mu + LOG(lamn) * (dndcb * f0 * lambda0 * r0 + n * lambda0 * (dr0f * f0 + r0 * df0))
        END IF
        x = mu / rt
    RETURN
END SUBROUTINE xmech
