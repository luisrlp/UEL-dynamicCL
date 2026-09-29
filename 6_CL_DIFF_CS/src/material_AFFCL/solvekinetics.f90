subroutine solveKinetics(root, args, nargs, rootOld, converged, atbound, dHdcb, dHdl)

    ! Implicit (Backward-Euler) update of the bound concentration of one direction:
    ! root of H(cb) = cb - cbt - dt*R(cb) (kineticsFunc), same algorithm as material.py.
    !
    ! Catch bonds can make H(cb) non-monotonic (dH/dcb -> 0 or < 0), with several roots.
    ! The solution must continue the previous state: starting at rootOld, H is marched in
    ! the direction of the reaction (upwards if H(rootOld) < 0) with growing multiplicative
    ! steps up to the FIRST sign change, and the root is then found in that bracket with a
    ! safeguarded Newton (Numerical Recipes RTSAFE). No intermediate root (other branch)
    ! can be skipped, up to the resolution of the marching steps.
    !
    ! Outputs: root, converged (.false. if H is NaN/Inf, root = rootOld, or RTSAFE did not
    ! converge; if H keeps its sign up to a bound, the root lies beyond it and the bound is
    ! returned as converged, with atbound = .true.), dHdcb and dHdl (dH/dcb and dH/dlambda_i)
    ! at root, used for the sub-stepping checks and the consistent tangent.
    ! No program stop: the caller decides (sub-stepping / PNEWDT cut-back).

    implicit none

    ! Dummy arguments
    integer, intent(in)     :: nargs
    real(8), intent(in)     :: args(nargs)
    real(8), intent(in)     :: rootOld
    real(8), intent(out)    :: root
    logical, intent(out)    :: converged
    logical, intent(out)    :: atbound
    real(8), intent(out)    :: dHdcb, dHdl

    ! Local variables
    integer :: j, imarch
    logical :: upwards, newton_ok
    real(8) :: rootMax, rootMin, x0, xend, a, b, s
    real(8) :: Ha, dHa, Hb, dHb, dum
    real(8) :: xl, xh, x, dx, dxold, Hx, dHx, xnew
    real(8) :: last_x, last_H, last_dH, last_dHdl

    ! Parameter declarations
    integer, parameter :: maxit = 50       ! RTSAFE iterations
    integer, parameter :: maxmarch = 200   ! marching steps (bound to bound needs ~20)
    real(8), parameter :: xacc  = 1.0d-12  ! tolerance on |H| and on the root step
    real(8), parameter :: s0    = 0.05d0   ! first relative marching step
    real(8), parameter :: sgrow = 1.5d0    ! growth factor of the marching step

    ! Safe bounds for the root (as material.py)
    rootMin = 1.0d-10
    rootMax = args(9) - 1.0d-10            ! cbmax - 1e-10

    last_x = -1.0d0
    converged = .false.
    atbound = .false.

    x0 = min(max(rootOld, rootMin), rootMax)
    root = x0
    call evalH(x0, Ha, dHa)

    if (abs(Ha) < xacc) then
        converged = .true.
    else if (abs(Ha) < huge(Ha)) then
        ! March from rootOld in the direction of the reaction up to the first sign change of H
        upwards = (Ha < 0.0d0)
        if (upwards) then
            xend = rootMax
        else
            xend = rootMin
        end if
        a = x0
        s = s0
        do imarch = 1, maxmarch
            if (upwards) then
                b = min(a * (1.0d0 + s), xend)
            else
                b = max(a / (1.0d0 + s), xend)
            end if
            call evalH(b, Hb, dHb)
            if (.not. (abs(Hb) <= huge(Hb))) exit          ! NaN/Inf: give up (not converged)
            ! Sign change (no product Ha*Hb: both can be huge and the product would overflow)
            if ((Ha <= 0.0d0 .and. Hb >= 0.0d0) .or. (Ha >= 0.0d0 .and. Hb <= 0.0d0)) then
                call rtsafe(min(a, b), max(a, b), a, Ha, dHa, merge(Ha, Hb, a < b))
                exit
            end if
            if (b == xend) then
                ! No sign change up to the bound: the root lies beyond it, between the bound and
                ! 0 (H(0) = -cbt - dt*Ron < 0) or cbmax (H -> +inf). This happens e.g. at very
                ! large forces (complete unbinding, root ~ cbt/(dt*koff) << rootMin). The bound
                ! is the solution to within 1e-10 (the bound distance): accept it as converged.
                root = xend
                converged = .true.
                atbound = .true.
                exit
            end if
            a = b
            Ha = Hb
            dHa = dHb
            s = sgrow * s
        end do
    end if

    ! dH/dcb and dH/dlambda_i at the root (reuse the last evaluation when it is at the root)
    if (last_x == root) then
        dHdcb = last_dH
        dHdl = last_dHdl
    else
        call evalH(root, dum, dum)
        dHdcb = last_dH
        dHdl = last_dHdl
    end if

    return

contains

    subroutine evalH(xin, H, dH)
        ! kineticsFunc wrapper that keeps the last evaluation (x, H, dH/dcb, dH/dlambda_i)
        real(8), intent(in)  :: xin
        real(8), intent(out) :: H, dH
        real(8) :: xloc
        xloc = xin
        call kineticsFunc(xloc, H, dH, last_dHdl, args, nargs)
        last_x = xin
        last_H = H
        last_dH = dH
    end subroutine evalH

    subroutine rtsafe(x1, x2, xstart, Hstart, dHstart, H1)
        ! Safeguarded Newton (Numerical Recipes RTSAFE) on [x1, x2], which brackets a sign
        ! change of H. Starts at xstart, a bracket end where H and dH/dcb are already known
        ! (Hstart, dHstart); H1 = H(x1) is also known, so no evaluation is repeated.
        ! Sets root and converged of the host.
        real(8), intent(in) :: x1, x2, xstart, Hstart, dHstart, H1

        if (H1 < 0.0d0) then
            xl = x1
            xh = x2
        else
            xl = x2
            xh = x1
        end if
        x = xstart
        Hx = Hstart
        dHx = dHstart
        dxold = abs(x2 - x1)
        dx = dxold

        do j = 1, maxit
            if (abs(Hx) < xacc) then
                root = x
                converged = .true.
                return
            end if
            ! Newton only if finite, inside the bracket and decreasing fast enough; else bisect.
            ! Written in overflow-safe form: the Newton step Hx/dHx is only formed once it is known
            ! to be shorter than the bracket; then xnew inside (xl, xh) is the RTSAFE range test.
            newton_ok = (abs(Hx) < huge(Hx)) .and. (abs(dHx) < huge(dHx)) .and. (dHx /= 0.0d0)
            if (newton_ok) newton_ok = (abs(Hx) < abs(dHx) * abs(xh - xl))
            if (newton_ok) then
                xnew = x - Hx / dHx
                newton_ok = ((xnew - xl) * (xnew - xh) < 0.0d0) .and. (abs(Hx) <= 0.5d0 * abs(dxold * dHx))
            end if
            dxold = dx
            if (newton_ok) then
                dx = Hx / dHx
                x = x - dx
            else
                dx = 0.5d0 * (xh - xl)
                x = xl + dx
            end if
            call evalH(x, Hx, dHx)
            if (abs(dx) < xacc) then
                root = x
                converged = (abs(Hx) < huge(Hx))
                return
            end if
            if (Hx < 0.0d0) then
                xl = x
            else
                xh = x
            end if
        end do

        ! Maximum iterations exceeded
        root = x
        converged = .false.
    end subroutine rtsafe

end subroutine solveKinetics
