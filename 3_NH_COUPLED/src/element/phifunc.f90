    subroutine phiFunc(phi, f, df, args, nargs)
        ! This subroutine serves as the function we would like to solve for
        ! the polymer volume fraction by finding phi such that f = 0

        implicit none

        integer,           intent(in)    :: nargs
        real(8),          intent(inout) :: phi
        real(8),          intent(out)   :: f, df
        real(8),          intent(in)    :: args(nargs)

        real(8) :: mu, mu0, Rgas, theta, chi, Vmol, Kbulk
        real(8) :: detF, RT
        real(8), parameter :: zero  = 0.0
        real(8), parameter :: one   = 1.0
        real(8), parameter :: two   = 2.0
        real(8), parameter :: three = 3.0
        real(8), parameter :: third = 1.0 / 3.0

        ! Obtain relevant quantities
        mu    = args(1)
        mu0   = args(2)
        Rgas  = args(3)
        theta = args(4)
        chi   = args(5)
        Vmol  = args(6)
        Kbulk = args(7)
        detF  = args(8)

        ! Compute the useful quantity
        RT = Rgas * theta

        ! Compute the residual
        f = (mu0 - mu) / RT &
                + log(one - phi) + phi + chi * phi * phi &
                - ((Kbulk * Vmol) / RT) * log(detF * phi) &
                + ((Kbulk * Vmol) / (two * RT)) * (log(detF * phi)**two)

        ! Compute the tangent
        if (phi > 0.999_8) then
            df = zero
        else
            df = one - (one / (one - phi)) + two * chi * phi &
                    - (Kbulk * Vmol) / (RT * phi) &
                    + ((Kbulk * Vmol) / (RT * phi)) * log(detF * phi)
        end if

    end subroutine phiFunc

