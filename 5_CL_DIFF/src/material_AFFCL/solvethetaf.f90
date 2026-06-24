subroutine solveThetaf(root, args, nargs, rootOld)

    ! This subroutine will numerically solve for the free
    ! crosslinker fraction (thetaf) based on the current
    ! chemical potential. See Numerical Recipes RTSAFE.

    implicit none

    ! 1. Dummy arguments explicitly strictly typed with INTENT
    integer, intent(in)     :: nargs
    real(8), intent(in)     :: args(nargs)
    real(8), intent(in) :: rootOld
    real(8), intent(out)    :: root

    ! 2. Local variables
    integer :: j
    real(8) :: f, df, fl, fh, xl, xh, x1, x2, swap, dxold
    real(8) :: dx, temp, rootMax, rootMin

    ! 3. Modern parameter declarations
    integer, parameter :: maxit = 50
    real(8), parameter :: xacc  = 1.0d-6
    real(8), parameter :: zero  = 0.0d0

    ! Set the safe bounds for thetaf (must be strictly between 0 and 1)
    rootMax = 1.0d0 - 1.0d-8
    rootMin = 1.0d-8

    x1 = rootMin
    x2 = rootMax
    call thetafFunc(x1, fl, df, args, nargs)
    call thetafFunc(x2, fh, df, args, nargs)

    ! Check if the root is safely bracketed
    if (fl * fh >= zero) then
        root = rootOld
        write(*,*) 'FYI, root not bracketed on thetaf'
        write(*,*) 'fl=', fl
        write(*,*) 'fh=', fh
        write(*,*) 'rootOld=', rootOld
        write(*,*) 'mu =', args(1)
        write(*,*) 'mu0=', args(2)
        write(*,*) 'Rgas=', args(3)
        write(*,*) 'theta=', args(4)
        write(*,*) 'chi=', args(5)
        write(*,*) 'Vmol=', args(6)
        write(*,*) 'Kbulk=', args(7)
        write(*,*) 'detF=', args(8)
        write(*,*) 'cb=', args(9)
        write(*,*) 'cfmax=', args(10)
        call exit
        return
    end if

    ! Orient the search so that f(xl) < 0
    if (fl < 0.0d0) then
        xl = x1
        xh = x2
    else
        xh = x1
        xl = x2
        swap = fl
        fl = fh
        fh = swap
    end if

    ! Initialize the guess for the root, the "step size before last", and the last step
    if (rootOld < rootMin) root = rootMin ! rootOld = rootMin
    if (rootOld > rootMax) root = rootMax ! rootOld = rootMax
    
    root  = rootOld
    dxold = abs(x2 - x1)
    dx    = dxold
    
    call thetafFunc(root, f, df, args, nargs)

    ! Loop over allowed iterations (Replaced old DO 10 loop)
    do j = 1, maxit
        
        ! Bisect if Newton is out of range, or not decreasing fast enough.
        if ( (((root - xh) * df - f) * ((root - xl) * df - f) >= 0.0d0) .or. &
             (abs(2.0d0 * f) > abs(dxold * df)) ) then

            dxold = dx
            dx    = 0.5d0 * (xh - xl)
            root  = xl + dx
            
            ! Change in root is negligible
            if (xl == root) return

        else
            ! Newton step is acceptable. Take it.
            dxold = dx
            dx    = f / df
            temp  = root
            root  = root - dx
            
            ! Change in root is negligible
            if (temp == root) return

        end if

        ! Convergence criterion
        if (abs(dx) < xacc) return

        ! The one new function evaluation per iteration
        call thetafFunc(root, f, df, args, nargs)

        ! Maintain the bracket on the root
        if (f < 0.0d0) then
            xl = root
            fl = f
        else
            xh = root
            fh = f
        end if

    end do

    ! If loop finishes without returning, maximum iterations were exceeded
    write(*, '(/1X,A)') 'solveThetaf EXCEEDING MAXIMUM ITERATIONS'
    
    return
end subroutine solveThetaf
