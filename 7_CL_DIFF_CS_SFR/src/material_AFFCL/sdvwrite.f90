    SUBROUTINE sdvwrite(det, statev, sigma, thetaf, dmudx, Vmol, &
        jfluid, cb, cb_tot, cfmax, lamn)                                           
    !>    WRITE ALL STATE VARIABLES TO STATEV AT END OF INCREMENT                                                       
    !>
    !>    STATEV layout (defined in global.f90):
    !>      Slot  1               : thetaf  (thetaf)
    !>      Slot  2               : det      (Jacobian J)
    !>      Slot  3               : c        (fluid content)
    !>      Slot  4               : cb_tot   (total bound CL concentration)
    !>      Slots 5-10            : sigma    (Cauchy stress, 6 components Voigt)
    !>      Slots 11-13           : -dmudx   (chemical potential gradient)
    !>      Slots 14-16           : jfluid   (fluid flux vector)
    !>      Slots 17-16+NDIR      : cb(i)    (bound CL concentration per unique direction)
    !>      Slots 17+NDIR-NSDV    : lamn(i)  (reference stretch of stretch-free rebinding per unique direction)
    use global
    implicit none
  
    DOUBLE PRECISION, INTENT(IN)  :: det
    DOUBLE PRECISION, INTENT(IN)  :: sigma(6)
    DOUBLE PRECISION, INTENT(IN)  :: thetaf
    DOUBLE PRECISION, INTENT(IN)  :: dmudx(3,1), jfluid(3,1)
    DOUBLE PRECISION, INTENT(IN)  :: Vmol
    DOUBLE PRECISION, INTENT(IN)  :: cb(ndir), cb_tot
    DOUBLE PRECISION, INTENT(IN)  :: cfmax
    DOUBLE PRECISION, INTENT(IN)  :: lamn(ndir)
    DOUBLE PRECISION, INTENT(OUT) :: statev(nsdv)
  
    INTEGER :: idir
  
    ! --- Macroscopic quantities (slots 1-15, fixed layout) ---
    statev(1)     = thetaf
    statev(2)     = det
    statev(3)     = thetaf * cfmax + cb_tot  ! fluid content c
    statev(4)     = cb_tot
    statev(5:10)   = sigma(1:6)        ! Cauchy stress (Voigt)
    statev(11:13) = -dmudx(1:3,1)    ! chemical potential gradient
    statev(14:16) = jfluid(1:3,1)    ! fluid flux vector


! Bound crosslinkers and reference stretch
DO idir = 1, ndir
    statev(nsdv - 2*ndir + idir) = cb(idir)
    statev(nsdv - ndir + idir) = lamn(idir)
END DO


! write(*,*) 'nsdv = ', nsdv
! write(*,*) 'statev = ', statev

RETURN

END SUBROUTINE sdvwrite

