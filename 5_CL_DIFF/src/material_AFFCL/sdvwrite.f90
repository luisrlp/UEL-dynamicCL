    SUBROUTINE sdvwrite(det, statev, sigma, phi_tau, dmudx, Vmol, jfluid, cb, cb_tot)                                           
    !>    WRITE ALL STATE VARIABLES TO STATEV AT END OF INCREMENT                                                       
    !>
    !>    STATEV layout (defined in global.f90):
    !>      Slot  1       : phi_tau  (polymer volume fraction)
    !>      Slot  2       : det      (Jacobian J)
    !>      Slot  3       : c        (fluid content)
    !>      Slots 4-9     : sigma    (Cauchy stress, 6 components Voigt)
    !>      Slots 10-12   : -dmudx   (chemical potential gradient)
    !>      Slots 13-15   : jfluid   (fluid flux vector)
    !>      Slots 16-NSDV : cb(i)    (bound CL concentration per unique direction)
    use global
    implicit none
  
    DOUBLE PRECISION, INTENT(IN)  :: det
    DOUBLE PRECISION, INTENT(IN)  :: sigma(6)
    DOUBLE PRECISION, INTENT(IN)  :: phi_tau
    DOUBLE PRECISION, INTENT(IN)  :: dmudx(3,1), jfluid(3,1)
    DOUBLE PRECISION, INTENT(IN)  :: Vmol
    DOUBLE PRECISION, INTENT(IN)  :: cb(ndir), cb_tot
    DOUBLE PRECISION, INTENT(OUT) :: statev(nsdv)
  
    INTEGER :: idir
  
    ! --- Macroscopic quantities (slots 1-15, fixed layout) ---
    statev(1)     = phi_tau
    statev(2)     = det
    statev(3)     = (1.0d0 - phi_tau) / (Vmol * phi_tau * det)  ! fluid content c
    statev(4)     = cb_tot
    statev(5:10)   = sigma(1:6)        ! Cauchy stress (Voigt)
    statev(11:13) = -dmudx(1:3,1)    ! chemical potential gradient
    statev(14:16) = jfluid(1:3,1)    ! fluid flux vector


! Slots 17 to NSDV: bound crosslinker concentrations
DO idir = 1, ndir
    statev(nsdv - ndir + idir) = cb(idir)
END DO


! write(*,*) 'nsdv = ', nsdv
! write(*,*) 'statev = ', statev

RETURN

END SUBROUTINE sdvwrite

