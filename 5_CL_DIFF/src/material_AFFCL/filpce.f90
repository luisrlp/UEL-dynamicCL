subroutine filpce(q0, force, dw, ddw)

    !!! Don't forget to apply exponential scaling to outputs !!!


    implicit none
    real(8), intent(in) :: q0
    real(8), intent(out) :: force, dw, ddw
    real(8) :: f_std, dw_std, ddw_std
    real(8), dimension(10) :: thresholds ! = (/1.001D0, 1.002D0, 1.003D0, 1.004D0, 1.005D0, 1.006D0, 1.007D0, 1.008D0, 1.009D0, 1.01D0/)
    real(8), dimension(10) :: force_values ! = (/0.00033D0, 0.00108D0, 0.00200D0, 0.00275D0, 0.00368D0, 0.00405D0, 0.00489D0, 0.00559D0, 0.00620D0, 0.007370D0/)
    real(8), dimension(10) :: dw_values ! = (/0.00053D0, 0.00174D0, 0.00325D0, 0.00450D0, 0.00610D0, 0.00670D0, 0.00808D0, 0.00910D0, 0.00980D0, 0.00120D0/)
    real(8), dimension(10) :: ddw_values ! = (/1.1300D0, 1.2000D0, 1.2900D0, 1.3100D0, 1.3100D0, 1.3100D0, 1.3100D0, 1.3100D0, 1.3100D0, 1.3100D0/)
    integer :: i

    f_std = 0.0
    dw_std = 0.0
    ddw_std = 0.0
    thresholds = (/1.001D0, 1.002D0, 1.003D0, 1.004D0, 1.005D0, 1.006D0, 1.007D0, 1.008D0, 1.009D0, 1.01D0/)
    force_values = (/0.00033D0, 0.00108D0, 0.00200D0, 0.00275D0, 0.00368D0, 0.00405D0, 0.00489D0, 0.00559D0, 0.00620D0, 0.007370D0/)
    dw_values = (/0.00053D0, 0.00174D0, 0.00325D0, 0.00450D0, 0.00610D0, 0.00670D0, 0.00808D0, 0.00910D0, 0.00980D0, 0.00120D0/)
    ddw_values = (/1.1300D0, 1.2000D0, 1.2900D0, 1.3100D0, 1.3100D0, 1.3100D0, 1.3100D0, 1.3100D0, 1.3100D0, 1.3100D0/)
    
    ! call linear_interpolate(q0, f_std, dw_std, ddw_std)

    ! if (q0 >= 1.0D0 .AND. q0 < 1.05D0) then 
    !     force = 3.143817878145808d0*q0**2-5.5450506678545555d0*q0+2.4010851592154583d0
    !     dw = 5.193416679990963d0*q0**2-9.165195610413837d0*q0+3.9715403544588344d0
    !     ddw = 54.99363826460525d0*q0**2-100.99997536261706d0*q0+47.22043323959344d0

    ! else if (q0 >= 1.05D0 .AND. q0 < 1.1D0) then
    !     force = 7.896205363263728d0*q0**2-15.54807109600342d0*q0+7.665454844613117d0
    !     dw = 13.229624521389654d0*q0**2-26.082696045354645d0*q0+12.876203946089644d0
    !     ddw = 212.43640751141814d0*q0**2-432.72558760677595d0*q0+221.982176930368d0

    ! else if (q0 >= 1.1D0 .AND. q0 < 1.12D0) then
    !     force = 11.011913045449417d0*q0**2-22.520969697374664d0*q0+11.565896009642872d0
    !     dw = 17.311746766183273d0*q0**2-35.2551371685696d0*q0+18.02686633962719d0
    !     ddw = 276.30227891320976d0*q0**2-578.2204744989647d0*q0+304.76232751727565d0
        
    ! else if (q0 >= 1.12D0 .AND. q0 < 1.14D0) then
    !     force = -36.830475981616196d0*q0**2+85.92882998048795d0*q0-49.88814711175428d0
    !     dw = -64.7313899344681d0*q0**2+150.75156065915147d0*q0-87.39208966212128d0
    !     ddw = -2438.956034176943d0*q0**2+5577.27060372108d0*q0-3183.5982113828027d0
        
    ! else if (q0 >= 1.14D0 .AND. q0 <= 1.15D0) then
    !     force = -163.7191055592126d0*q0**2+378.77618279620384d0*q0-218.8358100817485d0
    !     dw = -288.20701242487434d0*q0**2+666.4427935792925d0*q0-384.86181405032886d0
    !     ddw = -9407.226531397582d0*q0**2+21665.548471067206d0*q0-12468.657545978633d0

    if (q0 >= 1.0D0 .AND. q0 < 1.15D0) then 

        force = 28.386319169362196D0*q0**3-84.37306994347098D0*q0**2+84.38968434687419d0*q0-28.40360003276212d0
        dw = 46.610613481110704D0*q0**3-138.43961073313733D0*q0**2+138.3632113963927D0*q0-46.535279845964695D0
        ddw = 856.696285474769D0*q0**3-2586.113628014857D0*q0**2+2613.065702768421D0*q0-882.4573288455576D0

    else 
        force = 1.0
        dw = 1.0
        ddw = 1.0

    end if

    ! Uncomment to run simulation with "extreme" values (expected values +- std)
    ! force = force + f_std
    ! dw = dw + dw_std
    ! ddw = ddw + ddw_std

    ! To handle negative outputs
    ! if (force < 0.0) then
    !     force = 0.0D0
    !     do i = 1, size(thresholds)
    !         if (q0 < thresholds(i)) then
    !             force = force_values(i)
    !             exit
    !         end if
    !     end do
    ! end if

    ! if (dw < 0.0) then
    !     dw = 0.0D0
    !     do i = 1, size(thresholds)
    !         if (q0 < thresholds(i)) then
    !             dw = dw_values(i)
    !             exit
    !         end if
    !     end do
    ! end if

    ! if (ddw < 0.0) then
    !     ddw = 0.0D0
    !     do i = 1, size(thresholds)
    !         if (q0 < thresholds(i)) then
    !             ddw = ddw_values(i)
    !             exit
    !         end if
    !     end do
    ! end if

    ! if (force < 0.0) then
    !     force = 0.0
    ! end if

    ! if (dw < 0.0) then
    !     dw = 0.0
    ! end if

    ! if (ddw < 0.0) then
    !     ddw = 0.0
    ! end if

    ! write(*,*) "lambda, Force: ", q0, force
end subroutine filpce

