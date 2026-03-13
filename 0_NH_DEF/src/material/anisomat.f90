subroutine anisomat(sseaniso, daniso, diso, k1, k2, kdisp, i4, i1)
      use global
      implicit none

      ! Arguments
      double precision, intent(out) :: sseaniso
      double precision, intent(out) :: daniso(4)
      double precision, intent(inout) :: diso(5)
      double precision, intent(in) :: k1, k2, kdisp, i4, i1

      ! Local variables
      double precision :: dudi1, d2ud2i1
      double precision :: e1, ee2, ee3, dudi4, d2ud2i4, d2dudi1di4, d2dudi2di4

      dudi1 = diso(1)
      d2ud2i1 = diso(3)

      e1 = i4 * (one - three * kdisp) + i1 * kdisp - one
      sseaniso = (k1 / k2) * (dexp(k1 * e1 * e1) - one)

      if (e1 > zero) then
            ee2 = dexp(k2 * e1 * e1)
            ee3 = (one + two * k2 * e1 * e1)

            dudi1 = dudi1 + k1 * kdisp * e1 * ee2
            d2ud2i1 = d2ud2i1 + k1 * kdisp * kdisp * ee3 * ee2

            dudi4 = k1 * (one - three * kdisp) * e1 * ee2
            d2ud2i4 = k1 * ((one - three * kdisp)**two) * ee3 * ee2
            d2dudi1di4 = k1 * (one - three * kdisp) * kdisp * ee3 * ee2
            d2dudi2di4 = zero
      else
            dudi4 = zero
            d2ud2i4 = zero
            d2dudi1di4 = zero
            d2dudi2di4 = zero
            d2ud2i1 = zero
      end if

      ! First derivative of sseaniso with respect to i1
      daniso(1) = dudi4
      ! First derivative of sseaniso with respect to i2
      daniso(2) = d2ud2i4
      ! Second derivative of sseaniso with respect to i1
      daniso(3) = d2dudi1di4
      ! Second derivative of sseaniso with respect to i2
      daniso(4) = d2dudi2di4

      diso(1) = dudi1
      diso(3) = d2ud2i1

      return
end subroutine anisomat
