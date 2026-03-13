subroutine pk2anisomatfic(afic, daniso, cbar, inv4, st0, ndi)
      !> ANISOTROPIC MATRIX: 2PK 'FICTICIOUS' STRESS TENSOR
      !! INPUT:
      !! DANISO - ANISOTROPIC STRAIN-ENERGY DERIVATIVES
      !! CBAR - DEVIATORIC LEFT CAUCHY-GREEN TENSOR
      !! INV1, INV4 - CBAR INVARIANTS
      !! UNIT2 - 2ND ORDER IDENTITY TENSOR
      !! OUTPUT:
      !! AFIC - 2ND PIOLA KIRCHOOF 'FICTICIOUS' STRESS TENSOR
      use global
      implicit none

      integer, intent(in) :: ndi
      real(8), intent(in) :: daniso(3), cbar(3, 3), inv4
      real(8), intent(in) :: st0(3, 3)
      real(8), intent(out) :: afic(ndi, ndi)
      real(8) :: dudi4, di4dc(3, 3)

      ! FIRST DERIVATIVE OF SSEANISO IN ORDER TO I4
      dudi4 = daniso(1)

      di4dc = st0

      afic = 2.0d0 * (dudi4 * di4dc)

      return
end subroutine pk2anisomatfic
