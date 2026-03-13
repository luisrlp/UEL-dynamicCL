subroutine CMATANISOMATFIC(CMANISOMATFIC, M0, DANISO, UNIT2, DET, NDI)
      ! ANISOTROPIC MATRIX: MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
      use global
      implicit none

      integer :: NDI, I, J, K, L
      real(8) :: CMANISOMATFIC(NDI,NDI,NDI,NDI), UNIT2(NDI,NDI)
      real(8) :: M0(NDI,NDI), DANISO(3), DET
      real(8) :: CINV4(NDI,NDI,NDI,NDI), CINV14(NDI,NDI,NDI,NDI)
      real(8) :: D2UDI4, D2UDI1DI4
      real(8) :: IMM(NDI,NDI,NDI,NDI), MMI(NDI,NDI,NDI,NDI)
      real(8) :: MM0(NDI,NDI,NDI,NDI)

      ! 2ND DERIVATIVE OF SSEANISO IN ORDER TO I4
      D2UDI4 = DANISO(2)
      ! 2ND DERIVATIVE OF SSEANISO IN ORDER TO I1 AND I4
      D2UDI1DI4 = DANISO(3)

      call TENSORPROD2(M0, M0, MM0, NDI)
      call TENSORPROD2(UNIT2, M0, IMM, NDI)
      call TENSORPROD2(M0, UNIT2, MMI, NDI)

      do I = 1, NDI
            do J = 1, NDI
                  do K = 1, NDI
                        do L = 1, NDI
                              CINV4(I,J,K,L) = D2UDI4 * MM0(I,J,K,L)
                              CINV14(I,J,K,L) = D2UDI1DI4 * (IMM(I,J,K,L) + MMI(I,J,K,L))
                              CMANISOMATFIC(I,J,K,L) = FOUR * (CINV4(I,J,K,L) + CINV14(I,J,K,L))
                        end do
                  end do
            end do
      end do

      return
end subroutine CMATANISOMATFIC

