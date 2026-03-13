      subroutine integ(props, nprops, dtime, F_tau, mu_tau, phi_t, theta, T_tau, SpTanMod, phi_tau)

         ! This subroutine computes everything required for the time integration
         ! of the problem.
         !
         ! Inputs:
         !  1) material parameters, props(nprops)
         !  2) time increment, dtime
         !  3) deformation gradient, F_tau(3,3)
         !  4) chemical potential, mu_tau
         !  5) old polymer volume fraction, phi_t
         !  6) temperature, theta
         !
         ! Outputs:
         !  1) Cauchy stress, T_tau(3,3)
         !  2) spatial tangent modulus, SpTanMod(3,3,3,3)
         !  3) polymer volume fraction, phi_tau
         !  4) time rate of polymer volume fraction, dPdt
         !  5) derivative of the phi with mu, DphiDmu
         !  6) derivative of the time rate of phi with mu, DphidotDmu
         !  7) scalar fluid permeability, Mfluid
         !  8) derivative of permeability with chemical potential, DmDmu
         !  9) volume of a mole of fluid, Vmol
         ! 10) displacement - chemical potential modulus terms
         ! 11) chemical potential - displacement modulus terms

         use global
         implicit none

         integer :: i, j, k, l, m, n, nprops, nargs, stat
         parameter(nargs=8)

         real(8) :: Iden(3,3), props(nprops), F_tau(3,3), phi_tau, mu_tau, phi_t
         real(8) :: theta, TR_tau(3,3), T_tau(3,3), dTRdF(3,3,3,3), Gshear, Kbulk
         real(8) :: SpTanMod(3,3,3,3), chi, D, mu0, Vmol, Rgas, detF, FinvT(3,3), dPdt
         real(8) :: B_tau(3,3), trB_tau, C_tau(3,3), trC_tau, args(nargs), detFe
         real(8) :: deltaMU, DphiDmu, dPdt_per, dPdt_m, DphidotDmu, Mfluid, Finv(3,3)
         real(8) :: phi_per, phi_m, dtime, DmDmu, DphiDJ, SpUCMod(3,3), DmDphi, DmDJ
         real(8) :: SpCUModFac(3,3), detFs

         ! Identity tensor
         call onem0(Iden)

         ! Compute the inverse of F, its determinant, and its transpose
         call matInv3Dd(F_tau, Finv, detF, stat)
         if (stat == 0) then
            write(*,*) 'Problem: detF.lt.zero'
            call exit
         endif
         FinvT = transpose(Finv)

         ! Compute the left Cauchy-Green tensor and its trace
         B_tau = matmul(F_tau, transpose(F_tau))
         trB_tau = B_tau(1,1) + B_tau(2,2) + B_tau(3,3)

         ! Compute the right Cauchy-Green tensor and its trace
         C_tau = matmul(transpose(F_tau), F_tau)
         trC_tau = C_tau(1,1) + C_tau(2,2) + C_tau(3,3)

         ! Compute the Cauchy stress
         T_tau = (Gshear * (B_tau - Iden) + Kbulk * log(detF) * Iden) / detF

         ! Compute the 1st Piola stress
         TR_tau = Gshear * (F_tau - FinvT) + Kbulk * log(detF) * FinvT

         ! Compute dTRdF, the so-called material tangent modulus
         dTRdF = 0.0
         do i = 1, 3
            do j = 1, 3
               do k = 1, 3
                  do l = 1, 3
                     dTRdF(i,j,k,l) = dTRdF(i,j,k,l) + Gshear * Iden(i,k) * Iden(j,l) &
                                              + Gshear * Finv(l,i) * Finv(j,k) &
                                              + Kbulk * Finv(j,i) * Finv(l,k) &
                                              - Kbulk * log(detF) * Finv(l,i) * Finv(j,k)
                  end do
               end do
            end do
         end do

         ! Calculate the so-called spatial tangent modulus, based
         ! on the push forward of the material tangent modulus
         SpTanMod = 0.0
         do i = 1, 3
            do j = 1, 3
               do k = 1, 3
                  do l = 1, 3
                     do m = 1, 3
                        do n = 1, 3
                           SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) + &
                                                         (dTRdF(i,m,k,n) * F_tau(j,m) * F_tau(l,n)) / detF
                        end do
                     end do
                  end do
               end do
            end do
         end do

         phi_tau = detF
         return
      end subroutine integ
