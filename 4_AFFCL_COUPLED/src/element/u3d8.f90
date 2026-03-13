      subroutine U3D8(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, &
                  PROPS, NPROPS, coords, MCRD, NNODE, Uall, DUall, Vel, Accn, JTYPE, &
                  TIME, DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG, &
                  PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS, &
                  NJPROP, PERIOD, nDim, nIntt, nInttS)

         use global
         implicit none

         ! Variables defined in UEL, passed back to Abaqus
         real(8), intent(out) :: RHS(MLVARX, 1), AMATRX(NDOFEL, NDOFEL), SVARS(NSVARS), ENERGY(8)

         ! Variables passed into UEL
         real(8), intent(in) :: PROPS(NPROPS), coords(MCRD, NNODE), Uall(NDOFEL), DUall(MLVARX, 1), &
                           Vel(NDOFEL), Accn(NDOFEL), TIME(2), DTIME, PARAMS(1), ADLMAG(MDLOAD, 1), &
                           PREDEF(2, NPREDF, NNODE), DDLMAG(MDLOAD, 1), PERIOD
         integer, intent(in) :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE, KSTEP, KINC, &
                           JELEM, NDLOAD, JDLTYP(MDLOAD, 1), NPREDF, LFLAGS(4), MLVARX, MDLOAD, &
                           JPROPS(NJPROP), NJPROP

         ! Local variables
         real(8) :: u(NNODE, 3), du(NNODE, NDOFEL), thetaNew(NNODE), thetaOld(NNODE), &
                  dtheta(NNODE), muNew(NNODE), muOld(NNODE), dMU(NNODE), uNew(NNODE, NDOFEL), &
                  uOld(NNODE, NDOFEL), u_t(NNODE, NDOFEL), v(NNODE, 3), coordsC(MCRD, NNODE)
         integer :: i, j, k, l, m, n, nInttPt, nDim, intpt, pOrder, face, nIntt, ii, jj, pe, stat, q, &
                  nInttV, nInttPtV, p, ngSdv, nlSdv, kk, lenJobName, lenOutDir, nInttS, faceFlag, &
                  nshr, ntens
         real(8) :: statev(nsdv), prev_statev(nsdv), Iden(3, 3), Le, theta0, phi0, Ru(3 * NNODE, 1), Rc(NNODE, 1), &
                  body(3), Kuu(3 * NNODE, 3 * NNODE), Kcc(NNODE, NNODE), sh0(NNODE), detMapJ0, &
                  dshxi(NNODE, 3), dsh0(NNODE, 3), dshC0(NNODE, 3), detMapJ0C, Vmol, Fc_tau(3, 3), &
                  Fc_t(3, 3), detFc_tau, detFc_t, w(nIntt), DmDmu, DmDJ, sh(NNODE), detMapJ, phi_t, &
                  dsh(NNODE, 3), detMapJC, phiLmt, umeror, dshC(NNODE, 3), mu_tau, mu_t, dMUdX(3, 1), &
                  dMUdt, F_tau(3, 3), F_t(3, 3), detF_tau, xi(nIntt, 3), detF, TR_tau(3, 3), T_tau(3, 3), &
                  xi0(nIntt, 3), Ff_t(3, 3), Ff_tau(3, 3), SpTanMod(3, 3, 3, 3), phi_tau, dPdt, DphiDmu, &
                  DphidotDmu, Mfluid, Smat(6, 1), Bmat(6, 3 * NNODE), BodyForceRes(3 * NNODE, 1), flux, &
                  Gmat(9, 3 * NNODE), G0mat(9, 3 * NNODE), Amat(9, 9), Qmat(9, 9), dA, xLocal(nInttS), &
                  yLocal(nInttS), zLocal(nInttS), wS(nInttS), Kuc(3 * NNODE, NNODE), Kcu(NNODE, 3 * NNODE), &
                  Nvec(1, NNODE), ResFac, AmatUC(6, 1), TanFac, AmatCU(3, 9), SpUCMod(3, 3), &
                  SpCUMod(3, 3, 3), SpCUModFac(3, 3), pi, detF_t, PNEWDT
         character(len=256) :: jobName, outDir, fileName

         ! Get element parameters
         nlSdv = JPROPS(1) ! number of local sdv's per integ point
         ngSdv = JPROPS(2) ! number of global sdv's per integ point

         ! Allocate memory for the globalSdv's used for viewing results 
         ! on the dummy mesh
         pi = 4.0d0 * atan(1.0d0)
         xi0 = 0.0d0

         ! Initialize energy
         ENERGY = 0.0d0

         if (.not. allocated(globalSdv)) then
            ! Allocate memory for the globalSdv's
            !
            ! numElem needs to be set in the MODULE
            ! nIntt needs to be set in the UEL
            !
            allocate(globalSdv(numElem, nIntt, ngSdv), stat=err)
            if (err /= 0) then
               write(*, *) '//////////////////////////////////////////////'
               write(*, *) 'error when allocating globalSdv'
               write(*, *) '//////////////////////////////////////////////'
               write(*, *) '   stat=', stat
               write(*, *) '  ngSdv=', ngSdv
               write(*, *) '   nIntt=', nIntt
               write(*, *) 'numElem=', numElem
               write(*, *) '  nNode=', nNode
               write(*, *) 'lbound(globalSdv)=', lbound(globalSdv)
               write(*, *) 'ubound(globalSdv)=', ubound(globalSdv)
               write(*, *) "Error during allocation. Error code:", err
               write(*, *) '//////////////////////////////////////////////'
               call exit
            endif
            write(*, *) '-------------------------------------------------'
            write(*, *) '----------- globalSDV ALLOCATED -----------------'
            write(*, *) '-------------------------------------------------'
            write(*, *) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
            write(*, *) '---------- numElem=', numElem
            write(*, *) '---------- U3D8 ELEMENTS ------------------------'
            write(*, *) '-------------------------------------------------'
            write(*, *) '---------- YOU PUT NUMBER OF POINTS -------------'
            write(*, *) '---------- nIntt =', nIntt
            write(*, *) '---------- nInttS=', nInttS
            write(*, *) '-------------------------------------------------'
            write(*, *) '---------- YOU PUT NUMBER OF SDVs ---------------'
            write(*, *) '---------- ngSdv=', ngSdv
            write(*, *) '-------------------------------------------------'
         endif


      !      write(*,*) 'NDOFEL',NDOFEL
      !      write(*,*) 'MLVARX',MLVARX
      !      write(*,*) 'NRHS',NRHS
      !      write(*,*) 'NSVARS',NSVARS
      !      write(*,*) 'NPROPS',NPROPS
      !      write(*,*) 'NJPROP',NJPROP
      !      write(*,*) 'MCRD',MCRD
      !      write(*,*) 'NNODE',NNODE
      !      write(*,*) 'JTYPE',JTYPE
      !      write(*,*) 'KSTEP',KSTEP
      !      write(*,*) 'KINC',KINC
      !      write(*,*) 'JELEM',JELEM
      !      write(*,*) 'NDLOAD',NDLOAD
      !      write(*,*) 'MDLOAD',MDLOAD
      !      write(*,*) 'NPREDF',NPREDF
      !      write(*,*) '#################################'
         
      ! Identity tensor
      !
         call onem0(Iden)


         ! Obtain initial conditions
         !
         theta0 = props(20)
         phi0   = props(21)

         ! Initialize the residual and tangent matrices to zero.
         Ru = zero
         Rc = zero
         Kuu = zero
         Kcc = zero
         Kuc = zero
         Kcu = zero
!      Energy = zero

         ! Body forces
         body(1:3) = 0.0d0

         ! Obtain nodal displacements
         k = 0
         do i = 1, NNODE
            do j = 1, nDim
               k = k + 1
               u(i, j) = Uall(k)
               du(i, j) = DUall(k, 1)
               uOld(i, j) = u(i, j) - du(i, j)
            end do
            k = k + 1
            muNew(i) = Uall(k)
            dMU(i) = DUall(k,1)
            muOld(i) = muNew(i) - dMU(i)
         end do

         ! Obtain current nodal coordinates
         do i = 1, NNODE
            do j = 1, nDim
               coordsC(j, i) = coords(j, i) + u(i, j)
            end do
         end do

         ! Impose any time-stepping changes on the increments of chemical potential or displacement if you want
         !
         ! chemical potential increment
         !
         do i=1,nNode
            if(dabs(dMU(i)).gt.1.d6) then
               pnewdt = 0.5
               return
            endif
         enddo
         !
         ! displacement increment, based on element diagonal
         !
         Le = sqrt((coordsC(1, 1) - coordsC(1, 7))**2 + &
                 (coordsC(2, 1) - coordsC(2, 7))**2 + &
                 (coordsC(3, 1) - coordsC(3, 7))**2)
         ! add some kind of flag here???
         do i = 1, NNODE
            do j = 1, nDim
               if (abs(du(i, j)) > 10.0d0 * Le) then
                  PNEWDT = 0.5d0
                  return
               endif
            end do
         end do

            !----------------------------------------------------------------
            ! 
            ! Take this opportunity to perform calculations at the element
            !  centroid.  Here, check for hourglass stabilization and get
            !  the deformation gradient for use in the `F-bar' method.
            !
            ! Reference for the F-bar method:
            !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
            !  Design of simple low order finite elements for large strain
            !  analysis of nearly incompressible solids. International Journal
            !  of Solids and Structures, 33, 3277-3296.
            !
            !
            ! Obtain shape functions and their local gradients at the element
            !  centroid, that means xi=eta=zeta=0.0, and nInttPt=1
            !
            if (nNode == 8) then
               call calcShape3DLinear(1, xi0, 1, sh0, dshxi)
            else
               write(*, *) 'Incorrect number of nodes: nNode.ne.8'
               call exit
            endif

            ! Map shape functions from local to global reference coordinate system
            !
            call mapShape3D(nNode, dshxi, coords, dsh0, detMapJ0, stat)
            if (stat == 0) then
               PNEWDT = 0.5
               return
            endif

            ! Map shape functions from local to global current coordinate system
            !
            call mapShape3D(nNode, dshxi, coordsC, dshC0, detMapJ0C, stat)
            if (stat == 0) then
               PNEWDT = 0.5
               return
            endif

            ! Calculate the deformation gradient at the element centroid
            !  at the the beginning and end of the increment for use in 
            !  the `F-bar' method
            !
            Fc_tau = Iden
            Fc_t = Iden
            do i = 1, nDim
               do j = 1, nDim
                  do k = 1, nNode
                     ! F at the end of increment
                     Fc_tau(i, j) = Fc_tau(i, j) + dsh0(k, j) * u(k, i)
                     ! F at the beginning of increment
                     Fc_t(i, j) = Fc_t(i, j) + dsh0(k, j) * uOld(k, i)
               end do
               end do
            end do
            ! 
            call mdet(Fc_tau, detFc_tau)
            call mdet(Fc_t, detFc_t)
            !
            ! With the deformation gradient known at the element centroid
            !  we are now able to implement the `F-bar' method later
            !
            !----------------------------------------------------------------
            !----------------------------------------------------------------
            ! Begin the loop over integration points
            !
            ! Obtain integration point local coordinates and weights
            !
            if (nIntt == 1) then
               call xint3D1pt(xi, w, nInttPt) ! 1-pt integration, nIntt=1 above
            elseif (nIntt == 8) then
               call xint3D8pt(xi, w, nInttPt) ! 8-pt integration, nIntt=8 above
            else
               write(*, *) 'Invalid number of int points, nIntt=', nIntt
               call exit
            endif

            nshr = nDim
            ntens = nDim + nshr
            ! Loop over integration points
            !
            jj = 0 ! jj is used for tracking the state variables
            !
            ! STATEV = [phi, detf, [cauchy_stresses]]
            do intpt = 1, nInttPt

               ! Obtain state variables from previous increment
               !
               if ((KINC <= 1) .and. (KSTEP == 1)) then
               ! this is the first increment, of the first step
               !  give initial conditions (or just anything)
               ! statev = 0.9999d0 !initial determinant of the deformation gradient
               prev_statev(1) = phi0
               prev_statev(2) = one
               statev(1) = phi0
               statev(2) = one
               phi_t  = phi0
               else
               ! this is not the first increment, read old values
               statev = SVARS(1 + jj : nsdv + jj)
               prev_statev = SVARS(1 + jj : nsdv + jj)
               phi_t  = svars(1+jj)
               endif

               ! Obtain shape functions and their local gradients
               !
               if (nNode == 8) then
               call calcShape3DLinear(nInttPt, xi, intpt, sh, dshxi)
               else
               write(*, *) 'Incorrect number of nodes: nNode.ne.8'
               call exit
               endif

               ! Map shape functions from local to global reference coordinate system
               !
               call mapShape3D(nNode, dshxi, coords, dsh, detMapJ, stat)
               if (stat == 0) then
               PNEWDT = 0.5
               return
               endif

               ! Map shape functions from local to global current coordinate system
               !
               call mapShape3D(nNode, dshxi, coordsC, dshC, detMapJC, stat)
               if (stat == 0) then
               PNEWDT = 0.5
               return
               endif

               ! Obtain the chemical potential and its derivative's at
               !  this intPt at the begining and end of the incrment
               !
               mu_tau = zero
               mu_t = zero
               dMUdt = zero
               dMUdX = zero
               do k=1,nNode
                  mu_tau = mu_tau + muNew(k)*sh(k)
                  mu_t   = mu_t + muOld(k)*sh(k)
                  do i=1,nDim
                     dMUdX(i,1) = dMUdX(i,1) + muNew(k)*dshC(k,i)
                  enddo
               enddo
               dMUdt = (mu_tau - mu_t)/dtime


               ! Obtain, and modify the deformation gradient at this integration
               !  point.  Modify the deformation gradient for use in the `F-bar'
               !  method.  Also, take care of plane-strain or axisymmetric
               !
               F_tau = Iden
               F_t = Iden
               do i = 1, nDim
               do j = 1, nDim
                  do k = 1, nNode
                  F_tau(i, j) = F_tau(i, j) + dsh(k, j) * u(k, i)
                  F_t(i, j) = F_t(i, j) + dsh(k, j) * uOld(k, i)
                  end do
               end do
               end do
               !
               ! Modify the deformation gradient for the `F-bar' method
               !  only when using the 8 node fully integrated linear
               !  element, do not use the `F-bar' method for any other element
               !
               if ((nNode == 8) .and. (nIntt == 8)) then
                  call mdet(F_tau, detF_tau)
                  call mdet(F_t, detF_t)
                  F_tau = ((detFc_tau / detF_tau)**(1.0d0 / 3.0d0)) * F_tau
                  !!! ORIGINAL:
                  ! F_t = ((detFc_tau / detF_tau)**(1.0d0 / 3.0d0)) * F_t
                  !!! CHANGED TO USE DETFC_T
                  F_t = ((detFc_t / detF_t)**(1.0d0 / 3.0d0)) * F_t
               endif
               call mdet(F_tau, detF)

               !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               !
               ! Perform the time integration at this integ. point to compute
               !  all the specific forms and parameters needed for the solution
               !
               call material(T_tau, statev, SpTanMod, &
                       F_t, F_tau, detF_tau, &
                       TIME, DTIME, PREDEF, &
                       nDim, nshr, ntens, nsdv, PROPS, NPROPS, coords, PNEWDT, &
                       JELEM, intpt, KSTEP, KINC,MU_TAU,PHI_T,THETA0,PHI_TAU,DPDT, &
                       DPHIDMU,DPHIDOTDMU,MFLUID,DMDMU,DMDJ,VMOL,SPUCMOD,SPCUMODFAC)
               !
               !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               ! Previous determinant of the deformation gradient
               write(*, *) 'detF_t=', prev_statev(2)
               ! Current determinant of the deformation gradient
               write(*, *) 'detF_tau=', statev(2)
               ! Save the state variables at this integ point
               !  at the end of the increment
               !
               SVARS(1 + jj : nsdv + jj) = statev
               jj = jj + nlSdv 
         ! setup for the next intPt      
         ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         !
         globalSdv(jelem, intPt, 1:nsdv) = statev

         ! Time stepping algorithm based on the constitutive response
         phiLmt = 0.005d0
         phi_tau = statev(1)
         phi_t = prev_statev(1)
         umeror = abs((phi_tau - phi_t)/phiLmt)
         ! write(*, *) 'umeror=', umeror
         if (umeror <= 0.5d0) then
            pnewdt = 1.5d0
         elseif (umeror > 0.5d0 .and. umeror <= 0.8d0) then
            pnewdt = 1.25d0
         elseif (umeror > 0.8d0 .and. umeror <= 1.25d0) then
            pnewdt = 0.75d0
         else
            pnewdt = 0.5d0
         endif

         ! Compute/update the displacement residual vector
         Smat(1, 1) = T_tau(1, 1)
         Smat(2, 1) = T_tau(2, 2)
         Smat(3, 1) = T_tau(3, 3)
         Smat(4, 1) = T_tau(1, 2)
         Smat(5, 1) = T_tau(2, 3)
         Smat(6, 1) = T_tau(1, 3)

         Bmat = 0.0d0
         do kk = 1, nNode
            Bmat(1, 1 + nDim * (kk - 1)) = dshC(kk, 1)
            Bmat(2, 2 + nDim * (kk - 1)) = dshC(kk, 2)
            Bmat(3, 3 + nDim * (kk - 1)) = dshC(kk, 3)
            Bmat(4, 1 + nDim * (kk - 1)) = dshC(kk, 2)
            Bmat(4, 2 + nDim * (kk - 1)) = dshC(kk, 1)
            Bmat(5, 2 + nDim * (kk - 1)) = dshC(kk, 3)
            Bmat(5, 3 + nDim * (kk - 1)) = dshC(kk, 2)
            Bmat(6, 1 + nDim * (kk - 1)) = dshC(kk, 3)
            Bmat(6, 3 + nDim * (kk - 1)) = dshC(kk, 1)
         end do

         BodyForceRes = 0.0d0
         do kk = 1, nNode
            BodyForceRes(1 + nDim * (kk - 1), 1) = sh(kk) * body(1)
            BodyForceRes(2 + nDim * (kk - 1), 1) = sh(kk) * body(2)
            BodyForceRes(3 + nDim * (kk - 1), 1) = sh(kk) * body(3)
         end do

         Ru = Ru + detmapJC * w(intpt) * &
            (-matmul(transpose(Bmat), Smat) + BodyForceRes)      
         !
         ! Compute/update the chemical potential residual vector
         !
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo
         !
         ResFac = (dPdt)/(detF*Vmol*phi_tau*phi_tau)
         !
         Rc = Rc + detmapJC*w(intpt)*(transpose(Nvec)*ResFac - Mfluid*matmul(dshC,dMUdX))
         
         
         !
         ! Compute/update the displacement tangent matrix
         !
            Gmat = 0.0d0
          do kk = 1, nNode
            Gmat(1, 1 + nDim * (kk - 1)) = dshC(kk, 1)
            Gmat(2, 2 + nDim * (kk - 1)) = dshC(kk, 1)
            Gmat(3, 3 + nDim * (kk - 1)) = dshC(kk, 1)
            Gmat(4, 1 + nDim * (kk - 1)) = dshC(kk, 2)
            Gmat(5, 2 + nDim * (kk - 1)) = dshC(kk, 2)
            Gmat(6, 3 + nDim * (kk - 1)) = dshC(kk, 2)
            Gmat(7, 1 + nDim * (kk - 1)) = dshC(kk, 3)
            Gmat(8, 2 + nDim * (kk - 1)) = dshC(kk, 3)
            Gmat(9, 3 + nDim * (kk - 1)) = dshC(kk, 3)
          end do

          G0mat = 0.0d0
          do kk = 1, nNode
            G0mat(1, 1 + nDim * (kk - 1)) = dshC0(kk, 1)
            G0mat(2, 2 + nDim * (kk - 1)) = dshC0(kk, 1)
            G0mat(3, 3 + nDim * (kk - 1)) = dshC0(kk, 1)
            G0mat(4, 1 + nDim * (kk - 1)) = dshC0(kk, 2)
            G0mat(5, 2 + nDim * (kk - 1)) = dshC0(kk, 2)
            G0mat(6, 3 + nDim * (kk - 1)) = dshC0(kk, 2)
            G0mat(7, 1 + nDim * (kk - 1)) = dshC0(kk, 3)
            G0mat(8, 2 + nDim * (kk - 1)) = dshC0(kk, 3)
            G0mat(9, 3 + nDim * (kk - 1)) = dshC0(kk, 3)
          end do

          Amat = 0.0d0
          Amat(1, 1) = SpTanMod(1, 1, 1, 1)
          Amat(1, 2) = SpTanMod(1, 1, 2, 1)
          Amat(1, 3) = SpTanMod(1, 1, 3, 1)
          Amat(1, 4) = SpTanMod(1, 1, 1, 2)
          Amat(1, 5) = SpTanMod(1, 1, 2, 2)
          Amat(1, 6) = SpTanMod(1, 1, 3, 2)
          Amat(1, 7) = SpTanMod(1, 1, 1, 3)
          Amat(1, 8) = SpTanMod(1, 1, 2, 3)
          Amat(1, 9) = SpTanMod(1, 1, 3, 3)
          Amat(2, 1) = SpTanMod(2, 1, 1, 1)
          Amat(2, 2) = SpTanMod(2, 1, 2, 1)
          Amat(2, 3) = SpTanMod(2, 1, 3, 1)
          Amat(2, 4) = SpTanMod(2, 1, 1, 2)
          Amat(2, 5) = SpTanMod(2, 1, 2, 2)
          Amat(2, 6) = SpTanMod(2, 1, 3, 2)
          Amat(2, 7) = SpTanMod(2, 1, 1, 3)
          Amat(2, 8) = SpTanMod(2, 1, 2, 3)
          Amat(2, 9) = SpTanMod(2, 1, 3, 3)
          Amat(3, 1) = SpTanMod(3, 1, 1, 1)
          Amat(3, 2) = SpTanMod(3, 1, 2, 1)
          Amat(3, 3) = SpTanMod(3, 1, 3, 1)
          Amat(3, 4) = SpTanMod(3, 1, 1, 2)
          Amat(3, 5) = SpTanMod(3, 1, 2, 2)
          Amat(3, 6) = SpTanMod(3, 1, 3, 2)
          Amat(3, 7) = SpTanMod(3, 1, 1, 3)
          Amat(3, 8) = SpTanMod(3, 1, 2, 3)
          Amat(3, 9) = SpTanMod(3, 1, 3, 3)
          Amat(4, 1) = SpTanMod(1, 2, 1, 1)
          Amat(4, 2) = SpTanMod(1, 2, 2, 1)
          Amat(4, 3) = SpTanMod(1, 2, 3, 1)
          Amat(4, 4) = SpTanMod(1, 2, 1, 2)
          Amat(4, 5) = SpTanMod(1, 2, 2, 2)
          Amat(4, 6) = SpTanMod(1, 2, 3, 2)
          Amat(4, 7) = SpTanMod(1, 2, 1, 3)
          Amat(4, 8) = SpTanMod(1, 2, 2, 3)
          Amat(4, 9) = SpTanMod(1, 2, 3, 3)
          Amat(5, 1) = SpTanMod(2, 2, 1, 1)
          Amat(5, 2) = SpTanMod(2, 2, 2, 1)
          Amat(5, 3) = SpTanMod(2, 2, 3, 1)
          Amat(5, 4) = SpTanMod(2, 2, 1, 2)
          Amat(5, 5) = SpTanMod(2, 2, 2, 2)
          Amat(5, 6) = SpTanMod(2, 2, 3, 2)
          Amat(5, 7) = SpTanMod(2, 2, 1, 3)
          Amat(5, 8) = SpTanMod(2, 2, 2, 3)
          Amat(5, 9) = SpTanMod(2, 2, 3, 3)
          Amat(6, 1) = SpTanMod(3, 2, 1, 1)
          Amat(6, 2) = SpTanMod(3, 2, 2, 1)
          Amat(6, 3) = SpTanMod(3, 2, 3, 1)
          Amat(6, 4) = SpTanMod(3, 2, 1, 2)
          Amat(6, 5) = SpTanMod(3, 2, 2, 2)
          Amat(6, 6) = SpTanMod(3, 2, 3, 2)
          Amat(6, 7) = SpTanMod(3, 2, 1, 3)
          Amat(6, 8) = SpTanMod(3, 2, 2, 3)
          Amat(6, 9) = SpTanMod(3, 2, 3, 3)
          Amat(7, 1) = SpTanMod(1, 3, 1, 1)
          Amat(7, 2) = SpTanMod(1, 3, 2, 1)
          Amat(7, 3) = SpTanMod(1, 3, 3, 1)
          Amat(7, 4) = SpTanMod(1, 3, 1, 2)
          Amat(7, 5) = SpTanMod(1, 3, 2, 2)
          Amat(7, 6) = SpTanMod(1, 3, 3, 2)
          Amat(7, 7) = SpTanMod(1, 3, 1, 3)
          Amat(7, 8) = SpTanMod(1, 3, 2, 3)
          Amat(7, 9) = SpTanMod(1, 3, 3, 3)
          Amat(8, 1) = SpTanMod(2, 3, 1, 1)
          Amat(8, 2) = SpTanMod(2, 3, 2, 1)
          Amat(8, 3) = SpTanMod(2, 3, 3, 1)
          Amat(8, 4) = SpTanMod(2, 3, 1, 2)
          Amat(8, 5) = SpTanMod(2, 3, 2, 2)
          Amat(8, 6) = SpTanMod(2, 3, 3, 2)
          Amat(8, 7) = SpTanMod(2, 3, 1, 3)
          Amat(8, 8) = SpTanMod(2, 3, 2, 3)
          Amat(8, 9) = SpTanMod(2, 3, 3, 3)
          Amat(9, 1) = SpTanMod(3, 3, 1, 1)
          Amat(9, 2) = SpTanMod(3, 3, 2, 1)
          Amat(9, 3) = SpTanMod(3, 3, 3, 1)
          Amat(9, 4) = SpTanMod(3, 3, 1, 2)
          Amat(9, 5) = SpTanMod(3, 3, 2, 2)
          Amat(9, 6) = SpTanMod(3, 3, 3, 2)
          Amat(9, 7) = SpTanMod(3, 3, 1, 3)
          Amat(9, 8) = SpTanMod(3, 3, 2, 3)
          Amat(9, 9) = SpTanMod(3, 3, 3, 3)

          Qmat = 0.0d0
          Qmat(1, 1) = (1.0d0 / 3.0d0) * (Amat(1, 1) + Amat(1, 5) + Amat(1, 9)) - (2.0d0 / 3.0d0) * T_tau(1, 1)
          Qmat(2, 1) = (1.0d0 / 3.0d0) * (Amat(2, 1) + Amat(2, 5) + Amat(2, 9)) - (2.0d0 / 3.0d0) * T_tau(2, 1)
          Qmat(3, 1) = (1.0d0 / 3.0d0) * (Amat(3, 1) + Amat(3, 5) + Amat(3, 9)) - (2.0d0 / 3.0d0) * T_tau(3, 1)
          Qmat(4, 1) = (1.0d0 / 3.0d0) * (Amat(4, 1) + Amat(4, 5) + Amat(4, 9)) - (2.0d0 / 3.0d0) * T_tau(1, 2)
          Qmat(5, 1) = (1.0d0 / 3.0d0) * (Amat(5, 1) + Amat(5, 5) + Amat(5, 9)) - (2.0d0 / 3.0d0) * T_tau(2, 2)
          Qmat(6, 1) = (1.0d0 / 3.0d0) * (Amat(6, 1) + Amat(6, 5) + Amat(6, 9)) - (2.0d0 / 3.0d0) * T_tau(3, 2)
          Qmat(7, 1) = (1.0d0 / 3.0d0) * (Amat(7, 1) + Amat(7, 5) + Amat(7, 9)) - (2.0d0 / 3.0d0) * T_tau(1, 3)
          Qmat(8, 1) = (1.0d0 / 3.0d0) * (Amat(8, 1) + Amat(8, 5) + Amat(8, 9)) - (2.0d0 / 3.0d0) * T_tau(2, 3)
          Qmat(9, 1) = (1.0d0 / 3.0d0) * (Amat(9, 1) + Amat(9, 5) + Amat(9, 9)) - (2.0d0 / 3.0d0) * T_tau(3, 3)
          Qmat(1, 5) = Qmat(1, 1)
          Qmat(2, 5) = Qmat(2, 1)
          Qmat(3, 5) = Qmat(3, 1)
          Qmat(4, 5) = Qmat(4, 1)
          Qmat(5, 5) = Qmat(5, 1)
          Qmat(6, 5) = Qmat(6, 1)
          Qmat(7, 5) = Qmat(7, 1)
          Qmat(8, 5) = Qmat(8, 1)
          Qmat(9, 5) = Qmat(9, 1)
          Qmat(1, 9) = Qmat(1, 1)
          Qmat(2, 9) = Qmat(2, 1)
          Qmat(3, 9) = Qmat(3, 1)
          Qmat(4, 9) = Qmat(4, 1)
          Qmat(5, 9) = Qmat(5, 1)
          Qmat(6, 9) = Qmat(6, 1)
          Qmat(7, 9) = Qmat(7, 1)
          Qmat(8, 9) = Qmat(8, 1)
          Qmat(9, 9) = Qmat(9, 1)

          if ((nNode == 8) .and. (nIntt == 8)) then
            ! This is the tangent using the F-bar method with the 8 node fully integrated linear element
            ! Correction required to the tangent matrix for the `F-bar' method
            Kuu = Kuu + detMapJC * w(intpt) * &
                 (matmul(matmul(transpose(Gmat), Amat), Gmat) + &
                 matmul(transpose(Gmat), matmul(Qmat, (G0mat - Gmat))))
          else
            ! This is the tangent NOT using the F-bar method with all other elements
            Kuu = Kuu + detMapJC * w(intpt) * &
                 (matmul(matmul(transpose(Gmat), Amat), Gmat))
          end if

         !
         ! Compute/update the chemical potential tangent matrix
         !
         TanFac = (one/(detF*Vmol*phi_tau**two))* &
                  (two*(dPdt/phi_tau)*DphiDmu - DphidotDmu)
         !
         Kcc = Kcc + detmapJC*w(intPt)* &
                  (TanFac*matmul(transpose(Nvec),Nvec) &
                  + Mfluid*matmul(dshC,transpose(dshC)) &
                  + DmDmu*matmul(matmul(dshC,dMUdX),Nvec))

         ! Compute/update the chemical potential - displacement tangent matrix
         !  The F-bar method will have some effect, however we neglect that here.
         !
         SpCUMod = zero
         do i=1,nDim
            do k=1,nDim
               do l=1,nDim
                  SpCUMod(i,k,l) = SpCUMod(i,k,l) &
                                 + dMUdX(k,1)*SpCUModFac(i,l)
               enddo
            enddo
         enddo
         !
         AmatCU = zero
         AmatCU(1,1) = SpCUMod(1,1,1)
         AmatCU(1,2) = SpCUMod(1,2,1)
         AmatCU(1,3) = SpCUMod(1,3,1)
         AmatCU(1,4) = SpCUMod(1,1,2)
         AmatCU(1,5) = SpCUMod(1,2,2)
         AmatCU(1,6) = SpCUMod(1,3,2)
         AmatCU(1,7) = SpCUMod(1,1,3)
         AmatCU(1,8) = SpCUMod(1,2,3)
         AmatCU(1,9) = SpCUMod(1,3,3)
         AmatCU(2,1) = SpCUMod(2,1,1)
         AmatCU(2,2) = SpCUMod(2,2,1)
         AmatCU(2,3) = SpCUMod(2,3,1)
         AmatCU(2,4) = SpCUMod(2,1,2)
         AmatCU(2,5) = SpCUMod(2,2,2)
         AmatCU(2,6) = SpCUMod(2,3,2)
         AmatCU(2,7) = SpCUMod(2,1,3)
         AmatCU(2,8) = SpCUMod(2,2,3)
         AmatCU(2,9) = SpCUMod(2,3,3)
         AmatCU(3,1) = SpCUMod(3,1,1)
         AmatCU(3,2) = SpCUMod(3,2,1)
         AmatCU(3,3) = SpCUMod(3,3,1)
         AmatCU(3,4) = SpCUMod(3,1,2)
         AmatCU(3,5) = SpCUMod(3,2,2)
         AmatCU(3,6) = SpCUMod(3,3,2)
         AmatCU(3,7) = SpCUMod(3,1,3)
         AmatCU(3,8) = SpCUMod(3,2,3)
         AmatCU(3,9) = SpCUMod(3,3,3)
         !
         Kcu = Kcu - detMapJC*w(intpt)* &
               (matmul(matmul(dshC,AmatCU),Gmat))


         ! Compute/update the displacement - chemical potential tangent matrix
         !  The F-bar method will have some effect, however we neglect that here.
         !
         AmatUC = zero
         AmatUC(1,1) = SpUCMod(1,1)
         AmatUC(2,1) = SpUCMod(2,2)
         AmatUC(3,1) = SpUCMod(3,3)
         AmatUC(4,1) = SpUCMod(1,2)
         AmatUC(5,1) = SpUCMod(2,3)
         AmatUC(6,1) = SpUCMod(1,3)
         !
         Kuc = Kuc + detMapJC*w(intpt)* &
               (matmul(matmul(transpose(Bmat),AmatUC),Nvec))

         
      end do
      !
      ! End the loop over integration points
      !----------------------------------------------------------------

          

      !----------------------------------------------------------------
      ! Start loop over surface flux terms
      !
      if (ndload > 0) then
         !
         ! loop over faces and make proper modifications to
         !  residuals and tangents if needed
         !
         do i = 1, ndload
         !
         ! based on my convention the face which the flux
         !  acts on is the flux ``label''
         !
         face = jdltyp(i, 1)
         flux = adlmag(i, 1)
         !
         if ((face >= 1) .and. (face <= 6)) then
            !
            ! fluid flux applied
            !
            select case (face)
            case (1)
            faceFlag = 1
            case (2)
            faceFlag = 2
            case (3)
            faceFlag = 3
            case (4)
            faceFlag = 4
            case (5)
            faceFlag = 5
            case (6)
            faceFlag = 6
            end select
            !
            select case (nInttS)
            case (1)
            call xintSurf3D1pt(faceFlag, xLocal, yLocal, zLocal, wS)
            case (4)
            call xintSurf3D4pt(faceFlag, xLocal, yLocal, zLocal, wS)
            case default
            write(*, *) 'Invalid nInttS points, nInttS=', nInttS
            call exit
            end select
            !
            ! loop over integ points on this element face
            !
            do ii = 1, nInttS
            
            ! Compute shape functions, derivatives, and the 
            !  mapping jacobian (dA)
            !
            call computeSurf3D(xLocal(ii), yLocal(ii), zLocal(ii), &
               faceFlag, coordsC, sh, dA)
            !
            ! Modify the chemical potential residual, loop over nodes
            !
            do n = 1, nNode
               Rc(n, 1) = Rc(n, 1) - wS(ii) * dA * sh(n) * flux
            end do 
            !
            ! No change to the tangent matrix
            !
            end do ! end loop over integ points
            !
         else
            write(*, *) 'Unknown face=', face
            call exit
         end if

         end do ! loop over ndload
      end if ! ndload.gt.0 or not
      !
      ! End loop over surface flux terms
      !----------------------------------------------------------------  
   !
      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      !
      
      call AssembleElement(nDim, nNode, nDofEl, &
         Ru,Rc,Kuu,Kuc,Kcu,Kcc, &
         rhs, amatrx)
   !      write(*,*) rhs(:,1)
   !      write(*,*) amatrx     
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------
      return 
      end subroutine U3D8