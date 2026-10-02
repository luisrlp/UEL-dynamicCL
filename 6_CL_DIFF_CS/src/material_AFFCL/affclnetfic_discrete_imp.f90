! SUBROUTINE affclnetfic_discrete(sfic,cfic,f,filprops,affprops,  &
!           efi,noel,det,prefdir,ndi) ! (original)

SUBROUTINE affclnetfic_discrete(sfic,cfic,f,unit2,filprops,affprops,  &
  efi,noel,det,prefdir,ndi,cb,dtime,chemprops, &
  thetaf, cbtau_tot, dPK2ficdcb, dcbdc, pnewdt)



!>    AFFINE NETWORK: 'FICTICIOUS' CAUCHY STRESS AND ELASTICITY TENSOR
!> DISCRETE ANGULAR INTEGRATION SCHEME (icosahedron)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: dPK2ficdcb(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: dcbdc(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pnewdt
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit2(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: filprops(10)
DOUBLE PRECISION, INTENT(IN)             :: affprops(5)
DOUBLE PRECISION, INTENT(IN)             :: chemprops(10)
DOUBLE PRECISION, INTENT(IN OUT)         :: efi
INTEGER, INTENT(IN OUT)                  :: noel
DOUBLE PRECISION, INTENT(IN OUT)         :: det

DOUBLE PRECISION, INTENT(IN)             :: dtime
DOUBLE PRECISION, INTENT(IN)             :: thetaf
DOUBLE PRECISION, INTENT(OUT)            :: cbtau_tot
DOUBLE PRECISION, INTENT(IN OUT)         :: cb(ndir)

INTEGER :: i1,j1,k1,l1,m1, im1, isub, n_sub
INTEGER, PARAMETER :: nargs = 19
DOUBLE PRECISION :: args(nargs)
DOUBLE PRECISION :: sfilfic(ndi,ndi), cfilfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,lambdai,dwi,ddwi,rwi,lambdaic
DOUBLE PRECISION :: l,Lp,r0f,r0,mu0str,b0,beta,lambda0,lambda0f,rho,n,fi,ffi,aratio
DOUBLE PRECISION :: r0c,etac,lambdaif
DOUBLE PRECISION :: bdisp,ang, frac(4)
DOUBLE PRECISION :: prefdir(nelem,4), pd(3),lambda_pref,prefdir0(3)
DOUBLE PRECISION :: dx,dxc,kb,theta,na
DOUBLE PRECISION :: cactin, Mactin, rhoactin, cbmax, cfmax, chi, D, MU0, VMOL, Koff0, Keq, Kcatch0
DOUBLE PRECISION :: cbt_i, cbtau_i, thetab_i, Kon, Koff_i, R_i, dtime_sub, cb_sub
INTEGER :: iter
DOUBLE PRECISION :: cb_new, Res, cb_pert, dcb, r0f_p, l_p, r0_p
DOUBLE PRECISION :: DfDcb, DdwDcb, Dr0Dcb, auxdwdcb, dHdlambda, dHdcb, dRiDcb, auxchem, dKoffDf_i
DOUBLE PRECISION :: dHdcb_i, dHdl_i, sens_i, dHdcb_min, kh_i, gerr_i, cbx_i, Hx_i, dHdcbx_i, dHdlx_i
LOGICAL :: conv_i, atb_i, ok_i, hard_i
INTEGER :: nhard, nacc, nsub_max
! Kinetics sub-stepping (per direction). The whole increment is tried first (n_sub = 1); it is
! repeated from the start-of-increment state with n_sub = 2, 4, ... <= MAX_SUBSTEPS sub-steps of
! dtime/n_sub if, in any sub-step, the solve does not converge, dH/dcb at the root < DH_MIN
! (uniqueness of the Backward-Euler root being lost), or the estimated Backward-Euler error
! exceeds DCB_MAX*max(min(cb_start, cb_end), CB_SCALE*cbmax) (accuracy, relative to the smaller
! of the two values). Error estimate (exact order for linear kinetics):
! err = |cb change| * min(kh/2, 1/kh), with kh = dt_sub*|dR/dcb| = |dH/dcb - 1| the sub-step
! length over the kinetic time scale, evaluated at the start and at the end of the sub-step
! (larger factor used): slow kinetics (kh << 1) err ~ |dcb|*kh/2; fast kinetics (kh >> 1,
! quasi-equilibrium, e.g. complete unbinding) err ~ |dcb|/kh -> 0.
DOUBLE PRECISION, PARAMETER :: DH_MIN = 0.1d0, DCB_MAX = 0.1d0, CB_SCALE = 1.0d-3
INTEGER, PARAMETER :: MAX_SUBSTEPS = 64
DOUBLE PRECISION :: dPK2filficdcb(ndi,ndi), pfdlambdadcfil(ndi,ndi), cfilficchem(ndi,ndi,ndi,ndi)

! INTEGRATION SCHEME
  integer, parameter :: nfacedir = 2
  integer ( kind = 4 ) ifacedir
  integer :: f3_start(nfacedir), f3_end(nfacedir), f2_start(nfacedir)
  integer, dimension(3, nfacedir) :: off_a, off_b, off_c
  integer ( kind = 4 ) a, b, c
  real ( kind = 8 ) a_xyz(3), b_xyz(3), c_xyz(3)
  real ( kind = 8 ) a2_xyz(3), b2_xyz(3), c2_xyz(3)
  real ( kind = 8 ) area_total, ai !area of triangle i
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: edge_point
  integer ( kind = 4 ) f1, f2, f3
  integer ( kind = 4 ) face, face_num, face_order_max, node_num, edge_num, point_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: face_order
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: face_point
  real ( kind = 8 ) node_xyz(3)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_coord
  real ( kind = 8 ) rr, aa, v



!  Size the icosahedron.
!
  call icos_size ( point_num, edge_num, face_num, face_order_max )
!
!  Set the icosahedron.
!
  allocate ( point_coord(1:3,1:point_num) )
  allocate ( edge_point(1:2,1:edge_num) )
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )

  call icos_shape ( point_num, edge_num, face_num, face_order_max, &
    point_coord, edge_point, face_order, face_point )
!
!  Set aux variables for the integration scheme 
!
f3_start(1) = 1; f3_end(1) = 3 * factor - 2
f2_start(1) = 1
f3_start(2) = 2; f3_end(2) = 3 * factor - 4
f2_start(2) = 2
off_a(:,1) = [2, -1, -1];  off_b(:,1) = [-1, 2, -1];  off_c(:,1) = [-1, -1, 2]
off_a(:,2) = [-2, 1, 1];   off_b(:,2) = [1, -2, 1];   off_c(:,2) = [1, 1, -2]
!
!  Initialize the integral data.
!
  rr = 0.0D+00
  area_total = 0.0D+00
  node_num = 0

!! initialize the model data
  !     FILAMENT
  aratio   = filprops(1)
  r0c      = filprops(2)
  etac     = filprops(3)
  mu0str   = filprops(4)
  beta     = filprops(5)
  Lp       = filprops(6)
  theta    = filprops(7)
  dx       = filprops(8)
  kb       = filprops(9)
  NA       = filprops(10)
  b0       = Lp * theta * kb
  !     NETWORK
  bdisp    = affprops(1)
  lambda0  = affprops(2)                                                                                                                                                           
  cactin   = affprops(3)                                                                                              
  Mactin   = affprops(4)                                                                                            
  rhoactin = affprops(5)
  !     CHEMICAL
  cbmax    = chemprops(1)
  cfmax    = chemprops(2)
  CHI      = chemprops(3)
  D        = chemprops(4)
  MU0      = chemprops(5)
  VMOL     = chemprops(6)
  Koff0    = chemprops(7)
  Keq      = chemprops(8)
  Kcatch0  = chemprops(9)
  dxc      = chemprops(10)

  ! aux=n*(det**(-one))
  cfic=zero
  sfic=zero

  ! rho=one
  r0=r0f+r0c

  aa = zero

  cbtau_tot = zero
  thetab_i = zero
  dPK2ficdcb=zero
  dcbdc=zero
  ! Kinetics diagnostics over all directions (unresolved directions, max sub-steps, min dH/dcb)
  nhard = 0
  nacc = 0
  nsub_max = 1
  dHdcb_min = huge(one)
!----------------------------------------------------------------------
  
  ! preferred direction measures (macroscale measures)
  ! prefdir0=prefdir(noel,2:4)
  ! Currently assuming all elements have the same preferential direction
  prefdir0=prefdir(1,2:4)
  !calculate preferred direction in the deformed configuration
  CALL deffil(lambda_pref,pd,prefdir0,f,ndi)
  !update preferential direction - deformed configuration
  pd=pd/dsqrt(dot_product(pd,pd))

!  Pick a face of the icosahedron, and identify its vertices as A, B, C.
!
! Integrate only one hemisphere of the icosahedron (faces 1 to 10) 
! Remember to multiply each direction's contribution by 2 to account for the other hemisphere
do face = 1, face_num/2
!
    a = face_point(1,face)
    b = face_point(2,face)
    c = face_point(3,face)
!
    a_xyz(1:3) = point_coord(1:3,a)
    b_xyz(1:3) = point_coord(1:3,b)
    c_xyz(1:3) = point_coord(1:3,c)
!
!  Some subtriangles will have the same direction as the face.
!  Generate each in turn, by determining the barycentric coordinates
!  of the centroid (F1,F2,F3), from which we can also work out the barycentric
!  coordinates of the vertices of the subtriangle.
!
  do ifacedir = 1, nfacedir
    do f3 = f3_start(ifacedir), f3_end(ifacedir), 3
      do f2 = f2_start(ifacedir), 3 * factor - f3 - ifacedir, 3

        f1 = 3 * factor - f3 - f2

        node_num = node_num + 1

        call sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
          node_xyz )

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + off_a(1,ifacedir), f2 + off_a(2,ifacedir), f3 + off_a(3,ifacedir), a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + off_b(1,ifacedir), f2 + off_b(2,ifacedir), f3 + off_b(3,ifacedir), b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + off_c(1,ifacedir), f2 + off_c(2,ifacedir), f3 + off_c(3,ifacedir), c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, ai )
        
        !direction of the sphere triangle barycenter - direction i
        mf0i=node_xyz
        ! write(*,*) 'mf0i =', mf0i
        CALL deffil(lambdai,mfi,mf0i,f,ndi)
        
        CALL bangle(ang,f,mfi,noel,pd,ndi)
        
        CALL density(rho,ang,bdisp,efi)
        
        fi = zero
        DfDcb = zero
        DdwDcb = zero
        dPK2filficdcb=zero
          ! ================= KINETICS (IMPLICIT NEWTON-RAPHSON) =================
        cbt_i = MAX(cb(node_num), 1.0d-10)
        cbtau_i = cbt_i  ! Initial guess is the old state
        kon = (Koff0 + Kcatch0) * Keq * exp(CHI * (1.0d0 - 2.0d0 * thetaf))

        args(1)  = lambdai
        args(2)  = lambda0
        args(3)  = aratio
        args(4)  = etac
        args(5)  = mu0str
        args(6)  = beta
        args(7)  = b0
        args(8)  = r0c
        args(9)  = cbmax
        args(10) = cfmax
        args(11) = dx / (kb * theta)
        args(12) = dtime
        args(13) = kon
        args(14) = Koff0
        args(15) = thetaf
        args(16) = cbt_i
        args(17) = det
        args(18) = Kcatch0
        args(19) = dxc / (kb * theta)

        ! Backward-Euler kinetics with automatic sub-stepping (lambdai and thetaf held at their
        ! end-of-increment values, so every sub-step is Backward Euler; n_sub = 1 is the plain step).
        ! Consistent tangent: sub-step k solves H_k(cb_k, cb_k-1, lambda) = 0, hence the sensitivity
        !   sens = dcb/dlambda_i:  sens_k = (sens_k-1 - dH_k/dlambda) / (dH_k/dcb),  sens_0 = 0
        ! (for n_sub = 1: sens = -dHdlambda/dHdcb).
        n_sub = 1
        DO
          cb_sub = cbt_i
          sens_i = zero
          ok_i = .TRUE.
          hard_i = .FALSE.
          args(12) = dtime / DBLE(n_sub)
          DO isub = 1, n_sub
            args(16) = cb_sub
            CALL solveKinetics(cb_new, args, nargs, cb_sub, conv_i, atb_i, dHdcb_i, dHdl_i)
            IF (dHdcb_i /= zero) sens_i = (sens_i - dHdl_i) / dHdcb_i
            ! Hard failures (no converged root, or uniqueness being lost)
            IF ((.NOT. conv_i) .OR. (dHdcb_i < DH_MIN)) THEN
              ok_i = .FALSE.
              hard_i = .TRUE.
            END IF
            ! Error factor g(kh) = min(kh/2, 1/kh) at the end (root) and at the start of the
            ! sub-step; the larger one is used (strongly nonlinear kinetics, e.g. fast transients)
            kh_i = ABS(dHdcb_i - one)
            gerr_i = MIN(half * kh_i, one / MAX(kh_i, 1.0d-300))
            cbx_i = cb_sub
            CALL kineticsFunc(cbx_i, Hx_i, dHdcbx_i, dHdlx_i, args, nargs)
            kh_i = ABS(dHdcbx_i - one)
            gerr_i = MAX(gerr_i, MIN(half * kh_i, one / MAX(kh_i, 1.0d-300)))
            IF (ABS(cb_new - cb_sub) * gerr_i &
                > DCB_MAX * MAX(MIN(cb_sub, cb_new), CB_SCALE * cbmax)) ok_i = .FALSE.
            dHdcb_min = MIN(dHdcb_min, dHdcb_i)
            cb_sub = cb_new
            ! Abandon this attempt and refine (the last attempt is always completed)
            IF ((.NOT. ok_i) .AND. (n_sub < MAX_SUBSTEPS)) EXIT
          END DO
          IF (ok_i .OR. (n_sub >= MAX_SUBSTEPS)) EXIT
          n_sub = 2 * n_sub
        END DO
        ! Not resolved with MAX_SUBSTEPS: hard failure (-> PNEWDT cut-back) or accuracy only (warning)
        IF (.NOT. ok_i) THEN
          IF (hard_i) THEN
            nhard = nhard + 1
          ELSE
            nacc = nacc + 1
          END IF
        END IF
        nsub_max = MAX(nsub_max, n_sub)
        cbtau_i = cb_sub

        cb(node_num) = cbtau_i
        thetab_i = cbtau_i / cbmax

        ! FINAL STATE UPDATE
        r0f = 1.6 * (1.0d3 * cb(node_num))**(-two/5.0d0)
        l = aratio * r0f
        r0 = r0f + r0c
        n = l**(-1) * (cactin * NA * Mactin / rhoactin)
        aux = n * (det**(-one))

        IF((etac > zero).AND.(etac .LE. one)) THEN
          lambdaif = etac*(r0/r0f)*(lambdai-one)+one
          lambda0f = etac*(r0/r0f)*(lambda0-one)+one
        ELSE
          lambdaif = lambdai
          lambda0f = lambda0
        END IF

        fi = zero
        ! Loaded filament: total filament stretch (prestretch included) >= 1, as in kineticsFunc.
        ! With LAMBDA0 = 1 this is lambdai >= 1 (the reference state F = I keeps the filament tangent)
        IF(lambdaif*lambda0f .GE. one) THEN
          CALL fil(fi,ffi,dwi,ddwi,lambdai,lambdaif,lambda0,lambda0f,l,r0,r0f,mu0str,beta,b0,etac,cb(node_num),DfDcb)
          !CALL fil_inext(fi,dwi,ddwi,lambdai,lambdaif,lambda0,lambda0f,l,r0,r0f,beta,b0,etac,cb(node_num),DfDcb)
          CALL sigfilfic(sfilfic,rho,lambdai,dwi,mfi,ai,ndi)
          CALL csfilfic(cfilfic,rho,lambdai,dwi,ddwi,mfi,ai,ndi)
          CALL csfilficchem(cfilficchem,rho,lambdai,mfi,ai,ndi)
          ! Chemical contributions to the mechanical tangent (Kuu)
          !! Term 3 - volumetric
          ! kin_aux = 
          !! Term 3.1
          dr0dcb = - (two / 5.d0) * r0f / cb(node_num)
          ddwdcb = lambda0 * (r0 * dfdcb + fi * dr0dcb)
          auxdwdcb = two / 5.d0 * dwi / cb(node_num) + ddwdcb
          !!! Leverage sigfilfic to get dSfic/Dcb for current direction
          call sigfilfic(dPK2filficdcb,rho,lambdai,auxdwdcb,mf0i,ai,ndi)
          !! Term 3.2
          call pfdlambdadc(pfdlambdadcfil,rho,lambdai,unit2,mfi,ai,det,ndi)
          ! auxchem = -dcb/dlambda_i from the sub-step sensitivity recursion
          ! (n_sub = 1: auxchem = dHdlambda/dHdcb, as before)
          auxchem = - sens_i
          ! auxchem = aux * auxdwdcb * (dHdcb)**(-one) * dHdlambda
          ! write(*,*) 'auxchem =', auxchem
          ! write(*,*) 'auxdwdcb =', auxdwdcb
          ! write(*,*) 'aux =', aux
          ! write(*,*) 'auxs*cfilficchem / cfilfic', auxdwdcb*auxchem*cfilficchem / cfilfic
          DO j1=1,ndi
            DO k1=1,ndi
                sfic(j1,k1) = sfic(j1,k1) + aux*sfilfic(j1,k1)
                ! Kuu - 3.1
                dPK2ficdcb(j1,k1) = dPK2ficdcb(j1,k1) + n * dPK2filficdcb(j1,k1)
                ! Kuu - 3 - volumetric
                dcbdC(j1,k1) = dcbdC(j1,k1) + lambdai**(-one) * auxchem * pfdlambdadcfil(j1,k1)
                DO l1=1,ndi
                  DO m1=1,ndi
                    cfic(j1,k1,l1,m1) = cfic(j1,k1,l1,m1) + aux*cfilfic(j1,k1,l1,m1) &
                    - aux*auxdwdcb*auxchem*cfilficchem(j1,k1,l1,m1)
                  END DO
                END DO
            END DO
          END DO
        END IF
        
        ! Update total bound CL concentration for this integration point                                                       
        cbtau_tot = cbtau_tot + cb(node_num) * rho * ai
        ! ==============================================================   

        !v=dwi
        !rr = rr + ai * v
        !area_total = area_total + ai

      end do
    end do
  end do
  end do

  ! Directions whose kinetics were not resolved with MAX_SUBSTEPS sub-steps:
  !  - hard failure (no converged root or dH/dcb < DH_MIN): ask Abaqus to repeat the increment
  !    with a smaller time increment (dH/dcb = 1 + dt*(...) recovers as dt decreases);
  !  - accuracy check only: the result is kept (the per-sub-step error estimate is conservative),
  !    and a warning is written. PNEWDT < 1 would discard an otherwise converged increment.
  IF (nhard > 0) THEN
    pnewdt = MIN(pnewdt, half)
    write(*,*) 'WARNING: local kinetics failed (not converged or dH/dcb <', DH_MIN, ') with', MAX_SUBSTEPS, &
               'sub-steps in', nhard, 'directions. noel =', noel, ' min dH/dcb =', dHdcb_min, ' -> PNEWDT =', pnewdt
  END IF
  IF (nacc > 0) THEN
    write(*,*) 'WARNING: kinetics accuracy check not met with', MAX_SUBSTEPS, 'sub-steps in', nacc, &
               'directions (result kept). noel =', noel
  END IF
!
!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )

RETURN
END SUBROUTINE affclnetfic_discrete
