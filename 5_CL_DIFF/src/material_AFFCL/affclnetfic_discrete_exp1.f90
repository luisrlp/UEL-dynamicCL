! SUBROUTINE affclnetfic_discrete(sfic,cfic,f,filprops,affprops,  &
!           efi,noel,det,prefdir,ndi) ! (original)

SUBROUTINE affclnetfic_discrete(sfic,cfic,f,filprops,affprops,  &
  efi,noel,det,prefdir,ndi,cb,dtime,cabp,cfmax,cbmax,chi,Keq,Koff0, &
  thetaf, cb_tot_new)



!>    AFFINE NETWORK: 'FICTICIOUS' CAUCHY STRESS AND ELASTICITY TENSOR
!> DISCRETE ANGULAR INTEGRATION SCHEME (icosahedron)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: filprops(10)
DOUBLE PRECISION, INTENT(IN)             :: affprops(5)
DOUBLE PRECISION, INTENT(IN OUT)         :: efi
INTEGER, INTENT(IN OUT)                  :: noel
DOUBLE PRECISION, INTENT(IN OUT)         :: det

DOUBLE PRECISION, INTENT(IN)             :: dtime
DOUBLE PRECISION, INTENT(IN)             :: cabp
DOUBLE PRECISION, INTENT(IN)             :: cfmax
DOUBLE PRECISION, INTENT(IN)             :: cbmax
DOUBLE PRECISION, INTENT(IN)             :: CHI
DOUBLE PRECISION, INTENT(IN)             :: Keq
DOUBLE PRECISION, INTENT(IN)             :: Koff0
DOUBLE PRECISION, INTENT(IN)             :: thetaf
DOUBLE PRECISION, INTENT(OUT)            :: cb_tot_new
DOUBLE PRECISION, INTENT(IN OUT)         :: cb(ndir)

INTEGER :: i1,j1,k1,l1,m1, im1, isub, n_sub
DOUBLE PRECISION :: sfilfic(ndi,ndi), cfilfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,lambdai,dwi,ddwi,rwi,lambdaic
DOUBLE PRECISION :: l,Lp,r0f,r0,mu0str,b0,beta,lambda0,lambda0f,rho,n,fi,ffi,aratio
DOUBLE PRECISION :: r0c,etac,lambdaif,lambdaimax
DOUBLE PRECISION :: bdisp,ang, frac(4)
DOUBLE PRECISION :: prefdir(nelem,4), pd(3),lambda_pref,prefdir0(3)
DOUBLE PRECISION :: dx,kb,theta,na
DOUBLE PRECISION :: cactin, Mactin, rhoactin
DOUBLE PRECISION :: cb_i, thetab_i, Kon, Koff_i, R_i

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
  
    ! aux=n*(det**(-one))
    cfic=zero
    sfic=zero
  
    ! rho=one
    r0=r0f+r0c
  
    aa = zero
    lambdaimax=zero

    cb_tot_new = zero
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
        
        ! ================= DYNAMIC GEOMETRY =================
        cb_i = MAX(cb(node_num), 1.0d-8)
        r0f = 1.6 * (10.0d3 * cb_i)**(-two/5.0d0)
        l = aratio * r0f
        r0 = r0f + r0c
        n = l**(-1) * (cactin * NA * Mactin / rhoactin)
        aux = n * (det**(-one))
        ! ====================================================

        !direction of the sphere triangle barycenter - direction i
        mf0i=node_xyz
        CALL deffil(lambdai,mfi,mf0i,f,ndi)

        CALL bangle(ang,f,mfi,noel,pd,ndi)
  
        CALL density(rho,ang,bdisp,efi)

        fi = zero

        ! Comment following if statement when using filpce
        IF((etac > zero).AND.(etac .LE. one))THEN
          lambdaif=etac*(r0/r0f)*(lambdai-one)+one
          lambda0f=etac*(r0/r0f)*(lambda0-one)+one
          lambdaic=(lambdai*r0-lambdaif*r0f)/r0c
        ELSE
          lambdaif=lambdai ! False for a filament attached to a stiff crosslinker (etac = 1), only valid for etac = 0 (???)
          lambdaic=zero ! False for a stiff crosslinker (etac = 1), only valid for etac = 0 (???)
        END IF
        IF(lambdai > lambdaimax)THEN
          lambdaimax=lambdai
        END IF
        IF(lambdai .GE. 1.0d0)THEN 
          
          CALL fil(fi,ffi,dwi,ddwi,lambdaif,lambda0,lambda0f,l,r0,r0f,mu0str,beta,b0,etac)
          ! CALL filpce(lambdai, fi, dwi, ddwi)

          ! Factor of 2 accounts for the hemisphere not explicitly integrated.
          CALL sigfilfic(sfilfic,rho,lambdai,dwi,mfi,ai,ndi)

          CALL csfilfic(cfilfic,rho,lambdai,dwi,ddwi,mfi,ai,ndi)

          DO j1=1,ndi
            DO k1=1,ndi
                sfic(j1,k1)=sfic(j1,k1)+aux*sfilfic(j1,k1)
                DO l1=1,ndi
                  DO m1=1,ndi
                    cfic(j1,k1,l1,m1)=cfic(j1,k1,l1,m1)+aux*cfilfic(j1,k1,l1,m1)
                  END DO
                END DO
            END DO
          END DO

        END IF
        
        ! ================= KINETICS & ODE INTEGRATION =================                                            
        thetab_i = cb_i / cbmax
        
        ! To avoid numerical issues
        thetab_i = MIN(MAX(thetab_i, 1.0d-6), one - 1.0d-6)
        
        kon = Koff0 * Keq * exp(CHI * (1.0d0 - 2.0d0 * thetaf))                                                    
        koff_i = Koff0 * exp((fi * dx) / (kb * theta))                                                              
                                                                                                                    
        R_i = kon * cfmax * (thetaf / (1.0d0 - thetaf)) &                                                            
            - koff_i * cbmax * (thetab_i / (1.0d0 - thetab_i))                                                       
                                                                                                                    
        ! Explicit Euler Integration
        !!!!!! MAY BE REPLACED WITH A MORE STABLE INTEGRATION SCHEME !!!!!!                                                                            
        cb(node_num) = cb(node_num) + dtime * R_i                                                                    
                                                                                                                    
    ! Accumulate macroscopic pool for the NEXT time step                                                        
    ! (Note: ai is scaled by 2.0*pi because we only integrate one hemisphere)                                   
    cb_tot_new = cb_tot_new + cb(node_num) * rho * ai                                              
        ! ==============================================================   

        !v=dwi
        !rr = rr + ai * v
        !area_total = area_total + ai
        !write(*,*) etac

      end do
    end do
  end do
  end do
!
!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )
  ! IF (elem_num == 45) THEN
  !   IF(lambdaimax > 1.00d0)THEN
  !     ! write(*,*) 'WARNING (lambdamax > 1.15)!!!!!!!'
  !     write(*,*) 'lambdamax = ', lambdaimax
  !   END IF
  ! END IF

RETURN
END SUBROUTINE affclnetfic_discrete
