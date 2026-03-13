!****************************************************************************
!     Element subroutines
!****************************************************************************
      subroutine xint2D1pt(xi,w,nInttPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration
      !
      !  xi(nInttPt,2): xi,eta coordinates for the integration pts
      !  w(nInttPt):    corresponding integration weights
      !
      implicit none
      !
      integer :: nInttPt
      real(8) :: xi(1,2), w(1)


      ! Initialize
      !
      w = 0.0_8
      xi = 0.0_8


      ! Number of Gauss points
      !
      nInttPt = 1


      ! Gauss weights
      !
      w = 4.0_8
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 0.0_8
      xi(1,2) = 0.0_8


      return
      end subroutine xint2D1pt
      !
      !
!************************************************************************
      subroutine xint2D4pt(xi,w,nInttPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nInttPt,2): xi,eta coordinates for the integration pts
      !  w(nInttPt):    corresponding integration weights
      !
      implicit none
      !
      integer :: nInttPt
      real(8) :: xi(4,2), w(4)


      ! Initialize
      !
      w = 0.0_8
      xi = 0.0_8


      ! Number of Gauss points
      !
      nInttPt = 4


      ! Gauss weights
      !
      w(1) = 1.0_8
      w(2) = 1.0_8
      w(3) = 1.0_8
      w(4) = 1.0_8
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -sqrt(1.0_8/3.0_8)
      xi(1,2) = -sqrt(1.0_8/3.0_8)
      xi(2,1) = sqrt(1.0_8/3.0_8)
      xi(2,2) = -sqrt(1.0_8/3.0_8)
      xi(3,1) = -sqrt(1.0_8/3.0_8)
      xi(3,2) = sqrt(1.0_8/3.0_8)
      xi(4,1) = sqrt(1.0_8/3.0_8)
      xi(4,2) = sqrt(1.0_8/3.0_8)


      return
      end subroutine xint2D4pt
      !
      !
      subroutine xint3D1pt(xi, w, nInttPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 1 gauss point for integration
      !
      !  xi(nInttPt,3): xi, eta, zeta coordinates for the integration pts
      !  w(nInttPt):    corresponding integration weights
      
      implicit none

      integer :: nInttPt
      real(8) :: xi(1, 3), w(1)

      ! Initialize
      w = 0.0_8
      xi = 0.0_8

      ! Number of Gauss points
      nInttPt = 1

      ! Gauss weights
      w(1) = 8.0_8

      ! Gauss pt locations in master element
      xi(1, 1) = 0.0_8
      xi(1, 2) = 0.0_8
      xi(1, 3) = 0.0_8

      return
      end subroutine xint3D1pt
      !
      !
      subroutine xint3D8pt(xi, w, nInttPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 8 gauss points for integration
      !
      !  xi(nInttPt,3): xi, eta, zeta coordinates for the integration pts
      !  w(nInttPt):    corresponding integration weights
      
      implicit none

      integer :: nInttPt
      real(8) :: xi(8, 3), w(8)

      ! Initialize
      w = 0.0_8
      xi = 0.0_8

      ! Number of Gauss points
      nInttPt = 8

      ! Gauss weights
      w = 1.0_8

      ! Gauss pt locations in master element
      xi(1, 1) = -sqrt(1.0_8 / 3.0_8)
      xi(1, 2) = -sqrt(1.0_8 / 3.0_8)
      xi(1, 3) = -sqrt(1.0_8 / 3.0_8)
      xi(2, 1) = sqrt(1.0_8 / 3.0_8)
      xi(2, 2) = -sqrt(1.0_8 / 3.0_8)
      xi(2, 3) = -sqrt(1.0_8 / 3.0_8)
      xi(3, 1) = -sqrt(1.0_8 / 3.0_8)
      xi(3, 2) = sqrt(1.0_8 / 3.0_8)
      xi(3, 3) = -sqrt(1.0_8 / 3.0_8)
      xi(4, 1) = sqrt(1.0_8 / 3.0_8)
      xi(4, 2) = sqrt(1.0_8 / 3.0_8)
      xi(4, 3) = -sqrt(1.0_8 / 3.0_8)
      xi(5, 1) = -sqrt(1.0_8 / 3.0_8)
      xi(5, 2) = -sqrt(1.0_8 / 3.0_8)
      xi(5, 3) = sqrt(1.0_8 / 3.0_8)
      xi(6, 1) = sqrt(1.0_8 / 3.0_8)
      xi(6, 2) = -sqrt(1.0_8 / 3.0_8)
      xi(6, 3) = sqrt(1.0_8 / 3.0_8)
      xi(7, 1) = -sqrt(1.0_8 / 3.0_8)
      xi(7, 2) = sqrt(1.0_8 / 3.0_8)
      xi(7, 3) = sqrt(1.0_8 / 3.0_8)
      xi(8, 1) = sqrt(1.0_8 / 3.0_8)
      xi(8, 2) = sqrt(1.0_8 / 3.0_8)
      xi(8, 3) = sqrt(1.0_8 / 3.0_8)

      return
      end subroutine xint3D8pt
      !
      !
      subroutine xintSurf2D1pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for surface integration
      !
      !  xLocal(nInttPt): x coordinates for the integration pts
      !  yLocal(nInttPt): y coordinates for the integration pts
      !  w(nInttPt):    corresponding integration weights

      implicit none

      integer :: face
      real(8) :: xLocal(1), yLocal(1), w(1)
      real(8), parameter :: zero = 0.0_8, one = 1.0_8, two = 2.0_8

      ! Gauss weights
      w(1) = two

      ! Gauss pt locations in master element
      select case(face)
      case(1)
         xLocal(1) = zero
         yLocal(1) = -one
      case(2)
         xLocal(1) = one
         yLocal(1) = zero
      case(3)
         xLocal(1) = zero
         yLocal(1) = one
      case(4)
         xLocal(1) = -one
         yLocal(1) = zero
      case default
         write(*,*) 'face.ne.1,2,3,4'
         call exit
      end select

      end subroutine xintSurf2D1pt
      !
      !
      subroutine xintSurf2D2pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nInttPt): x coordinates for the integration pts
      !  yLocal(nInttPt): y coordinates for the integration pts
      !  w(nInttPt):    corresponding integration weights

      implicit none

      integer :: face
      real(8) :: xLocal(2), yLocal(2), w(2)
      real(8), parameter :: one = 1.0_8, three = 3.0_8

      ! Gauss weights
      w(1) = one
      w(2) = one

      ! Gauss pt locations in master element
      select case(face)
      case(1)
         xLocal(1) = -sqrt(one/three)
         yLocal(1) = -one
         xLocal(2) = sqrt(one/three)
         yLocal(2) = -one
      case(2)
         xLocal(1) = one
         yLocal(1) = -sqrt(one/three)
         xLocal(2) = one
         yLocal(2) = sqrt(one/three)
      case(3)
         xLocal(1) = -sqrt(one/three)
         yLocal(1) = one
         xLocal(2) = sqrt(one/three)
         yLocal(2) = one
      case(4)
         xLocal(1) = -one
         yLocal(1) = sqrt(one/three)
         xLocal(2) = -one
         yLocal(2) = -sqrt(one/three)
      case default
         write(*,*) 'face.ne.1,2,3,4'
         call exit
      end select

      end subroutine xintSurf2D2pt
      !
      !
      subroutine xintSurf2D3pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 3 gauss points for surface integration
      !
      !  xLocal(nInttPt): x coordinates for the integration pts
      !  yLocal(nInttPt): y coordinates for the integration pts
      !  w(nInttPt):    corresponding integration weights

      implicit none

      integer :: face
      real(8) :: xLocal(3), yLocal(3), w(3)
      real(8), parameter :: zero = 0.0_8, one = 1.0_8, two = 2.0_8, &
                   three = 3.0_8, five = 5.0_8, eight = 8.0_8, &
                   nine = 9.0_8

      ! Gauss weights
      w(1) = five/nine
      w(2) = eight/nine
      w(3) = five/nine

      ! Gauss pt locations in master element
      select case(face)
      case(1)
         xLocal(1) = -sqrt(three/five)
         yLocal(1) = -one
         xLocal(2) = zero
         yLocal(2) = -one
         xLocal(3) = sqrt(three/five)
         yLocal(3) = -one
      case(2)
         xLocal(1) = one
         yLocal(1) = -sqrt(three/five)
         xLocal(2) = one
         yLocal(2) = zero
         xLocal(3) = one
         yLocal(3) = sqrt(three/five)
      case(3)
         xLocal(1) = -sqrt(three/five)
         yLocal(1) = one
         xLocal(2) = zero
         yLocal(2) = one
         xLocal(3) = sqrt(three/five)
         yLocal(3) = one
      case(4)
         xLocal(1) = -one
         yLocal(1) = sqrt(three/five)
         xLocal(2) = -one
         yLocal(2) = zero
         xLocal(3) = -one
         yLocal(3) = -sqrt(three/five)
      case default
         write(*,*) 'face.ne.1,2,3,4'
         call exit
      end select

      end subroutine xintSurf2D3pt
      !
      !
      subroutine xintSurf3D1pt(face, xLocal, yLocal, zLocal, w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 1 gauss point for surface integration
      !
      !  xLocal(nInttPt): x coordinates for the integration pts
      !  yLocal(nInttPt): y coordinates for the integration pts
      !  zLocal(nInttPt): z coordinates for the integration pts
      !  w(nInttPt):    corresponding integration weights

      implicit none

      integer :: face
      real(8) :: xLocal(1), yLocal(1), zLocal(1), w(1)
      real(8), parameter :: zero = 0.0_8, one = 1.0_8, four = 4.0_8

      ! Gauss weights
      w(1) = four

      ! Gauss pt locations in master element
      select case(face)
      case(1)
         xLocal(1) = zero
         yLocal(1) = zero
         zLocal(1) = -one
      case(2)
         xLocal(1) = zero
         yLocal(1) = zero
         zLocal(1) = one
      case(3)
         xLocal(1) = zero
         yLocal(1) = -one
         zLocal(1) = zero
      case(4)
         xLocal(1) = one
         yLocal(1) = zero
         zLocal(1) = zero
      case(5)
         xLocal(1) = zero
         yLocal(1) = one
         zLocal(1) = zero
      case(6)
         xLocal(1) = -one
         yLocal(1) = zero
         zLocal(1) = zero
      case default
         write(*,*) 'face.ne.1,2,3,4,5,6'
         call exit
      end select

      end subroutine xintSurf3D1pt
      !
      !
      subroutine xintSurf3D4pt(face, xLocal, yLocal, zLocal, w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 4 gauss points for surface integration
      !
      !  xLocal(nInttPt): x coordinates for the integration pts
      !  yLocal(nInttPt): y coordinates for the integration pts
      !  zLocal(nInttPt): z coordinates for the integration pts
      !  w(nInttPt):    corresponding integration weights

      implicit none

      integer :: face
      real(8) :: xLocal(4), yLocal(4), zLocal(4), w(4)
      real(8), parameter :: one = 1.0_8, three = 3.0_8

      ! Gauss weights
      w = one

      ! Gauss pt locations in master element
      select case(face)
      case(1)
         xLocal = [-sqrt(one/three), sqrt(one/three), sqrt(one/three), -sqrt(one/three)]
         yLocal = [-sqrt(one/three), -sqrt(one/three), sqrt(one/three), sqrt(one/three)]
         zLocal = [-one, -one, -one, -one]
      case(2)
         xLocal = [-sqrt(one/three), sqrt(one/three), sqrt(one/three), -sqrt(one/three)]
         yLocal = [-sqrt(one/three), -sqrt(one/three), sqrt(one/three), sqrt(one/three)]
         zLocal = [one, one, one, one]
      case(3)
         xLocal = [-sqrt(one/three), sqrt(one/three), sqrt(one/three), -sqrt(one/three)]
         yLocal = [-one, -one, -one, -one]
         zLocal = [-sqrt(one/three), -sqrt(one/three), sqrt(one/three), sqrt(one/three)]
      case(4)
         xLocal = [one, one, one, one]
         yLocal = [-sqrt(one/three), sqrt(one/three), sqrt(one/three), -sqrt(one/three)]
         zLocal = [-sqrt(one/three), -sqrt(one/three), sqrt(one/three), sqrt(one/three)]
      case(5)
         xLocal = [-sqrt(one/three), sqrt(one/three), sqrt(one/three), -sqrt(one/three)]
         yLocal = [one, one, one, one]
         zLocal = [-sqrt(one/three), -sqrt(one/three), sqrt(one/three), sqrt(one/three)]
      case(6)
         xLocal = [-one, -one, -one, -one]
         yLocal = [-sqrt(one/three), sqrt(one/three), sqrt(one/three), -sqrt(one/three)]
         zLocal = [-sqrt(one/three), -sqrt(one/three), sqrt(one/three), sqrt(one/three)]
      case default
         write(*,*) 'face.ne.1,2,3,4,5,6'
         call exit
      end select

      end subroutine xintSurf3D4pt
      !
      !
      subroutine calcShape2DLinear(nInttPt, xi_int, intpt, sh, dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      !                          eta
      !   4-----------3          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          O--------- xi
      !   1-----------2        origin at center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer :: intpt, nInttPt
      real(8) :: xi_int(nInttPt, 2), sh(4), dshxi(4, 2), xi, eta
      real(8), parameter :: zero = 0.0_8, one = 1.0_8, fourth = 1.0_8 / 4.0_8
      

      ! Location in the master element
      !
      xi = xi_int(intpt, 1)
      eta = xi_int(intpt, 2)
      
      
      ! The shape functions
      !
      sh(1) = fourth * (one - xi) * (one - eta)
      sh(2) = fourth * (one + xi) * (one - eta)
      sh(3) = fourth * (one + xi) * (one + eta)
      sh(4) = fourth * (one - xi) * (one + eta)
      
      
      ! The first derivatives
      !
      dshxi(1, 1) = -fourth * (one - eta)
      dshxi(1, 2) = -fourth * (one - xi)
      dshxi(2, 1) = fourth * (one - eta)
      dshxi(2, 2) = -fourth * (one + xi)
      dshxi(3, 1) = fourth * (one + eta)
      dshxi(3, 2) = fourth * (one + xi)
      dshxi(4, 1) = -fourth * (one + eta)
      dshxi(4, 2) = fourth * (one - xi)
      

      return
      end subroutine calcShape2DLinear
      !
      !
      subroutine calcShape3DLinear(nInttPt, xi_int, intpt, sh, dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses an 8-node linear 3D element as shown
      !
      !      8-----------7
      !     /|          /|       zeta
      !    / |         / |       
      !   5-----------6  |       |     eta
      !   |  |        |  |       |   /
      !   |  |        |  |       |  /
      !   |  4--------|--3       | /
      !   | /         | /        |/
      !   |/          |/         O--------- xi
      !   1-----------2        origin at cube center
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer :: intpt, nInttPt
      real(8) :: xi_int(nInttPt, 3), sh(8), dshxi(8, 3)
      real(8) :: d2shxi(8,3,3)
      real(8) :: xi, eta, zeta
      real(8), parameter :: zero = 0.0_8, one = 1.0_8, eighth = 1.0_8 / 8.0_8

      ! Location in the master element
      xi = xi_int(intpt, 1)
      eta = xi_int(intpt, 2)
      zeta = xi_int(intpt, 3)

      ! The shape functions
      sh(1) = eighth * (one - xi) * (one - eta) * (one - zeta)
      sh(2) = eighth * (one + xi) * (one - eta) * (one - zeta)
      sh(3) = eighth * (one + xi) * (one + eta) * (one - zeta)
      sh(4) = eighth * (one - xi) * (one + eta) * (one - zeta)
      sh(5) = eighth * (one - xi) * (one - eta) * (one + zeta)
      sh(6) = eighth * (one + xi) * (one - eta) * (one + zeta)
      sh(7) = eighth * (one + xi) * (one + eta) * (one + zeta)
      sh(8) = eighth * (one - xi) * (one + eta) * (one + zeta)

      ! The first derivatives
      dshxi(1, 1) = -eighth * (one - eta) * (one - zeta)
      dshxi(1, 2) = -eighth * (one - xi) * (one - zeta)
      dshxi(1, 3) = -eighth * (one - xi) * (one - eta)
      dshxi(2, 1) = eighth * (one - eta) * (one - zeta)
      dshxi(2, 2) = -eighth * (one + xi) * (one - zeta)
      dshxi(2, 3) = -eighth * (one + xi) * (one - eta)
      dshxi(3, 1) = eighth * (one + eta) * (one - zeta)
      dshxi(3, 2) = eighth * (one + xi) * (one - zeta)
      dshxi(3, 3) = -eighth * (one + xi) * (one + eta)
      dshxi(4, 1) = -eighth * (one + eta) * (one - zeta)
      dshxi(4, 2) = eighth * (one - xi) * (one - zeta)
      dshxi(4, 3) = -eighth * (one - xi) * (one + eta)
      dshxi(5, 1) = -eighth * (one - eta) * (one + zeta)
      dshxi(5, 2) = -eighth * (one - xi) * (one + zeta)
      dshxi(5, 3) = eighth * (one - xi) * (one - eta)
      dshxi(6, 1) = eighth * (one - eta) * (one + zeta)
      dshxi(6, 2) = -eighth * (one + xi) * (one + zeta)
      dshxi(6, 3) = eighth * (one + xi) * (one - eta)
      dshxi(7, 1) = eighth * (one + eta) * (one + zeta)
      dshxi(7, 2) = eighth * (one + xi) * (one + zeta)
      dshxi(7, 3) = eighth * (one + xi) * (one + eta)
      dshxi(8, 1) = -eighth * (one + eta) * (one + zeta)
      dshxi(8, 2) = eighth * (one - xi) * (one + zeta)
      dshxi(8, 3) = eighth * (one - xi) * (one + eta)

      !
      ! The second derivatives
      !
      d2shxi = zero
      d2shxi(1,1,2) = eighth*(one - zeta)
      d2shxi(1,2,1) = d2shxi(1,1,2)
      d2shxi(1,1,3) = eighth*(one - eta)
      d2shxi(1,3,1) = d2shxi(1,1,3)
      d2shxi(1,2,3) = eighth*(one - xi)
      d2shxi(1,3,2) = d2shxi(1,2,3)
      d2shxi(2,1,2) = -eighth*(one - zeta)
      d2shxi(2,2,1) = d2shxi(2,1,2)
      d2shxi(2,1,3) = -eighth*(one - eta)
      d2shxi(2,3,1) = d2shxi(2,1,3)
      d2shxi(2,2,3) = eighth*(one + xi)
      d2shxi(2,3,2) = d2shxi(2,2,3)
      d2shxi(3,1,2) = eighth*(one - zeta)
      d2shxi(3,2,1) = d2shxi(2,1,2)
      d2shxi(3,1,3) = -eighth*(one + eta)
      d2shxi(3,3,1) = d2shxi(2,1,3)
      d2shxi(3,2,3) = -eighth*(one + xi)
      d2shxi(3,3,2) = d2shxi(2,2,3)
      d2shxi(4,1,2) = -eighth*(one - zeta)
      d2shxi(4,2,1) = d2shxi(2,1,2)
      d2shxi(4,1,3) = eighth*(one + eta)
      d2shxi(4,3,1) = d2shxi(2,1,3)
      d2shxi(4,2,3) = -eighth*(one - xi)
      d2shxi(4,3,2) = d2shxi(2,2,3)
      d2shxi(5,1,2) = eighth*(one + zeta)
      d2shxi(5,2,1) = d2shxi(2,1,2)
      d2shxi(5,1,3) = -eighth*(one - eta)
      d2shxi(5,3,1) = d2shxi(2,1,3)
      d2shxi(5,2,3) = -eighth*(one - xi)
      d2shxi(5,3,2) = d2shxi(2,2,3)
      d2shxi(6,1,2) = eighth*(one + zeta)
      d2shxi(6,2,1) = d2shxi(2,1,2)
      d2shxi(6,1,3) = eighth*(one - eta)
      d2shxi(6,3,1) = d2shxi(2,1,3)
      d2shxi(6,2,3) = -eighth*(one + xi)
      d2shxi(6,3,2) = d2shxi(2,2,3)
      d2shxi(7,1,2) = eighth*(one + zeta)
      d2shxi(7,2,1) = d2shxi(2,1,2)
      d2shxi(7,1,3) = eighth*(one + eta)
      d2shxi(7,3,1) = d2shxi(2,1,3)
      d2shxi(7,2,3) = eighth*(one + xi)
      d2shxi(7,3,2) = d2shxi(2,2,3)
      d2shxi(8,1,2) = -eighth*(one + zeta)
      d2shxi(8,2,1) = d2shxi(2,1,2)
      d2shxi(8,1,3) = -eighth*(one + eta)
      d2shxi(8,3,1) = d2shxi(2,1,3)
      d2shxi(8,2,3) = eighth*(one - xi)
      d2shxi(8,3,2) = d2shxi(2,2,3)

      return
      end subroutine calcShape3DLinear
      !
      !
      subroutine computeSurf(xLocal, yLocal, face, coords, sh, ds)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the length ds, so that one can
      !  do the numerical integration on the boundary for fluxes 
      !  on the 4-node quadrilateral elements

      implicit none

      integer :: face
      real(8) :: xLocal, yLocal, ds, dshxi(4, 2), sh(4)
      real(8) :: dXdXi, dXdEta, dYdXi, dYdEta
      real(8) :: coords(2, 4)
      real(8), parameter :: one = 1.0_8, fourth = 1.0_8 / 4.0_8
      real(8) :: normal(2, 1)

      sh(1) = fourth * (one - xLocal) * (one - yLocal)
      sh(2) = fourth * (one + xLocal) * (one - yLocal)
      sh(3) = fourth * (one + xLocal) * (one + yLocal)
      sh(4) = fourth * (one - xLocal) * (one + yLocal)
      
      dshxi(1, 1) = -fourth * (one - yLocal)
      dshxi(1, 2) = -fourth * (one - xLocal)
      dshxi(2, 1) = fourth * (one - yLocal)
      dshxi(2, 2) = -fourth * (one + xLocal)
      dshxi(3, 1) = fourth * (one + yLocal)
      dshxi(3, 2) = fourth * (one + xLocal)
      dshxi(4, 1) = -fourth * (one + yLocal)
      dshxi(4, 2) = fourth * (one - xLocal)

      dXdXi = sum(dshxi(:, 1) * coords(1, :))
      dXdEta = sum(dshxi(:, 2) * coords(1, :))
      dYdXi = sum(dshxi(:, 1) * coords(2, :))
      dYdEta = sum(dshxi(:, 2) * coords(2, :))

      ! Jacobian of the mapping
      if (face == 2 .or. face == 4) then
         ds = sqrt(dXdEta**2 + dYdEta**2)
      elseif (face == 1 .or. face == 3) then
         ds = sqrt(dXdXi**2 + dYdXi**2)
      else
         write(*, *) 'never should get here'
         call exit
      endif

      ! Surface normal, outward pointing in this case. Useful for
      !  ``follower'' type loads. The normal is referential or spatial
      !  depending on which coords were supplied to this subroutine
      !  (NOT fully tested)
      if (face == 2 .or. face == 4) then
         normal(1, 1) = dYdEta / sqrt(dXdEta**2 + dYdEta**2)
         normal(2, 1) = -dXdEta / sqrt(dXdEta**2 + dYdEta**2)
         if (face == 4) normal = -normal
      elseif (face == 1 .or. face == 3) then
         normal(1, 1) = dYdXi / sqrt(dXdXi**2 + dYdXi**2)
         normal(2, 1) = -dXdXi / sqrt(dXdXi**2 + dYdXi**2)
         if (face == 3) normal = -normal
      else
         write(*, *) 'never should get here'
         call exit
      endif

      return
      end subroutine computeSurf
      !
      !
      subroutine computeSurf3D(xLocal, yLocal, zLocal, face, coords, sh, dA)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the area dA, so that one can
      !  do the numerical integration on the boundary for fluxes 
      !  on the 8-node brick elements

      implicit none

      integer :: face, stat, i, j, k
      real(8) :: xLocal, yLocal, zLocal, dA, dshxi(8, 3), sh(8)
      real(8) :: coords(3, 8), mapJ(3, 3), mag, normal(3)
      real(8) :: dXdXi, dXdEta, dXdZeta, dYdXi, dYdEta, dYdZeta
      real(8) :: dZdXi, dZdEta, dZdZeta
      real(8), parameter :: one = 1.0_8, two = 2.0_8, eighth = 1.0_8 / 8.0_8, zero = 0.0_8

      ! The shape functions
      sh(1) = eighth * (one - xLocal) * (one - yLocal) * (one - zLocal)
      sh(2) = eighth * (one + xLocal) * (one - yLocal) * (one - zLocal)
      sh(3) = eighth * (one + xLocal) * (one + yLocal) * (one - zLocal)
      sh(4) = eighth * (one - xLocal) * (one + yLocal) * (one - zLocal)
      sh(5) = eighth * (one - xLocal) * (one - yLocal) * (one + zLocal)
      sh(6) = eighth * (one + xLocal) * (one - yLocal) * (one + zLocal)
      sh(7) = eighth * (one + xLocal) * (one + yLocal) * (one + zLocal)
      sh(8) = eighth * (one - xLocal) * (one + yLocal) * (one + zLocal)

      ! Shape function derivatives
      dshxi(1, 1) = -eighth * (one - yLocal) * (one - zLocal)
      dshxi(1, 2) = -eighth * (one - xLocal) * (one - zLocal)
      dshxi(1, 3) = -eighth * (one - xLocal) * (one - yLocal)
      dshxi(2, 1) = eighth * (one - yLocal) * (one - zLocal)
      dshxi(2, 2) = -eighth * (one + xLocal) * (one - zLocal)
      dshxi(2, 3) = -eighth * (one + xLocal) * (one - yLocal)
      dshxi(3, 1) = eighth * (one + yLocal) * (one - zLocal)
      dshxi(3, 2) = eighth * (one + xLocal) * (one - zLocal)
      dshxi(3, 3) = -eighth * (one + xLocal) * (one + yLocal)
      dshxi(4, 1) = -eighth * (one + yLocal) * (one - zLocal)
      dshxi(4, 2) = eighth * (one - xLocal) * (one - zLocal)
      dshxi(4, 3) = -eighth * (one - xLocal) * (one + yLocal)
      dshxi(5, 1) = -eighth * (one - yLocal) * (one + zLocal)
      dshxi(5, 2) = -eighth * (one - xLocal) * (one + zLocal)
      dshxi(5, 3) = eighth * (one - xLocal) * (one - yLocal)
      dshxi(6, 1) = eighth * (one - yLocal) * (one + zLocal)
      dshxi(6, 2) = -eighth * (one + xLocal) * (one + zLocal)
      dshxi(6, 3) = eighth * (one + xLocal) * (one - yLocal)
      dshxi(7, 1) = eighth * (one + yLocal) * (one + zLocal)
      dshxi(7, 2) = eighth * (one + xLocal) * (one + zLocal)
      dshxi(7, 3) = eighth * (one + xLocal) * (one + yLocal)
      dshxi(8, 1) = -eighth * (one + yLocal) * (one + zLocal)
      dshxi(8, 2) = eighth * (one - xLocal) * (one + zLocal)
      dshxi(8, 3) = eighth * (one - xLocal) * (one + yLocal)

      dXdXi = zero
      dXdEta = zero
      dXdZeta = zero
      dYdXi = zero
      dYdEta = zero
      dYdZeta = zero
      dZdXi = zero
      dZdEta = zero
      dZdZeta = zero
      do k = 1, 8
         dXdXi = dXdXi + dshxi(k, 1) * coords(1, k)
         dXdEta = dXdEta + dshxi(k, 2) * coords(1, k)
         dXdZeta = dXdZeta + dshxi(k, 3) * coords(1, k)
         dYdXi = dYdXi + dshxi(k, 1) * coords(2, k)
         dYdEta = dYdEta + dshxi(k, 2) * coords(2, k)
         dYdZeta = dYdZeta + dshxi(k, 3) * coords(2, k)
         dZdXi = dZdXi + dshxi(k, 1) * coords(3, k)
         dZdEta = dZdEta + dshxi(k, 2) * coords(3, k)
         dZdZeta = dZdZeta + dshxi(k, 3) * coords(3, k)
      end do

      ! Jacobian of the mapping
      select case (face)
      case (1, 2)
         ! zeta = constant on this face
         dA = sqrt((dYdXi * dZdEta - dYdEta * dZdXi)**two + &
             (dXdXi * dZdEta - dXdEta * dZdXi)**two + &
             (dXdXi * dYdEta - dXdEta * dYdXi)**two)
      case (3, 5)
         ! eta = constant on this face
         dA = sqrt((dYdXi * dZdZeta - dYdZeta * dZdXi)**two + &
             (dXdXi * dZdZeta - dXdZeta * dZdXi)**two + &
             (dXdXi * dYdZeta - dXdZeta * dYdXi)**two)
      case (4, 6)
         ! xi = constant on this face
         dA = sqrt((dYdEta * dZdZeta - dYdZeta * dZdEta)**two + &
             (dXdEta * dZdZeta - dXdZeta * dZdEta)**two + &
             (dXdEta * dYdZeta - dXdZeta * dYdEta)**two)
      case default
         write(*, *) 'never should get here'
         call exit
      end select

      ! Surface normal, outward pointing in this case. Useful for
      !  ``follower'' type loads. The normal is referential or spatial
      !  depending on which coords were supplied to this subroutine
      !  (NOT fully tested)
      select case (face)
      case (1, 2)
         ! zeta = constant on this face
         normal(1) = dYdXi * dZdEta - dYdEta * dZdXi
         normal(2) = dXdXi * dZdEta - dXdEta * dZdXi
         normal(3) = dXdXi * dYdEta - dXdEta * dYdXi
         if (face == 1) normal = -normal
      case (3, 5)
         ! eta = constant on this face
         normal(1) = dYdXi * dZdZeta - dYdZeta * dZdXi
         normal(2) = dXdXi * dZdZeta - dXdZeta * dZdXi
         normal(3) = dXdXi * dYdZeta - dXdZeta * dYdXi
         if (face == 5) normal = -normal
      case (4, 6)
         ! xi = constant on this face
         normal(1) = dYdEta * dZdZeta - dYdZeta * dZdEta
         normal(2) = dXdEta * dZdZeta - dXdZeta * dZdEta
         normal(3) = dXdEta * dYdZeta - dXdZeta * dYdEta
         if (face == 6) normal = -normal
      case default
         write(*, *) 'never should get here'
         call exit
      end select
      mag = sqrt(normal(1)**two + normal(2)**two + normal(3)**two)
      normal = normal / mag

      end subroutine computeSurf3D
      !
      !
      subroutine mapShape2D(nNode, dshxi, coords, dsh, detMapJ, stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer :: i, j, k, nNode, stat
      !
      real(8) :: dshxi(nNode, 2), dsh(nNode, 2), coords(3, nNode), mapJ(2, 2), &
              mapJ_inv(2, 2), detMapJ
      !
      real(8), parameter :: zero = 0.0_8, one = 1.0_8, two = 2.0_8, &
                   half = 0.5_8, fourth = 0.25_8, eighth = 1.0_8 / 8.0_8

      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i = 1, 2
         do j = 1, 2
         do k = 1, nNode
            mapJ(i, j) = mapJ(i, j) + dshxi(k, i) * coords(j, k)
         end do
         end do
      end do

      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ, mapJ_inv, detMapJ, stat)
      if (stat == 0) then
         write(*, *) 'Problem: detF.lt.zero in mapShape2D'
         call exit
      end if

      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv, transpose(dshxi)))

      return
      end subroutine mapShape2D
      !
      !
      subroutine mapShape2Da(nNode, dshxi, coords, dsh, detMapJ, stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      ! This subroutine is exactly the same as the regular mapShape2D
      !  with the exception that coords(2,nNode) here and coords(3,nNode)
      !  in the regular.  I have noticed that a "heat transfer" and 
      !  "static" step uses MCRD=2, but for "coupled-temperature-displacement"
      !  you will get MCRD=3, even for a plane analysis.
      !
      implicit none
      !
      integer :: i, j, k, nNode, stat
      !
      real(8) :: dshxi(nNode, 2), dsh(nNode, 2), coords(2, nNode), mapJ(2, 2), &
              mapJ_inv(2, 2), detMapJ
      !
      real(8), parameter :: zero = 0.0_8, one = 1.0_8, two = 2.0_8, &
                   half = 0.5_8, fourth = 0.25_8, eighth = 1.0_8 / 8.0_8

      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i = 1, 2
         do j = 1, 2
         do k = 1, nNode
            mapJ(i, j) = mapJ(i, j) + dshxi(k, i) * coords(j, k)
         end do
         end do
      end do

      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ, mapJ_inv, detMapJ, stat)
      if (stat == 0) then
         write(*, *) 'Problem: detF.lt.zero in mapShape2Da'
         call exit
      end if

      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv, transpose(dshxi)))

      return
      end subroutine mapShape2Da
      !
      !
      subroutine mapShape3D(nNode, dshxi, coords, dsh, detMapJ, stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain. This subroutine works for both 8-node
      !  linear and 20-node quadratic 3D elements.
      !
      implicit none

      integer :: i, j, k, nNode, stat
      real(8) :: dshxi(nNode, 3), dsh(nNode, 3), coords(3, nNode)
      real(8) :: mapJ(3, 3), mapJ_inv(3, 3), detMapJ
      real(8), parameter :: zero = 0.0_8, one = 1.0_8, two = 2.0_8, &
                   half = 0.5_8, fourth = 0.25_8, eighth = 1.0_8 / 8.0_8

      ! Calculate the mapping Jacobian matrix:
      mapJ = zero
      do i = 1, 3
         do j = 1, 3
         do k = 1, nNode
            mapJ(i, j) = mapJ(i, j) + dshxi(k, i) * coords(j, k)
         end do
         end do
      end do

      ! Calculate the inverse and the determinant of Jacobian
      call matInv3Dd(mapJ, mapJ_inv, detMapJ, stat)
      if (stat == 0) then
         write(*, *) 'Problem: detF.lt.zero in mapShape3D'
         call exit
      end if

      ! Calculate first derivatives wrt x, y, z
      dsh = transpose(matmul(mapJ_inv, transpose(dshxi)))

      ! The second derivatives may be calculated.
      !

      return
      end subroutine mapShape3D
      
