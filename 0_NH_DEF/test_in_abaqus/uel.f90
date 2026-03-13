
      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and 
      !   the dummy mesh needs to have the same number of integ
      !   points.  You must set that parameter value here.
      !
      !  ElemOffset
      !   Offset between element numbers on the real mesh and
      !    dummy mesh.  That is set in the input file, and 
      !    that value must be set here the same.

      integer numElem,ElemOffset,err
      INTEGER NWP,NELEM,NCH,NSDV
      DOUBLE PRECISION  ONE, TWO, THREE, FOUR, SIX, ZERO
      DOUBLE PRECISION HALF,THIRD
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of UEL elements used here
      parameter(numElem=1)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ElemOffset=1000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PARAMETER(NELEM=1, NSDV=1)
      PARAMETER(ZERO=0.D0, ONE=1.0D0,TWO=2.0D0)
      PARAMETER(THREE=3.0D0,FOUR=4.0D0,SIX=6.0D0)
      PARAMETER(HALF=0.5d0,THIRD=1.d0/3.d0)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable :: globalSdv(:,:,:)

      end module global
subroutine AssembleElement(nDim, nNode, ndofel, Ru, Kuu, rhs, amatrx)

   ! Subroutine to assemble the local elements residual and tangent

   implicit none

   integer :: i, j, nDim, nNode, ndofel, nDofN
   integer :: A11, A12, B11, B12

   real(8), intent(in) :: Ru(nDim*nNode, 1), Kuu(nDim*nNode, nDim*nNode)
   real(8), intent(out) :: rhs(ndofel, 1), amatrx(ndofel, ndofel)

   ! Total number of degrees of freedom per node
   nDofN = ndofel / nNode

   ! Initialize
   rhs = 0.d0
   amatrx = 0.d0

   if (nDim == 2) then
      !
      ! Assemble the element level residual
      !
      do i = 1, nNode
         A11 = nDofN * (i - 1) + 1
         A12 = nDim * (i - 1) + 1
         !
         ! Displacement
         !
         rhs(A11, 1) = Ru(A12, 1)
         rhs(A11 + 1, 1) = Ru(A12 + 1, 1)
      end do
      !
      ! Assemble the element level tangent matrix
      !
      do i = 1, nNode
         do j = 1, nNode
            A11 = nDofN * (i - 1) + 1
            A12 = nDim * (i - 1) + 1
            B11 = nDofN * (j - 1) + 1
            B12 = nDim * (j - 1) + 1
            ! Displacement
            amatrx(A11, B11) = Kuu(A12, B12)
            amatrx(A11, B11 + 1) = Kuu(A12, B12 + 1)
            amatrx(A11 + 1, B11) = Kuu(A12 + 1, B12)
            amatrx(A11 + 1, B11 + 1) = Kuu(A12 + 1, B12 + 1)
         end do
      end do

   elseif (nDim == 3) then
      !
      ! Assemble the element level residual
      !
      do i = 1, nNode
         A11 = nDofN * (i - 1) + 1
         A12 = nDim * (i - 1) + 1
         !
         ! Displacement
         !
         rhs(A11, 1) = Ru(A12, 1)
         rhs(A11 + 1, 1) = Ru(A12 + 1, 1)
         rhs(A11 + 2, 1) = Ru(A12 + 2, 1)
      end do
      !
      ! Assemble the element level tangent matrix
      !
      do i = 1, nNode
         do j = 1, nNode
            A11 = nDofN * (i - 1) + 1
            A12 = nDim * (i - 1) + 1
            B11 = nDofN * (j - 1) + 1
            B12 = nDim * (j - 1) + 1
            !
            ! Displacement
            !
            amatrx(A11, B11) = Kuu(A12, B12)
            amatrx(A11, B11 + 1) = Kuu(A12, B12 + 1)
            amatrx(A11, B11 + 2) = Kuu(A12, B12 + 2)
            amatrx(A11 + 1, B11) = Kuu(A12 + 1, B12)
            amatrx(A11 + 1, B11 + 1) = Kuu(A12 + 1, B12 + 1)
            amatrx(A11 + 1, B11 + 2) = Kuu(A12 + 1, B12 + 2)
            amatrx(A11 + 2, B11) = Kuu(A12 + 2, B12)
            amatrx(A11 + 2, B11 + 1) = Kuu(A12 + 2, B12 + 1)
            amatrx(A11 + 2, B11 + 2) = Kuu(A12 + 2, B12 + 2)
         end do
      end do

   else
      write(*, *) 'How did you get nDim=', nDim
      call exit
   end if

   return
end subroutine AssembleElement
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

      return
      end subroutine mapShape3D
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
              write(*, *) 'stat, error:', stat, err
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
         ! XXXXXXXXXXXXXXXXXXXXXX

         ! Initialize the residual and tangent matrices to zero.
         Ru = 0.0d0
         Kuu = 0.0d0

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
         end do

         ! Obtain current nodal coordinates
         do i = 1, NNODE
            do j = 1, nDim
               coordsC(j, i) = coords(j, i) + u(i, j)
            end do
         end do

         ! Impose any time-stepping changes on the increments of chemical potential or displacement if you want
         ! displacement increment, based on element diagonal
         Le = sqrt((coordsC(1, 1) - coordsC(1, 7))**2 + &
                 (coordsC(2, 1) - coordsC(2, 7))**2 + &
                 (coordsC(3, 1) - coordsC(3, 7))**2)
         ! add some kinf of flag here???
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
            !  centriod, that means xi=eta=zeta=0.0, and nInttPt=1
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
            do intpt = 1, nInttPt

               ! Obtain state variables from previous increment
               !
               if ((KINC <= 1) .and. (KSTEP == 1)) then
               ! this is the first increment, of the first step
               !  give initial conditions (or just anything)
               statev = 0.9999d0 !initial determinant of the deformation gradient
               prev_statev = 0.9999d0
               else
               ! this is not the first increment, read old values
               statev = SVARS(1 + jj : nsdv + jj)
               prev_statev = SVARS(1 + jj : nsdv + jj)
               endif
               write(*, *) 'SVARS=', SVARS

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
               F_t = ((detFc_tau / detF_tau)**(1.0d0 / 3.0d0)) * F_t
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
                       JELEM, intpt, KSTEP, KINC)
               !
               !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               write(*, *) 'detF_t=', prev_statev(1)
               write(*, *) 'detF_tau=', statev(1)
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
         write(*, *) 'umeror=', umeror
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
            (-matmul(transpose(Bmat), Smat) + BodyForceRes)      ! Compute/update the displacement tangent matrix
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
            Kuu = Kuu + detMapJC * w(intpt) * &
                 (matmul(matmul(transpose(Gmat), Amat), Gmat) + &
                 matmul(transpose(Gmat), matmul(Qmat, (G0mat - Gmat))))
          else
            ! This is the tangent NOT using the F-bar method with all other elements
            Kuu = Kuu + detMapJC * w(intpt) * &
                 (matmul(matmul(transpose(Gmat), Amat), Gmat))
          end if
         
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
         Ru, Kuu, &
         rhs, amatrx)
   !      write(*,*) rhs(:,1)
   !      write(*,*) amatrx     
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------
      return 
      end subroutine U3D8!************************************************************************
!
! User element for transient fluid permeation, and large 
!  elastic deformation in 2D or 3D.  This is for plane strain,
!  axisymetric, and 3D.
!
! Solution variables (or nodal variables) are the displacements and the
!  chemical potential.
! 
! This subroutine is for the following element types
!  > two-dimensional 4 node isoparametric element as shown below
!       with 1pt (reduced) or 4pt (full) gauss integration.
!  > three-dimensional 8 node isoparametric element as shown below
!       with 1pt (reduced) or 8pt (full) gauss integration.
!
! In order to avoid locking for the fully-integrated element, we
!  use the F-bar method of de Souza Neto (1996).
!
!  Mechanical, traction- and pressure-type boundary conditions 
!   may be applied to the dummy mesh using the Abaqus built-in 
!   commands *Dload or *Dsload.
!
! Surface flux boundary conditions are supported in the following
!  elements.  Based on our convention, the face on which the fliud
!  flux is applied is the "label", i.e.
!  - U1,U2,U3,U4,... refer to fluid fluxes applied to faces 
!                     1,2,3,4,... respectively,
!
!     
!              A eta (=xi_2)
!  4-node      |
!   quad       |Face 3
!        4-----------3
!        |     |     |
!        |     |     |
!  Face 4|     ------|---> xi (=xi_1)
!        |           | Face2
!        |           |
!        1-----------2
!          Face 1
!
!
!  8-node     8-----------7
!  brick     /|          /|       zeta
!           / |         / |       
!          5-----------6  |       |     eta
!          |  |        |  |       |   /
!          |  |        |  |       |  /
!          |  4--------|--3       | /
!          | /         | /        |/
!          |/          |/         O--------- xi
!          1-----------2        origin at cube center
!
!     Face numbering follows:
!       Face 1 = nodes 1,2,3,4
!       Face 2 = nodes 5,8,7,6
!       Face 3 = nodes 1,5,6,2
!       Face 4 = nodes 2,6,7,3
!       Face 5 = nodes 3,7,8,4
!       Face 6 = nodes 4,8,5,1
!
! Shawn A. Chester, December 2010 -- as used in my prior publications
! Shawn A. Chester, December 2013 -- modified for public distribution
!
!***********************************************************************
!
! User element statement in the input file (set ? values as needed):
!
!  2D elements
!  *User Element,Nodes=4,Type=U?,Iproperties=2,Properties=9,Coordinates=2,Variables=?,Unsymm
!  1,2,11
!
!  3D elements
!  *User Element,Nodes=8,Type=U3,Iproperties=2,Properties=9,Coordinates=3,Variables=?,Unsymm
!  1,2,3,11
!
!
!     State Variables
!     --------------------------------------------------------------
!     Global SDV's (used for visualization)
!       1) polymer volume fraction (phi)
!
!     Local SDV's (used for the solution procedure)
!       j = 0
!       do k = 1,nInttPt
!          svars(1+j) = phi ---- polymer volume fraction at integ pt k
!          j = j + nlSdv
!       end loop over k
!
!     In the input file, set 'User output variables'= number of global SDV's
!
!     In the input file, set 'ngSdv'= number of global SDV's
!
!     In the input file, set 'nlSdv'= number of local SDV's
!
!     In the input file, set 'varibles'=(nlSdv*nInttPt)
!
!
!     Material Properties Vector
!     --------------------------------------------------------------
!     Kbulk  = props(1) ! Bulk modulus
!     C10    = props(2) ! Chi parameter
!     C01     = props(3) ! Coefficient of permeability
!     nlSdv  = jprops(1) ! Number of local sdv's per integ pt
!     ngSdv  = jprops(2) ! Number of global sdv's per integ pt
!
!***********************************************************************

subroutine UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, &
            PROPS, NPROPS, coords, MCRD, NNODE, Uall, DUall, Vel, &
            Accn, JTYPE, TIME, DTIME, KSTEP, KINC, JELEM, PARAMS, &
            NDLOAD, JDLTYP, ADLMAG, PREDEF, NPREDF, LFLAGS, MLVARX, &
            DDLMAG, MDLOAD, PNEWDT, JPROPS, NJPROP, PERIOD)
  use global
  implicit none

  ! Variables defined in UEL, passed back to Abaqus
  real(8), intent(out) :: RHS(MLVARX, 1), AMATRX(NDOFEL, NDOFEL), &
                    SVARS(NSVARS), ENERGY(8)
  
  ! Variables passed into UEL
  real(8), intent(in) :: PROPS(NPROPS), coords(MCRD, NNODE), &
                   Uall(NDOFEL), DUall(MLVARX, 1), Vel(NDOFEL), &
                   Accn(NDOFEL), TIME(2), DTIME, PARAMS(1), &
                   ADLMAG(MDLOAD, 1), PREDEF(2, NPREDF, NNODE), &
                   DDLMAG(MDLOAD, 1), PNEWDT, PERIOD
  integer, intent(in) :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, &
                   JTYPE, KSTEP, KINC, JELEM, NDLOAD, JDLTYP(MDLOAD, 1), &
                   NPREDF, LFLAGS(4), MLVARX, MDLOAD, JPROPS(NJPROP), &
                   NJPROP

  integer :: lenJobName, lenOutDir, nDim, nIntt, nInttS
  character(len=256) :: jobName, outDir, fileName

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nIntt = 8  ! number of volume integration points
  nInttS = 1 ! number of surface integration points
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !----------------------------------------------------------------
  ! Perform initial checks
  ! Open the debug/error message file
  ! call getJobName(jobName, lenJobName)
  ! call getOutDir(outDir, lenOutDir)
  ! fileName = trim(outDir(1:lenOutDir)) // '\aaMSGS_' // &
  !            trim(jobName(1:lenJobName)) // '.dat'
  ! open(unit=80, file=fileName, status='unknown')

  ! Check the procedure type, this should be a coupled displacement
  ! which are any of the following (64, 65, 72, 73)
  if ((lflags(1) == 1) .or. (lflags(1) == 65) .or. & !64,65
     (lflags(1) == 72) .or. (lflags(1) == 73)) then
   ! all is good
  else
   write(*,*) 'Abaqus does not have the right procedure'
   write(*,*) 'go back and check the procedure type'
   write(*,*) 'lflags(1)=', lflags(1)
   ! write(80,*) 'Abaqus does not have the right procedure'
   ! write(80,*) 'go back and check the procedure type'
   ! write(80,*) 'lflags(1)=', lflags(1)
   call exit
  endif

  ! Make sure Abaqus knows you are doing a large deformation problem
  ! I think this only matters when it comes to output in viewer
  if (lflags(2) == 0) then
   ! lflags(2)=0 -> small disp.
   ! lflags(2)=1 -> large disp.
   write(*,*) 'Abaqus thinks you are doing'
   write(*,*) 'a small displacement analysis'
   write(*,*) 'go in and set nlgeom=yes'
   ! write(80,*) 'Abaqus thinks you are doing'
   ! write(80,*) 'a small displacement analysis'
   ! write(80,*) 'go in and set nlgeom=yes'
   call exit
  endif

  ! Check to see if you are doing a general step or a linear perturbation step
  if (lflags(4) == 1) then
   ! lflags(4)=0 -> general step
   ! lflags(4)=1 -> linear perturbation step
   write(*,*) 'Abaqus thinks you are doing'
   write(*,*) 'a linear perturbation step'
   ! write(80,*) 'Abaqus thinks you are doing'
   ! write(80,*) 'a linear perturbation step'
   call exit
  endif

  ! Do nothing if a ``dummy'' step
  if (dtime == 0.0) return

  ! Done with initial checks
  !----------------------------------------------------------------
  if (jtype == 3) then
   ! This is a 3D analysis
   nDim = 3
   call U3D8(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, &
           PROPS, NPROPS, coords, MCRD, NNODE, Uall, DUall, Vel, &
           Accn, JTYPE, TIME, DTIME, KSTEP, KINC, JELEM, PARAMS, &
           NDLOAD, JDLTYP, ADLMAG, PREDEF, NPREDF, LFLAGS, MLVARX, &
           DDLMAG, MDLOAD, PNEWDT, JPROPS, NJPROP, PERIOD, &
           nDim, nIntt, nInttS)
  else
   ! We have a problem...
   write(*,*) 'Element type not supported, jtype=', jtype
   ! write(80,*) 'Element type not supported, jtype=', jtype
   call exit
  endif

  ! Done with this element, RHS and AMATRX already returned
  ! as output from the specific element routine called
  !----------------------------------------------------------------
  return
end subroutine UEL






!****************************************************************************
!     Utility subroutines
!****************************************************************************

!***************************************************************************
subroutine matInv3Dd(A, A_inv, det_A, istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      use global
      implicit none
      !
      integer :: istat
      real(8) :: A(3,3), A_inv(3,3), det_A, det_A_inv

      istat = 1

      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) - &
                              A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) + &
                              A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))

      if (det_A <= 0.0_8) then
            write(*,*) 'WARNING: subroutine matInv3D:'
            write(*,*) 'WARNING: det of mat=', det_A
            istat = 0
            return
      end if

      det_A_inv = 1.0_8 / det_A

      A_inv(1,1) = det_A_inv * (A(2,2)*A(3,3) - A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv * (A(3,2)*A(1,3) - A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv * (A(1,2)*A(2,3) - A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv * (A(3,1)*A(2,3) - A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv * (A(1,1)*A(3,3) - A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv * (A(2,1)*A(1,3) - A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv * (A(2,1)*A(3,2) - A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv * (A(3,1)*A(1,2) - A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv * (A(1,1)*A(2,2) - A(2,1)*A(1,2))

      return
end subroutine matInv3Dd

!***************************************************************************
subroutine matInv2D(A, A_inv, det_A, istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      use global
      implicit none
      !
      integer :: istat
      real(8) :: A(2,2), A_inv(2,2), det_A, det_A_inv

      istat = 1

      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)

      if (det_A <= 0.0_8) then
            write(*,*) 'WARNING: subroutine matInv2D:'
            write(*,*) 'WARNING: det of mat=', det_A
            istat = 0
            return
      end if

      det_A_inv = 1.0_8 / det_A

      A_inv(1,1) =  det_A_inv * A(2,2)
      A_inv(1,2) = -det_A_inv * A(1,2)
      A_inv(2,1) = -det_A_inv * A(2,1)
      A_inv(2,2) =  det_A_inv * A(1,1)

      return
end subroutine matInv2D

!****************************************************************************
subroutine mdet(A, det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      use global
      implicit none
      !
      real(8) :: A(3,3), det

      det = A(1,1)*A(2,2)*A(3,3) + &
                        A(1,2)*A(2,3)*A(3,1) + &
                        A(1,3)*A(2,1)*A(3,2) - &
                        A(3,1)*A(2,2)*A(1,3) - &
                        A(3,2)*A(2,3)*A(1,1) - &
                        A(3,3)*A(2,1)*A(1,2)

      return
end subroutine mdet

!****************************************************************************
subroutine onem0(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      use global
      implicit none
      !
      integer :: i, j
      real(8) :: A(3,3)

      do i = 1, 3
            do j = 1, 3
                  if (i == j) then
                        A(i,j) = 1.0_8
                  else
                        A(i,j) = 0.0_8
                  end if
            end do
      end do

      return
end subroutine onem0
      subroutine UVARM(UVAR, DIRECT, T, TIME, DTIME, CMNAME, ORNAME, &
                               NUVARM, NOEL, NPT, LAYER, KSPT, KSTEP, KINC, NDI, NSHR, COORD, &
                               JMAC, JMATYP, MATLAYO, LACCFLA)

        ! This subroutine is used to transfer SDV's from the UEL
        ! onto the dummy mesh for viewing. Note that an offset of
        ! ElemOffset is used between the real mesh and the dummy mesh.
        ! If your model has more than ElemOffset UEL elements, then
        ! this will need to be modified.

        use global

        include 'aba_param.inc'

        character(len=80) :: CMNAME, ORNAME
        character(len=3) :: FLGRAY(15)
        real :: UVAR(NUVARM), DIRECT(3,3), T(3,3), TIME(2)
        integer :: ARRAY(15), JARRAY(15), JMAC(*), JMATYP(*), COORD(*)
        integer :: i1

        ! The dimensions of the variables FLGRAY, ARRAY and JARRAY
        ! must be set equal to or greater than 15.

        ! uvar(1) = globalSdv(noel-ElemOffset,npt,1)
        ! for example
        ! uvar(2) = globalSdv(noel-ElemOffset,npt,2)
        do i1 = 1, nsdv
            UVAR(i1) = globalSdv(noel-ElemOffset, npt, i1)
        end do
        ! uvar(3) = globalSdv(noel-ElemOffset,npt,3)
        ! uvar(4) = globalSdv(noel-ElemOffset,npt,4)

        return
      end subroutine UVARM
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
SUBROUTINE csisomatfic(cisomatfic,cmisomatfic,distgr,det,ndi)



!>    ISOTROPIC MATRIX: SPATIAL 'FICTICIOUS' ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: cisomatfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: cmisomatfic(ndi,ndi,ndi,ndi) !(ndi,ndi) in previous version
DOUBLE PRECISION, INTENT(IN OUT)         :: distgr(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det



call push4(cisomatfic,cmisomatfic,distgr,det,ndi)

RETURN
END SUBROUTINE csisomatfic
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
SUBROUTINE cmatisomatfic(cmisomatfic,cbar,cbari1,cbari2,  &
        diso,unit2,unit4,det,ndi)



!>    ISOTROPIC MATRIX: MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: cmisomatfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: cbar(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: cbari1
DOUBLE PRECISION, INTENT(IN OUT)         :: cbari2
DOUBLE PRECISION, INTENT(IN)             :: diso(5)
DOUBLE PRECISION, INTENT(IN)             :: unit2(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit4(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1
    

DOUBLE PRECISION :: dudi1,dudi2,d2ud2i1,d2ud2i2,d2udi1i2
DOUBLE PRECISION :: aux,aux1,aux2,aux3,aux4
DOUBLE PRECISION :: uij,ukl,cij,ckl

dudi1=diso(1)
dudi2=diso(2)
d2ud2i1=diso(3)
d2ud2i2=diso(4)
d2udi1i2=diso(5)

aux1=four*(d2ud2i1+two*cbari1*d2udi1i2+ dudi2+cbari1*cbari1*d2ud2i2)
aux2=-four*(d2udi1i2+cbari1*d2ud2i2)
aux3=four*d2ud2i2
aux4=-four*dudi2

DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        uij=unit2(i1,j1)
        ukl=unit2(k1,l1)
        cij=cbar(i1,j1)
        ckl=cbar(k1,l1)
        aux=aux1*uij*ukl+ aux2*(uij*ckl+cij*ukl)+aux3*cij*ckl+  &
            aux4*unit4(i1,j1,k1,l1)
        cmisomatfic(i1,j1,k1,l1)=aux * det**(-four/three)
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE cmatisomatfic
SUBROUTINE contraction22(aux,lt,rt,ndi)
!>       DOUBLE CONTRACTION BETWEEN 2nd ORDER AND 2ND ORDER  TENSOR
!>      INPUT:
!>       LT - RIGHT 2ND ORDER TENSOR
!>       RT - LEFT  2nd ODER TENSOR
!>      OUTPUT:
!>       aux - DOUBLE CONTRACTED TENSOR (scalar)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: lt(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rt(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: aux
INTEGER :: i1,j1


    aux=zero
    DO i1=1,ndi
      DO j1=1,ndi
        aux=aux+lt(i1,j1)*rt(j1,i1)
      END DO
    END DO
RETURN
END SUBROUTINE contraction22
SUBROUTINE contraction24(s,LT,rt,ndi)



!>       DOUBLE CONTRACTION BETWEEN 4TH ORDER AND 2ND ORDER  TENSOR
!>      INPUT:
!>       LT - RIGHT 2ND ORDER TENSOR
!>       RT - LEFT  4TH ODER TENSOR
!>      OUTPUT:
!>       S - DOUBLE CONTRACTED TENSOR (2ND ORDER)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: s(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: lt(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rt(ndi,ndi,ndi,ndi)



INTEGER :: i1,j1,k1,l1


DOUBLE PRECISION :: aux



DO k1=1,ndi
  DO l1=1,ndi
    aux=zero
    DO i1=1,ndi
      DO j1=1,ndi
        aux=aux+lt(k1,l1)*rt(i1,j1,k1,l1)
      END DO
    END DO
    s(k1,l1)=aux
  END DO
END DO
RETURN
END SUBROUTINE contraction24
SUBROUTINE contraction42(s,LT,rt,ndi)



!>       DOUBLE CONTRACTION BETWEEN 4TH ORDER AND 2ND ORDER  TENSOR
!>      INPUT:
!>       LT - left 4TH ORDER TENSOR
!>       RT - right  2ND ODER TENSOR
!>      OUTPUT:
!>       S - DOUBLE CONTRACTED TENSOR (2ND ORDER)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: s(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: LT(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rt(ndi,ndi)


INTEGER :: i1,j1,k1,l1


DOUBLE PRECISION :: aux



DO i1=1,ndi
  DO j1=1,ndi
    aux=zero
    DO k1=1,ndi
      DO l1=1,ndi
        aux=aux+LT(i1,j1,k1,l1)*rt(k1,l1)
      END DO
    END DO
    s(i1,j1)=aux
  END DO
END DO
RETURN
END SUBROUTINE contraction42
SUBROUTINE contraction44(s,LT,rt,ndi)



!>       DOUBLE CONTRACTION BETWEEN 4TH ORDER TENSORS
!>      INPUT:
!>       LT - RIGHT 4TH ORDER TENSOR
!>       RT - LEFT  4TH ORDER TENSOR
!>      OUTPUT:
!>       S - DOUBLE CONTRACTED TENSOR (4TH ORDER)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: s(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: LT(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rt(ndi,ndi,ndi,ndi)



INTEGER :: i1,j1,k1,l1,m1,n1


DOUBLE PRECISION :: aux



DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        aux=zero
        DO m1=1,ndi
          DO n1=1,ndi
            aux=aux+LT(i1,j1,m1,n1)*rt(m1,n1,k1,l1)
          END DO
        END DO
        s(i1,j1,k1,l1)=aux
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE contraction44
SUBROUTINE deformation(f,c,b,ndi)



!>     RIGHT AND LEFT CAUCHY-GREEN DEFORMATION TENSORS
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: b(ndi,ndi)


!     RIGHT CAUCHY-GREEN DEFORMATION TENSOR
c=matmul(transpose(f),f)
!     LEFT CAUCHY-GREEN DEFORMATION TENSOR
b=matmul(f,transpose(f))
RETURN
END SUBROUTINE deformation
SUBROUTINE fslip(f,fbar,det,ndi)



!>     DISTORTION GRADIENT
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: fbar(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: det


INTEGER :: i1,j1

DOUBLE PRECISION :: scale1

!     JACOBIAN DETERMINANT/VOLUME RATIO (J = det(F))
det = f(1,1) * f(2,2) * f(3,3) - f(1,2) * f(2,1) * f(3,3)

IF (ndi == 3) THEN
  det = det + f(1,2) * f(2,3) * f(3,1) + f(1,3) * f(3,2) * f(2,1)  &
      - f(1,3) * f(3,1) * f(2,2) - f(2,3) * f(3,2) * f(1,1)
END IF

scale1=det**(-one /three)

DO i1=1,ndi
  DO j1=1,ndi
    fbar(i1,j1)=scale1*f(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE fslip
SUBROUTINE getoutdir(outdir, lenoutdir)



!>     GET CURRENT WORKING DIRECTORY
INCLUDE 'aba_param.inc'


CHARACTER (LEN=256), INTENT(IN OUT)      :: outdir
INTEGER, INTENT(OUT)                     :: lenoutdir



CALL getcwd(outdir)
!        OUTDIR=OUTDIR(1:SCAN(OUTDIR,'\',BACK=.TRUE.)-1)
lenoutdir=len_trim(outdir)

RETURN
END SUBROUTINE getoutdir
SUBROUTINE onem(a,aa,aas,ndi)



!>      THIS SUBROUTINE GIVES:
!>          2ND ORDER IDENTITY TENSORS - A
!>          4TH ORDER IDENTITY TENSOR - AA
!>          4TH ORDER SYMMETRIC IDENTITY TENSOR - AAS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: aa(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: aas(ndi,ndi,ndi,ndi)



INTEGER :: i,j,k,l

a = zero
aa = zero
aas = zero

DO i = 1, ndi
  a(i,i) = one
END DO

DO i=1,ndi
  DO j=1,ndi
    DO k=1,ndi
      DO l=1,ndi
        IF (i == k .and. j == l) then
          aa(i,j,k,l) = one
        END IF
        aas(i,j,k,l) = (one/two)*(a(i,k)*a(j,l)+a(i,l)*a(j,k))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE onem
SUBROUTINE indexx(stress,ddsdde,sig,tng,ntens,ndi)



!>    INDEXATION: FULL SIMMETRY  IN STRESSES AND ELASTICITY TENSORS
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
INTEGER, INTENT(IN)                      :: ntens
DOUBLE PRECISION, INTENT(OUT)            :: stress(ntens)
DOUBLE PRECISION, INTENT(OUT)            :: ddsdde(ntens,ntens)
DOUBLE PRECISION, INTENT(IN)             :: sig(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: tng(ndi,ndi,ndi,ndi)



INTEGER :: ii1(6),ii2(6), i1,j1


DOUBLE PRECISION :: pp1,pp2

ii1(1)=1
ii1(2)=2
ii1(3)=3
ii1(4)=1
ii1(5)=1
ii1(6)=2

ii2(1)=1
ii2(2)=2
ii2(3)=3
ii2(4)=2
ii2(5)=3
ii2(6)=3

DO i1=1,ntens
!       STRESS VECTOR
  stress(i1)=sig(ii1(i1),ii2(i1))
  DO j1=1,ntens
!       DDSDDE - FULLY SIMMETRY IMPOSED
    pp1=tng(ii1(i1),ii2(i1),ii1(j1),ii2(j1))
    pp2=tng(ii1(i1),ii2(i1),ii2(j1),ii1(j1))
    ddsdde(i1,j1)=(one/two)*(pp1+pp2)
  END DO
END DO

RETURN

END SUBROUTINE indexx
SUBROUTINE initialize(statev)
use global
IMPLICIT NONE

!      DOUBLE PRECISION TIME(2),KSTEP
INTEGER :: pos1, i
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)


pos1=0
!       DETERMINANT
statev(pos1+1)=one
!       CL RELATIVE STIFFNESS
DO i = pos1+2, nsdv
    statev(i)=zero
END DO
!        CONTRACTION VARIANCE
!statev(pos1+2)=zero

RETURN

END SUBROUTINE initialize
SUBROUTINE invariants(a,inv1,inv2,ndi)



!>    1ST AND 2ND INVARIANTS OF A TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: inv1
DOUBLE PRECISION, INTENT(OUT)            :: inv2



INTEGER :: i1
DOUBLE PRECISION :: aa(ndi,ndi)
DOUBLE PRECISION :: inv1aa

inv1=zero
inv1aa=zero
aa=matmul(a,a)
DO i1=1,ndi
  inv1=inv1+a(i1,i1)
  inv1aa=inv1aa+aa(i1,i1)
END DO
inv2=(one/two)*(inv1*inv1-inv1aa)

RETURN
END SUBROUTINE invariants
SUBROUTINE isomat(sseiso,diso,c10,cbari1)



!>     ISOTROPIC MATRIX : ISOCHORIC SEF AND DERIVATIVES
use global
IMPLICIT NONE


DOUBLE PRECISION, INTENT(OUT)            :: sseiso
DOUBLE PRECISION, INTENT(OUT)            :: diso(5)
DOUBLE PRECISION, INTENT(IN)             :: c10
!DOUBLE PRECISION, INTENT(IN)             :: c01
DOUBLE PRECISION, INTENT(IN OUT)         :: cbari1
!DOUBLE PRECISION, INTENT(IN OUT)         :: cbari2


SSEISO=C10*(CBARI1-THREE)

!FIRST DERIVATIVE OF SSEISO IN ORDER TO I1
DISO(1)=C10
!FIRST DERIVATIVE OF SSEISO IN ORDER TO I2
DISO(2)=ZERO
!SECOND DERIVATIVE OF SSEISO IN ORDER TO I1
DISO(3)=ZERO
!SECOND DERIVATIVE OF SSEISO IN ORDER TO I2
DISO(4)=ZERO
!SECOND DERIVATIVE OF SSEISO IN ORDER TO I1 AND I2
DISO(5)=ZERO

RETURN
END SUBROUTINE isomat
SUBROUTINE metiso(cmiso,cmfic,pl,pkiso,pkfic,c,unit2,det,ndi)



!>    ISOCHORIC MATERIAL ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cmiso(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: cmfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pl(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pkiso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pkfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit2(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: det



INTEGER :: i1,j1,k1,l1
DOUBLE PRECISION :: cisoaux(ndi,ndi,ndi,ndi), cisoaux1(ndi,ndi,ndi,ndi),  &
    plt(ndi,ndi,ndi,ndi),cinv(ndi,ndi), pll(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: trfic,xx,yy,zz, aux,aux1

CALL matinv3d(c,cinv,ndi)
cisoaux1=zero
cisoaux=zero
CALL contraction44(cisoaux1,pl,cmfic,ndi)
DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        plt(i1,j1,k1,l1)=pl(k1,l1,i1,j1)
      END DO
    END DO
  END DO
END DO

CALL contraction44(cisoaux,cisoaux1,plt,ndi)

trfic=zero
aux=det**(-two/three)
aux1=aux**two
CALL contraction22(trfic,aux*pkfic,c,ndi)

DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        xx=aux1*cisoaux(i1,j1,k1,l1)
        pll(i1,j1,k1,l1)=(one/two)*(cinv(i1,k1)*cinv(j1,l1)+  &
            cinv(i1,l1)*cinv(j1,k1))- (one/three)*cinv(i1,j1)*cinv(k1,l1)
        yy=trfic*pll(i1,j1,k1,l1)
        zz=pkiso(i1,j1)*cinv(k1,l1)+cinv(i1,j1)*pkiso(k1,l1)
        
        cmiso(i1,j1,k1,l1)=xx+(two/three)*yy-(two/three)*zz
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE metiso
SUBROUTINE metvol(cvol,c,pv,ppv,det,ndi)



!>    VOLUMETRIC MATERIAL ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cvol(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pv
DOUBLE PRECISION, INTENT(IN OUT)         :: ppv
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1
DOUBLE PRECISION :: cinv(ndi,ndi)


CALL matinv3d(c,cinv,ndi)

DO i1 = 1, ndi
  DO j1 = 1, ndi
    DO k1 = 1, ndi
      DO l1 = 1, ndi
        cvol(i1,j1,k1,l1)= det*ppv*cinv(i1,j1)*cinv(k1,l1)  &
            -det*pv*(cinv(i1,k1)*cinv(j1,l1) +cinv(i1,l1)*cinv(j1,k1))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE metvol
SUBROUTINE matinv3d(a,a_inv,ndi)
!>    INVERSE OF A 3X3 MATRIX
!     RETURN THE INVERSE OF A(3,3) - A_INV
use global

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN)             :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: a_inv(ndi,ndi)

DOUBLE PRECISION :: det_a,det_a_inv

det_a = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) -  &
    a(2,1)*(a(1,2)*a(3,3) - a(3,2)*a(1,3)) +  &
    a(3,1)*(a(1,2)*a(2,3) - a(2,2)*a(1,3))

IF (det_a <= 0.d0) THEN
  WRITE(*,*) 'WARNING: SUBROUTINE MATINV3D:'
  WRITE(*,*) 'WARNING: DET OF MAT=',det_a
  RETURN
END IF

det_a_inv = 1.d0/det_a

a_inv(1,1) = det_a_inv*(a(2,2)*a(3,3)-a(3,2)*a(2,3))
a_inv(1,2) = det_a_inv*(a(3,2)*a(1,3)-a(1,2)*a(3,3))
a_inv(1,3) = det_a_inv*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
a_inv(2,1) = det_a_inv*(a(3,1)*a(2,3)-a(2,1)*a(3,3))
a_inv(2,2) = det_a_inv*(a(1,1)*a(3,3)-a(3,1)*a(1,3))
a_inv(2,3) = det_a_inv*(a(2,1)*a(1,3)-a(1,1)*a(2,3))
a_inv(3,1) = det_a_inv*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
a_inv(3,2) = det_a_inv*(a(3,1)*a(1,2)-a(1,1)*a(3,2))
a_inv(3,3) = det_a_inv*(a(1,1)*a(2,2)-a(2,1)*a(1,2))

RETURN
END SUBROUTINE matinv3d
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
SUBROUTINE pk2iso(pkiso,pkfic,pl,det,ndi)



!>    ISOCHORIC PK2 STRESS TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pkiso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pkfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pl(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: det



INTEGER :: i1,j1

DOUBLE PRECISION :: scale2

CALL contraction42(pkiso,pl,pkfic,ndi)

scale2=det**(-two/three)
DO i1=1,ndi
  DO j1=1,ndi
    pkiso(i1,j1)=scale2*pkiso(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE pk2iso
SUBROUTINE pk2isomatfic(fic,diso,cbar,cbari1,unit2,ndi)



!>     ISOTROPIC MATRIX: 2PK 'FICTICIOUS' STRESS TENSOR
!      INPUT:
!       DISO - STRAIN-ENERGY DERIVATIVES
!       CBAR - DEVIATORIC LEFT CAUCHY-GREEN TENSOR
!       CBARI1,CBARI2 - CBAR INVARIANTS
!       UNIT2 - 2ND ORDER IDENTITY TENSOR
!      OUTPUT:
!       FIC - 2ND PIOLA KIRCHOOF 'FICTICIOUS' STRESS TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: fic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: diso(5)
DOUBLE PRECISION, INTENT(IN)             :: cbar(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: cbari1
DOUBLE PRECISION, INTENT(IN)             :: unit2(ndi,ndi)



INTEGER :: i1,j1

DOUBLE PRECISION :: dudi1,dudi2
DOUBLE PRECISION :: aux1,aux2

dudi1=diso(1)
dudi2=diso(2)

aux1=two*(dudi1+cbari1*dudi2)
aux2=-two*dudi2

DO i1=1,ndi
  DO j1=1,ndi
    fic(i1,j1)=aux1*unit2(i1,j1)+aux2*cbar(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE pk2isomatfic
SUBROUTINE pk2vol(pkvol,pv,c,ndi)
!SUBROUTINE pk2vol(pkvol,pv,c,ndi, det)

! Shouldn't det be included??

!>    VOLUMETRIC PK2 STRESS
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pkvol(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pv
DOUBLE PRECISION, INTENT(IN OUT)         :: c(ndi,ndi)
!DOUBLE PRECISION, INTENT(IN)             :: det

INTEGER :: i1,j1
DOUBLE PRECISION :: cinv(ndi,ndi)


CALL matinv3d(c,cinv,ndi)

DO i1=1,ndi
  DO j1=1,ndi
    !pkvol(i1,j1)=det*pv*cinv(i1,j1)
    pkvol(i1,j1)=pv*cinv(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE pk2vol
SUBROUTINE projeul(a,aa,pe,ndi)



!>    EULERIAN PROJECTION TENSOR
!      INPUTS:
!          IDENTITY TENSORS - A, AA
!      OUTPUTS:
!          4TH ORDER SYMMETRIC EULERIAN PROJECTION TENSOR - PE
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: aa(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: pe(ndi,ndi,ndi,ndi)



INTEGER :: i,j,k,l



DO i=1,ndi
  DO j=1,ndi
    DO k=1,ndi
      DO l=1,ndi
        pe(i,j,k,l)=aa(i,j,k,l)-(one/three)*(a(i,j)*a(k,l))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE projeul
SUBROUTINE projlag(c,aa,pl,ndi)



!>    LAGRANGIAN PROJECTION TENSOR
!      INPUTS:
!          IDENTITY TENSORS - A, AA
!          ISOCHORIC LEFT CAUCHY GREEN TENSOR - C
!          INVERSE OF C - CINV
!      OUTPUTS:
!          4TH ORDER SYMMETRIC LAGRANGIAN PROJECTION TENSOR - PL
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: aa(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: pl(ndi,ndi,ndi,ndi)



INTEGER :: i,j,k,l

DOUBLE PRECISION :: cinv(ndi,ndi)

CALL matinv3d(c,cinv,ndi)

DO i=1,ndi
  DO j=1,ndi
    DO k=1,ndi
      DO l=1,ndi
        pl(i,j,k,l)=aa(i,j,k,l)-(one/three)*(cinv(i,j)*c(k,l))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE projlag
SUBROUTINE pull2(pk,sig,finv,det,ndi)



!>       PULL-BACK TIMES DET OF A 2ND ORDER TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pk(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: sig(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: finv(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: det



INTEGER :: i1,j1,ii1,jj1


DOUBLE PRECISION :: aux


DO i1=1,ndi
  DO j1=1,ndi
    aux=zero
    DO ii1=1,ndi
      DO jj1=1,ndi
        aux=aux+det*finv(i1,ii1)*finv(j1,jj1)*sig(ii1,jj1)
      END DO
    END DO
    pk(i1,j1)=aux
  END DO
END DO

RETURN
END SUBROUTINE pull2
SUBROUTINE pull4(mat,spatial,finv,det,ndi)



!>        PULL-BACK TIMES DET OF 4TH ORDER TENSOR

use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: mat(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: spatial(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: finv(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: det



INTEGER :: i1,j1,k1,l1,ii1,jj1,kk1,ll1


DOUBLE PRECISION :: aux


DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        aux=zero
        DO ii1=1,ndi
          DO jj1=1,ndi
            DO kk1=1,ndi
              DO ll1=1,ndi
                aux=aux+det* finv(i1,ii1)*finv(j1,jj1)*  &
                    finv(k1,kk1)*finv(l1,ll1)*spatial(ii1,jj1,kk1,ll1)
              END DO
            END DO
          END DO
        END DO
        mat(i1,j1,k1,l1)=aux
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE pull4
SUBROUTINE push2(sig,pk,f,det,ndi)



!>        PIOLA TRANSFORMATION
!>      INPUT:
!>       PK - 2ND PIOLA KIRCHOOF STRESS TENSOR
!>       F - DEFORMATION GRADIENT
!>       DET - DEFORMATION DETERMINANT
!>      OUTPUT:
!>       SIG - CAUCHY STRESS TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sig(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pk(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det


INTEGER :: i1,j1,ii1,jj1


DOUBLE PRECISION :: aux

DO i1=1,ndi
  DO j1=1,ndi
    aux=zero
    DO ii1=1,ndi
      DO jj1=1,ndi
        aux=aux+(det**(-one))*f(i1,ii1)*f(j1,jj1)*pk(ii1,jj1)
      END DO
    END DO
    sig(i1,j1)=aux
  END DO
END DO

RETURN
END SUBROUTINE push2
SUBROUTINE push4(spatial,mat,f,det,ndi)



!>        PIOLA TRANSFORMATION
!>      INPUT:
!>       MAT - MATERIAL ELASTICITY TENSOR
!>       F - DEFORMATION GRADIENT
!>       DET - DEFORMATION DETERMINANT
!>      OUTPUT:
!>       SPATIAL - SPATIAL ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: spatial(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: mat(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det


INTEGER :: i1,j1,k1,l1,ii1,jj1,kk1,ll1


DOUBLE PRECISION :: aux


DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        aux=zero
        DO ii1=1,ndi
          DO jj1=1,ndi
            DO kk1=1,ndi
              DO ll1=1,ndi
                aux=aux+(det**(-one))* f(i1,ii1)*f(j1,jj1)*  &
                    f(k1,kk1)*f(l1,ll1)*mat(ii1,jj1,kk1,ll1)
              END DO
            END DO
          END DO
        END DO
        spatial(i1,j1,k1,l1)=aux
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE push4
subroutine resetdfgrd(dfgrd, ndi)
      use global
      implicit none

      integer, intent(in) :: ndi
      real(8), intent(out) :: dfgrd(ndi, ndi)

      dfgrd = 0.0
      dfgrd(1,1) = 1.0
      dfgrd(2,2) = 1.0
      dfgrd(3,3) = 1.0

end subroutine resetdfgrd
SUBROUTINE rotation(f,r,u,ndi)



!>    COMPUTES ROTATION TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: r(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: u(ndi,ndi)




DOUBLE PRECISION :: uinv(ndi,ndi)

CALL matinv3d(u,uinv,ndi)

r = matmul(f,uinv)
RETURN
END SUBROUTINE rotation
SUBROUTINE sdvread(statev)
use global
implicit none
!>    VISCOUS DISSIPATION: READ STATE VARS
DOUBLE PRECISION, INTENT(IN)             :: statev(nsdv)




RETURN

END SUBROUTINE sdvread
SUBROUTINE sdvwrite(det,statev)
! VISCOUS DISSIPATION: WRITE STATE VARS
    use global
    IMPLICIT NONE

    DOUBLE PRECISION STATEV(NSDV),DET
    !write your sdvs here. they should be allocated 
    !after the viscous terms (check hvwrite)
     STATEV(1)=DET
    RETURN

END SUBROUTINE sdvwrite
SUBROUTINE setiso(ciso,cfic,pe,siso,sfic,unit2,ndi)


use global
IMPLICIT NONE

!>    ISOCHORIC SPATIAL ELASTICITY TENSOR

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: ciso(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pe(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: siso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: unit2(ndi,ndi)



INTEGER :: i1,j1,k1,l1
DOUBLE PRECISION :: cisoaux(ndi,ndi,ndi,ndi), cisoaux1(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: trfic,xx,yy,zz

cisoaux1=zero
cisoaux=zero

CALL contraction44(cisoaux1,pe,cfic,ndi)
CALL contraction44(cisoaux,cisoaux1,pe,ndi)

trfic=zero
DO i1=1,ndi
  trfic=trfic+sfic(i1,i1)
END DO

DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        xx=cisoaux(i1,j1,k1,l1)
        yy=trfic*pe(i1,j1,k1,l1)
        zz=siso(i1,j1)*unit2(k1,l1)+unit2(i1,j1)*siso(k1,l1)
        
        ciso(i1,j1,k1,l1)=xx+(two/three)*yy-(two/three)*zz
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE setiso
SUBROUTINE setjr(cjr,sigma,unit2,ndi)


use global
IMPLICIT NONE
!>    JAUMAN RATE CONTRIBUTION FOR THE SPATIAL ELASTICITY TENSOR

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cjr(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: sigma(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit2(ndi,ndi)



INTEGER :: i1,j1,k1,l1


DO i1 = 1, ndi
  DO j1 = 1, ndi
    DO k1 = 1, ndi
      DO l1 = 1, ndi
        
        cjr(i1,j1,k1,l1)= (one/two)*(unit2(i1,k1)*sigma(j1,l1)  &
            +sigma(i1,k1)*unit2(j1,l1)+unit2(i1,l1)*sigma(j1,k1)  &
            +sigma(i1,l1)*unit2(j1,k1))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE setjr
SUBROUTINE setvol(cvol,pv,ppv,unit2,unit4s,ndi)



!>    VOLUMETRIC SPATIAL ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cvol(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pv
DOUBLE PRECISION, INTENT(IN OUT)         :: ppv
DOUBLE PRECISION, INTENT(IN OUT)         :: unit2(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit4s(ndi,ndi,ndi,ndi)


INTEGER :: i1,j1,k1,l1



DO i1 = 1, ndi
  DO j1 = 1, ndi
    DO k1 = 1, ndi
      DO l1 = 1, ndi
        cvol(i1,j1,k1,l1)= ppv*unit2(i1,j1)*unit2(k1,l1)  &
            -two*pv*unit4s(i1,j1,k1,l1)
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE setvol
SUBROUTINE sigiso(siso,sfic,pe,ndi)



!>    ISOCHORIC CAUCHY STRESS
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: siso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pe(ndi,ndi,ndi,ndi)


CALL contraction42(siso,pe,sfic,ndi)

RETURN
END SUBROUTINE sigiso
SUBROUTINE sigisomatfic(sfic,pkfic,f,det,ndi)



!>    ISOTROPIC MATRIX:  ISOCHORIC CAUCHY STRESS
use global
IMPLICIT NONE


INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pkfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det





CALL push2(sfic,pkfic,f,det,ndi)

RETURN
END SUBROUTINE sigisomatfic
SUBROUTINE sigvol(svol,pv,unit2,ndi)



!>    VOLUMETRIC CAUCHY STRESS

use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: svol(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pv
DOUBLE PRECISION, INTENT(IN)             :: unit2(ndi,ndi)



INTEGER :: i1,j1



DO i1=1,ndi
  DO j1=1,ndi
    svol(i1,j1)=pv*unit2(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE sigvol
SUBROUTINE spectral(a,d,v)



!>    EIGENVALUES AND EIGENVECTOR OF A 3X3 MATRIX
!     THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
!     A SYMMETRIC 3X3 MATRIX A.

!     THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
!     EIGENVALUES IN ASCENDING ORDER, AND A MATRIX V WHOSE
!     COLUMNS CONTAIN THE CORRESPONDING EIGENVECTORS.

use global

DOUBLE PRECISION, INTENT(IN OUT)         :: a(3,3)
DOUBLE PRECISION                         :: e(3,3)
DOUBLE PRECISION, INTENT(IN OUT)         :: d(3)
DOUBLE PRECISION, INTENT(IN OUT)         :: v(3,3)

INTEGER :: nrot
INTEGER :: np=3



e = a

CALL jacobi(e,3,np,d,v,nrot)
CALL eigsrt(d,v,3,np)

RETURN
END SUBROUTINE spectral

!***********************************************************************

SUBROUTINE jacobi(a,n,np,d,v,nrot)

! COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
!  MATRIX A, WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL
!  NP BY NP ARRAY.  ON OUTPUT, ELEMENTS OF A ABOVE THE DIAGONAL
!  ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
!  AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
!  VECTOR D RETURNS THE EIGENVALUES OF A IN ITS FIRST N ELEMENTS.
!  V IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
!  A WHOSE COLUMNS CONTAIN, UPON OUTPUT, THE NORMALIZED
!  EIGENVECTORS OF A.  NROT RETURNS THE NUMBER OF JACOBI ROTATION
!  WHICH WERE REQUIRED.

! THIS SUBROUTINE IS TAKEN FROM 'NUMERICAL RECIPES.'
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: np
DOUBLE PRECISION, INTENT(IN OUT)         :: a(np,np)
INTEGER, INTENT(IN)                      :: n
DOUBLE PRECISION, INTENT(OUT)            :: d(np)
DOUBLE PRECISION, INTENT(OUT)            :: v(np,np)
INTEGER, INTENT(OUT)                     :: nrot

INTEGER :: ip,iq,  i,j
INTEGER, PARAMETER :: nmax=100

DOUBLE PRECISION :: b(nmax),z(nmax), sm,tresh,g,t,h,theta,s,c,tau


! INITIALIZE V TO THE IDENTITY MATRIX
DO i=1,3
  v(i,i)=one
  DO j=1,3
    IF (i /= j)THEN
      v(i,j)=zero
    END IF
  END DO
END DO
! INITIALIZE B AND D TO THE DIAGONAL OF A, AND Z TO ZERO.
!  THE VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
!  IN EQUATION (11.1.14)

DO ip = 1,n
  b(ip) = a(ip,ip)
  d(ip) = b(ip)
  z(ip) = 0.d0
END DO


! BEGIN ITERATION

nrot = 0
DO i=1,50
  
!         SUM OFF-DIAGONAL ELEMENTS
  
  sm = 0.d0
  DO ip=1,n-1
    DO iq=ip+1,n
      sm = sm + DABS(a(ip,iq))
    END DO
  END DO
  
!          IF SM = 0., THEN RETURN.  THIS IS THE NORMAL RETURN,
!          WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE
!          UNDERFLOW.
  
  IF (sm == 0.d0) RETURN
  
!          IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
!           |A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE,
!           SEE EQUATION (11.1.25).  THEREAFTER TRESH = 0.
  
  IF (i < 4) THEN
    tresh = 0.2D0*sm/n**2
  ELSE
    tresh = 0.d0
  END IF
  
  DO ip=1,n-1
    DO iq=ip+1,n
      g = 100.d0*DABS(a(ip,iq))
      
!              AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE
!               OFF-DIAGONAL ELEMENT IS SMALL.
      
      IF ((i > 4).AND.(DABS(d(ip))+g == DABS(d(ip)))  &
            .AND.(DABS(d(iq))+g == DABS(d(iq)))) THEN
        a(ip,iq) = 0.d0
      ELSE IF (DABS(a(ip,iq)) > tresh) THEN
        h = d(iq) - d(ip)
        IF (DABS(h)+g == DABS(h)) THEN
          
!                  T = 1./(2.*THETA), EQUATION (11.1.10)
          
          t =a(ip,iq)/h
        ELSE
          theta = 0.5D0*h/a(ip,iq)
          t =1.d0/(DABS(theta)+DSQRT(1.d0+theta**2.d0))
          IF (theta < 0.d0) t = -t
        END IF
        c = 1.d0/DSQRT(1.d0 + t**2.d0)
        s = t*c
        tau = s/(1.d0 + c)
        h = t*a(ip,iq)
        z(ip) = z(ip) - h
        z(iq) = z(iq) + h
        d(ip) = d(ip) - h
        d(iq) = d(iq) + h
        a(ip,iq) = 0.d0
        
!               CASE OF ROTATIONS 1 <= J < P
        
        DO j=1,ip-1
          g = a(j,ip)
          h = a(j,iq)
          a(j,ip) = g - s*(h + g*tau)
          a(j,iq) = h + s*(g - h*tau)
        END DO
        
!                CASE OF ROTATIONS P < J < Q
        
        DO j=ip+1,iq-1
          g = a(ip,j)
          h = a(j,iq)
          a(ip,j) = g - s*(h + g*tau)
          a(j,iq) = h + s*(g - h*tau)
        END DO
        
!                 CASE OF ROTATIONS Q < J <= N
        
        DO j=iq+1,n
          g = a(ip,j)
          h = a(iq,j)
          a(ip,j) = g - s*(h + g*tau)
          a(iq,j) = h + s*(g - h*tau)
        END DO
        DO j = 1,n
          g = v(j,ip)
          h = v(j,iq)
          v(j,ip) = g - s*(h + g*tau)
          v(j,iq) = h + s*(g - h*tau)
        END DO
        nrot = nrot + 1
      END IF
    END DO
  END DO
  
!          UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z
  
  DO ip=1,n
    b(ip) = b(ip) + z(ip)
    d(ip) = b(ip)
    z(ip) = 0.d0
  END DO
END DO

! IF THE ALGORITHM HAS REACHED THIS STAGE, THEN THERE
!  ARE TOO MANY SWEEPS.  PRINT A DIAGNOSTIC AND CUT THE
!  TIME INCREMENT.

WRITE (*,'(/1X,A/)') '50 ITERATIONS IN JACOBI SHOULD NEVER HAPPEN'

RETURN
END SUBROUTINE jacobi

!**********************************************************************

SUBROUTINE eigsrt(d,v,n,np)

!     GIVEN THE EIGENVALUES D AND EIGENVECTORS V AS OUTPUT FROM
!     JACOBI, THIS SUBROUTINE SORTS THE EIGENVALUES INTO ASCENDING
!     ORDER AND REARRANGES THE COLMNS OF V ACCORDINGLY.

!     THE SUBROUTINE WAS TAKEN FROM 'NUMERICAL RECIPES.'
use global

DOUBLE PRECISION, INTENT(IN OUT)         :: d(np)
DOUBLE PRECISION, INTENT(IN OUT)         :: v(np,np)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: np


INTEGER :: i,j,k

DOUBLE PRECISION :: p

DO i=1,n-1
  k = i
  p = d(i)
  DO j=i+1,n
    IF (d(j) >= p) THEN
      k = j
      p = d(j)
    END IF
  END DO
  IF (k /= i) THEN
    d(k) = d(i)
    d(i) = p
    DO j=1,n
      p = v(j,i)
      v(j,i) = v(j,k)
      v(j,k) = p
    END DO
  END IF
END DO

RETURN
END SUBROUTINE eigsrt
SUBROUTINE stretch(c,b,u,v,ndi)



!>    STRETCH TENSORS

use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi

DOUBLE PRECISION, INTENT(IN OUT)         :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: b(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: u(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: v(ndi,ndi)




DOUBLE PRECISION :: eigval(ndi),omega(ndi),eigvec(ndi,ndi)

CALL spectral(c,omega,eigvec)

eigval(1) = DSQRT(omega(1))
eigval(2) = DSQRT(omega(2))
eigval(3) = DSQRT(omega(3))

u(1,1) = eigval(1)
u(2,2) = eigval(2)
u(3,3) = eigval(3)

u = matmul(matmul(eigvec,u),transpose(eigvec))

CALL spectral(b,omega,eigvec)

eigval(1) = DSQRT(omega(1))
eigval(2) = DSQRT(omega(2))
eigval(3) = DSQRT(omega(3))
!      write(*,*) eigvec(1,1),eigvec(2,1),eigvec(3,1)

v(1,1) = eigval(1)
v(2,2) = eigval(2)
v(3,3) = eigval(3)

v = matmul(matmul(eigvec,v),transpose(eigvec))
RETURN
END SUBROUTINE stretch
subroutine tensorprod2(a, b, c, ndi)
       use global
       implicit none
       integer, intent(in) :: ndi
       real(8), intent(in) :: a(ndi, ndi), b(ndi, ndi)
       real(8), intent(out) :: c(ndi, ndi, ndi, ndi)
       integer :: i, j, k, l

       do i = 1, ndi
              do j = 1, ndi
                     do k = 1, ndi
                            do l = 1, ndi
                                   c(i, j, k, l) = a(i, j) * b(k, l)
                            end do
                     end do
              end do
       end do

end subroutine tensorprod2
!********************************************************************
! Record of revisions:                                              |
!        Date        Programmer        Description of change        |
!        ====        ==========        =====================        |
!     15/11/2017    Joao Ferreira      cont mech general eqs        |
!     01/11/2018    Joao Ferreira      comments added               |
!--------------------------------------------------------------------
!     Description:
!     UMAT: IMPLEMENTATION OF THE CONSTITUTIVE EQUATIONS BASED UPON 
!     A STRAIN-ENERGY FUNCTION (SEF).
!     THIS CODE, AS IS, EXPECTS A SEF BASED ON THE INVARIANTS OF THE 
!     CAUCHY-GREEN TENSORS. A VISCOELASTIC COMPONENT IS ALSO 
!     INCLUDED IF NEEDED. 
!     YOU CAN CHOOSE TO COMPUTE AT THE MATERIAL FRAME AND THEN 
!     PUSHFORWARD OR  COMPUTE AND THE SPATIAL FRAME DIRECTLY.
!--------------------------------------------------------------------
!     IF YOU WANT TO ADAPT THE CODE ACCORDING TO YOUR SEF:
!    ISOMAT - DERIVATIVES OF THE SEF IN ORDER TO THE INVARIANTS
!    ADD OTHER CONTRIBUTIONS: STRESS AND TANGENT MATRIX
!-------------------------------------------------------------------- 
!      STATE VARIABLES: CHECK ROUTINES - INITIALIZE, WRITESDV, READSDV.
!--------------------------------------------------------------------              
!     UEXTERNALDB: READ FILAMENTS ORIENTATION AND PREFERED DIRECTION
!--------------------------------------------------------------------
!---------------------------------------------------------------------
      SUBROUTINE MATERIAL(SIGMA,STATEV,DDSIGDDE,DFGRD0,DFGRD1,DET, &
      TIME,DTIME,PREDEF,NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS, &
      PNEWDT,NOEL,NPT,KSTEP,KINC)
!----------------------------------------------------------------------
!--------------------------- DECLARATIONS -----------------------------
!----------------------------------------------------------------------
      use global
      IMPLICIT NONE

!     ADD COMMON BLOCKS HERE IF NEEDED (and in uexternal)
!      COMMON /KBLOCK/KBLOCK

      CHARACTER(len=8) :: CMNAME

      INTEGER :: NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT, &
                 LAYER, KSPT, KSTEP, KINC

      REAL(KIND=8) :: STRESS(NTENS), STATEV(NSTATEV), &
                      DDSDDE(NTENS,NTENS), DDSDDT(NTENS), DRPLDE(NTENS), &
                      STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1), DPRED(1), &
                      PROPS(NPROPS), COORDS(3), DROT(3,3), DFGRD0(3,3), DFGRD1(3,3), &
                      FIBORI(NELEM,4)

      REAL(KIND=8) :: SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP, &
                      DTEMP, PNEWDT, CELENT

      INTEGER :: NTERM

!     FLAGS
!      INTEGER :: FLAG1

!     UTILITY TENSORS
      REAL(KIND=8) :: UNIT2(NDI,NDI), UNIT4(NDI,NDI,NDI,NDI), &
                      UNIT4S(NDI,NDI,NDI,NDI), &
                      PROJE(NDI,NDI,NDI,NDI), PROJL(NDI,NDI,NDI,NDI)

!     KINEMATICS
      REAL(KIND=8) :: DISTGR(NDI,NDI), C(NDI,NDI), B(NDI,NDI), &
                      CBAR(NDI,NDI), BBAR(NDI,NDI), DISTGRINV(NDI,NDI), &
                      UBAR(NDI,NDI), VBAR(NDI,NDI), ROT(NDI,NDI), &
                      DFGRD1INV(NDI,NDI)
      REAL(KIND=8) :: DET, CBARI1, CBARI2

!     VOLUMETRIC CONTRIBUTION
      REAL(KIND=8) :: PKVOL(NDI,NDI), SVOL(NDI,NDI), &
                      CVOL(NDI,NDI,NDI,NDI), CMVOL(NDI,NDI,NDI,NDI)
      REAL(KIND=8) :: KBULK, PV, PPV, SSEV

!     ISOCHORIC CONTRIBUTION
      REAL(KIND=8) :: SISO(NDI,NDI), PKISO(NDI,NDI), PK2(NDI,NDI), &
                      CISO(NDI,NDI,NDI,NDI), CMISO(NDI,NDI,NDI,NDI), &
                      SFIC(NDI,NDI), CFIC(NDI,NDI,NDI,NDI), &
                      PKFIC(NDI,NDI), CMFIC(NDI,NDI,NDI,NDI)

!     ISOCHORIC ISOTROPIC CONTRIBUTION
      REAL(KIND=8) :: C10, C01, SSEISO, DISO(5), PKMATFIC(NDI,NDI), &
                      SMATFIC(NDI,NDI), SISOMATFIC(NDI,NDI), &
                      CMISOMATFIC(NDI,NDI,NDI,NDI), &
                      CISOMATFIC(NDI,NDI,NDI,NDI)
      REAL(KIND=8) :: VORIF(3), VD(3), M0(3,3), MM(3,3), &
                      VORIF2(3), VD2(3), N0(3,3), NN(3,3)

!     LIST VARS OF OTHER CONTRIBUTIONS HERE

!     JAUMMAN RATE CONTRIBUTION (REQUIRED FOR ABAQUS UMAT)
      REAL(KIND=8) :: CJR(NDI,NDI,NDI,NDI)

!     CAUCHY STRESS AND ELASTICITY TENSOR
      REAL(KIND=8) :: SIGMA(NDI,NDI), DDSIGDDE(NDI,NDI,NDI,NDI), &
                      DDPKDDE(NDI,NDI,NDI,NDI)

!     TESTING/DEBUG VARS
      REAL(KIND=8) :: STEST(NDI,NDI), CTEST(NDI,NDI,NDI,NDI)
      INTEGER :: I1, J1, K1, L1
!----------------------------------------------------------------------
!-------------------------- INITIALIZATIONS ---------------------------
!----------------------------------------------------------------------
!     IDENTITY AND PROJECTION TENSORS
      UNIT2 = 0.0
      UNIT4 = 0.0
      UNIT4S = 0.0
      PROJE = 0.0
      PROJL = 0.0

!     KINEMATICS
      DISTGR = 0.0
      C = 0.0
      B = 0.0
      CBAR = 0.0
      BBAR = 0.0
      UBAR = 0.0
      VBAR = 0.0
      ROT = 0.0
      DET = 0.0
      CBARI1 = 0.0
      CBARI2 = 0.0

!     VOLUMETRIC
      PKVOL = 0.0
      SVOL = 0.0
      CVOL = 0.0
      KBULK = 0.0
      PV = 0.0
      PPV = 0.0
      SSEV = 0.0

!     ISOCHORIC
      SISO = 0.0
      PKISO = 0.0
      PK2 = 0.0
      CISO = 0.0
      CFIC = 0.0
      SFIC = 0.0
      PKFIC = 0.0

!     ISOTROPIC
      C10 = 0.0
      C01 = 0.0
      SSEISO = 0.0
      DISO = 0.0
      PKMATFIC = 0.0
      SMATFIC = 0.0
      SISOMATFIC = 0.0
      CMISOMATFIC = 0.0
      CISOMATFIC = 0.0

!     INITIALIZE OTHER CONT HERE

!     JAUMANN RATE
      CJR = 0.0

!     TOTAL CAUCHY STRESS AND ELASTICITY TENSORS
      SIGMA = 0.0
      DDSIGDDE = 0.0
!----------------------------------------------------------------------
!------------------------ IDENTITY TENSORS ----------------------------
!----------------------------------------------------------------------
      CALL ONEM(UNIT2, UNIT4, UNIT4S, NDI)
!----------------------------------------------------------------------
!------------------- MATERIAL CONSTANTS AND DATA ----------------------
!----------------------------------------------------------------------
!     VOLUMETRIC
      KBULK = PROPS(1)
!     ISOCHORIC ISOTROPIC NEO HOOKE
      C10 = PROPS(2)
!     NUMERICAL COMPUTATIONS
      NTERM = 60
!
!     STATE VARIABLES
!
      IF ((TIME(1) == 0.0) .AND. (KSTEP == 1)) THEN
         CALL INITIALIZE(STATEV)
      END IF
!        READ STATEV
      CALL SDVREAD(STATEV)
!      
!----------------------------------------------------------------------
!---------------------------- KINEMATICS ------------------------------
!----------------------------------------------------------------------
!     DISTORTION GRADIENT
      CALL FSLIP(DFGRD1, DISTGR, DET, NDI)
!     INVERSE OF DISTORTION GRADIENT
      CALL MATINV3D(DFGRD1, DFGRD1INV, NDI)
!     INVERSE OF DISTORTION GRADIENT
      CALL MATINV3D(DISTGR, DISTGRINV, NDI)
!     CAUCHY-GREEN DEFORMATION TENSORS
      CALL DEFORMATION(DFGRD1, C, B, NDI)
      CALL DEFORMATION(DISTGR, CBAR, BBAR, NDI)      
!     INVARIANTS OF DEVIATORIC DEFORMATION TENSORS
      CALL INVARIANTS(CBAR, CBARI1, CBARI2, NDI)
!     STRETCH TENSORS
      CALL STRETCH(CBAR, BBAR, UBAR, VBAR, NDI)
!     ROTATION TENSORS
      CALL ROTATION(DISTGR, ROT, UBAR, NDI)
!     DEVIATORIC PROJECTION TENSORS
      CALL PROJEUL(UNIT2, UNIT4S, PROJE, NDI)
!
      CALL PROJLAG(C, UNIT4, PROJL, NDI)
!----------------------------------------------------------------------
!--------------------- CONSTITUTIVE RELATIONS  ------------------------
!----------------------------------------------------------------------
!
!---- VOLUMETRIC ------------------------------------------------------
!     STRAIN-ENERGY AND DERIVATIVES (CHANGE HERE ACCORDING TO YOUR MODEL)
      CALL VOL(SSEV, PV, PPV, KBULK, DET)
      CALL ISOMAT(SSEISO, DISO, C10, CBARI1)
!
!---- ISOCHORIC ISOTROPIC ---------------------------------------------
!     PK2 'FICTICIOUS' STRESS TENSOR
      CALL PK2ISOMATFIC(PKMATFIC, DISO, CBAR, CBARI1, UNIT2, NDI)
!     CAUCHY 'FICTICIOUS' STRESS TENSOR
      CALL SIGISOMATFIC(SISOMATFIC, PKMATFIC, DISTGR, DET, NDI)
!     'FICTICIOUS' MATERIAL ELASTICITY TENSOR
      CALL CMATISOMATFIC(CMISOMATFIC, CBAR, CBARI1, CBARI2, DISO, UNIT2, UNIT4, DET, NDI)
!     'FICTICIOUS' SPATIAL ELASTICITY TENSOR
      CALL CSISOMATFIC(CISOMATFIC, CMISOMATFIC, DISTGR, DET, NDI)
!
!----------------------------------------------------------------------
!     SUM OF ALL ELASTIC CONTRIBUTIONS
!----------------------------------------------------------------------
!     STRAIN-ENERGY
      SSE = SSEV + SSEISO
!     PK2 'FICTICIOUS' STRESS
      PKFIC = PKMATFIC
!     CAUCHY 'FICTICIOUS' STRESS
      SFIC = SISOMATFIC
!     MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
      CMFIC = CMISOMATFIC
!     SPATIAL 'FICTICIOUS' ELASTICITY TENSOR
      CFIC = CISOMATFIC
!
!----------------------------------------------------------------------
!-------------------------- STRESS MEASURES ---------------------------
!----------------------------------------------------------------------
!
!---- VOLUMETRIC ------------------------------------------------------
!      PK2 STRESS
      CALL PK2VOL(PKVOL, PV, C, NDI)
!      CAUCHY STRESS
      CALL SIGVOL(SVOL, PV, UNIT2, NDI)
!
!---- ISOCHORIC -------------------------------------------------------
!      PK2 STRESS
      CALL PK2ISO(PKISO, PKFIC, PROJL, DET, NDI)
!      CAUCHY STRESS
      CALL SIGISO(SISO, SFIC, PROJE, NDI)
!
!---- VOLUMETRIC + ISOCHORIC ------------------------------------------
!      PK2 STRESS
      PK2 = PKVOL + PKISO
!      CAUCHY STRESS
      SIGMA = SVOL + SISO
!----------------------------------------------------------------------
!-------------------- MATERIAL ELASTICITY TENSOR ----------------------
!----------------------------------------------------------------------
!
!---- VOLUMETRIC ------------------------------------------------------
!
      CALL METVOL(CMVOL, C, PV, PPV, DET, NDI)
!
!---- ISOCHORIC -------------------------------------------------------
!
      CALL METISO(CMISO, CMFIC, PROJL, PKISO, PKFIC, C, UNIT2, DET, NDI)
!
!----------------------------------------------------------------------
!
      DDPKDDE = CMVOL + CMISO
!
!----------------------------------------------------------------------
!--------------------- SPATIAL ELASTICITY TENSOR ----------------------
!----------------------------------------------------------------------
!
!---- VOLUMETRIC ------------------------------------------------------
!
      CALL SETVOL(CVOL, PV, PPV, UNIT2, UNIT4S, NDI)
!
!---- ISOCHORIC -------------------------------------------------------
!
      CALL SETISO(CISO, CFIC, PROJE, SISO, SFIC, UNIT2, NDI)
!
!-----JAUMMAN RATE ----------------------------------------------------
!
      CALL SETJR(CJR, SIGMA, UNIT2, NDI)
!
!----------------------------------------------------------------------
!
!     ELASTICITY TENSOR
      DDSIGDDE = CVOL + CISO + CJR
!
!----------------------------------------------------------------------
!------------------------- INDEX ALLOCATION ---------------------------
!----------------------------------------------------------------------
!     VOIGT NOTATION  - FULLY SYMMETRY IMPOSED
      CALL INDEXX(STRESS, DDSDDE, SIGMA, DDSIGDDE, NTENS, NDI)
!
!----------------------------------------------------------------------
!--------------------------- STATE VARIABLES --------------------------
!----------------------------------------------------------------------
!     DO K1 = 1, NTENS
!      STATEV(1:27) = VISCOUS TENSORS
       CALL SDVWRITE(DET, STATEV)
!     END DO
!----------------------------------------------------------------------
      !write(*,*) 'F0'
      !write(*,*) DFGRD0
      !write(*,*) 'F1'
      !write(*,*) DFGRD1(1,1), DFGRD1(1,2), DFGRD1(1,3)
      !write(*,*) DFGRD1(2,1), DFGRD1(2,2), DFGRD1(2,3)
      !write(*,*) DFGRD1(3,1), DFGRD1(3,2), DFGRD1(3,3)
      !write(*,*) 'Cauchy stress'
      !write(*,*) STRESS(1), STRESS(2), STRESS(3), STRESS(4), STRESS(5), STRESS(6)
      RETURN
      END
!----------------------------------------------------------------------
!--------------------------- END OF UMAT ------------------------------
!----------------------------------------------------------------------
!
SUBROUTINE vol(ssev,pv,ppv,k,det)

! Code converted using TO_F90 by Alan Miller
! Date: 2020-12-12  Time: 12:08:12

!>     VOLUMETRIC CONTRIBUTION :STRAIN ENERGY FUNCTION AND DERIVATIVES
use global
implicit none


DOUBLE PRECISION :: g, aux
DOUBLE PRECISION, INTENT(OUT)            :: ssev
DOUBLE PRECISION, INTENT(OUT)            :: pv
DOUBLE PRECISION, INTENT(OUT)            :: ppv
DOUBLE PRECISION, INTENT(IN)             :: k
DOUBLE PRECISION, INTENT(IN)             :: det


g=(one/four)*(det*det-one-two*LOG(det))

ssev=k*g

pv=k*(one/two)*(det-one/det)
aux=k*(one/two)*(one+one/(det*det))
ppv=pv+det*aux

RETURN
END SUBROUTINE vol
