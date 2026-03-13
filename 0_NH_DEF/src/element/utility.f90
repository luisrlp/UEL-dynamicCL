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
