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
